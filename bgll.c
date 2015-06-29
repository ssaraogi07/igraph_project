#include <igraph.h>
#include "igraph_community.h"
#include "igraph_constructors.h"
#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_arpack.h"
#include "igraph_arpack_internal.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_interrupt_internal.h"
#include "igraph_components.h"
#include "igraph_dqueue.h"
#include "igraph_progress.h"
#include "igraph_stack.h"
#include "igraph_spmatrix.h"
#include "igraph_statusbar.h"
#include "igraph_types_internal.h"
#include "igraph_conversion.h"
#include "igraph_centrality.h"
#include "config.h"

#include <string.h>
#include <math.h>


/* Structure storing a community */
typedef struct {
  igraph_integer_t size;           /* Size of the community */
  igraph_real_t weight_inside;     /* Sum of edge weights inside community */
  igraph_real_t weight_all;        /* Sum of edge weights starting/ending
                                      in the community */
} igraph_i_multilevel_community;

/* Global community list structure */
typedef struct {
  long int communities_no, vertices_no;  /* Number of communities, number of vertices */
  igraph_real_t weight_sum;              /* Sum of edges weight in the whole graph */
  igraph_i_multilevel_community *item;   /* List of communities */
  igraph_vector_t *membership;           /* Community IDs */
  igraph_vector_t *weights;        /* Graph edge weights */
} igraph_i_multilevel_community_list;

/* Computes the modularity of a community partitioning */
igraph_real_t igraph_i_multilevel_community_modularity(
  const igraph_i_multilevel_community_list *communities) {
  igraph_real_t result = 0;
  long int i;
  igraph_real_t m = communities->weight_sum;
  
  for (i = 0; i < communities->vertices_no; i++) {
    if (communities->item[i].size > 0) {
      result += (communities->item[i].weight_inside - communities->item[i].weight_all*communities->item[i].weight_all/m)/m;
    } 
  }

  return result;
}

typedef struct {
  long int from;
  long int to;
  long int id;
} igraph_i_multilevel_link;

int igraph_i_multilevel_link_cmp(const void *a, const void *b)
{
  long int r = (((igraph_i_multilevel_link*)a)->from -
                ((igraph_i_multilevel_link*)b)->from);
  if (r != 0) return (int) r;
  
  return (int) (((igraph_i_multilevel_link*)a)->to -
		((igraph_i_multilevel_link*)b)->to);
}

/* removes multiple edges and returns new edge id's for each edge in |E|log|E| */
int igraph_i_multilevel_simplify_multiple(igraph_t *graph, igraph_vector_t *eids) {
  long int ecount = igraph_ecount(graph);
  long int i, l = -1, last_from = -1, last_to = -1;
  igraph_bool_t directed = igraph_is_directed(graph);
  igraph_integer_t from, to;
  igraph_vector_t edges;
  igraph_i_multilevel_link *links;

  /* Make sure there's enough space in eids to store the new edge IDs */
  IGRAPH_CHECK(igraph_vector_resize(eids, ecount));

  links = igraph_Calloc(ecount, igraph_i_multilevel_link);
  if (links == 0) {
    IGRAPH_ERROR("multi-level community structure detection failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, links);

  for (i = 0; i < ecount; i++) {
    igraph_edge(graph, (igraph_integer_t) i, &from, &to);
    links[i].from = from;
    links[i].to = to;
    links[i].id = i;
  }  

  qsort((void*)links, (size_t) ecount, sizeof(igraph_i_multilevel_link),
      igraph_i_multilevel_link_cmp);

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  for (i = 0; i < ecount; i++) {
    if (links[i].from == last_from && links[i].to == last_to) {
      VECTOR(*eids)[links[i].id] = l;
      continue;
    }

    last_from = links[i].from;
    last_to = links[i].to;

    igraph_vector_push_back(&edges, last_from);
    igraph_vector_push_back(&edges, last_to);

    l++;

    VECTOR(*eids)[links[i].id] = l;
  }

  free(links);
  IGRAPH_FINALLY_CLEAN(1);

  igraph_destroy(graph);
  IGRAPH_CHECK(igraph_create(graph, &edges, igraph_vcount(graph), directed));

  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

typedef struct {
  long int community;
  igraph_real_t weight;
} igraph_i_multilevel_community_link;

int igraph_i_multilevel_community_link_cmp(const void *a, const void *b)
{
  return (int) (((igraph_i_multilevel_community_link*)a)->community -
		((igraph_i_multilevel_community_link*)b)->community);
}

/**
 * Given a graph, a community structure and a vertex ID, this method
 * calculates:
 *
 * - edges: the list of edge IDs that are incident on the vertex
 * - weight_all: the total weight of these edges
 * - weight_inside: the total weight of edges that stay within the same
 *   community where the given vertex is right now, excluding loop edges
 * - weight_loop: the total weight of loop edges
 * - links_community and links_weight: together these two vectors list the
 *   communities incident on this vertex and the total weight of edges
 *   pointing to these communities
 */
int igraph_i_multilevel_community_links(const igraph_t *graph,
  const igraph_i_multilevel_community_list *communities,
  igraph_integer_t vertex, igraph_vector_t *edges,
  igraph_real_t *weight_all, igraph_real_t *weight_inside, igraph_real_t *weight_loop,
  igraph_vector_t *links_community, igraph_vector_t *links_weight) {
  
  long int i, n, last = -1, c = -1;
  igraph_real_t weight = 1;
  long int to, to_community;
  long int community = (long int) VECTOR(*(communities->membership))[(long int)vertex];
  igraph_i_multilevel_community_link *links;

  *weight_all = *weight_inside = *weight_loop = 0;

  igraph_vector_clear(links_community);
  igraph_vector_clear(links_weight);

  /* Get the list of incident edges */
  igraph_incident(graph, edges, vertex, IGRAPH_ALL);

  n = igraph_vector_size(edges);
  links = igraph_Calloc(n, igraph_i_multilevel_community_link);
  if (links == 0) {
    IGRAPH_ERROR("multi-level community structure detection failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, links);

  for (i = 0; i < n; i++) {
    long int eidx = (long int) VECTOR(*edges)[i];
    weight = VECTOR(*communities->weights)[eidx];

    to = IGRAPH_OTHER(graph, eidx, vertex);

    *weight_all += weight;
    if (to == vertex) {
      *weight_loop += weight;
      
      links[i].community = community;
      links[i].weight = 0;
      continue;
    }

    to_community = (long int)VECTOR(*(communities->membership))[to];
    if (community == to_community)
      *weight_inside += weight;

    /* debug("Link %ld (C: %ld) <-> %ld (C: %ld)\n", vertex, community, to, to_community); */
    
    links[i].community = to_community;
    links[i].weight = weight;
  }

  /* Sort links by community ID and merge the same */
  qsort((void*)links, (size_t) n, sizeof(igraph_i_multilevel_community_link),
      igraph_i_multilevel_community_link_cmp);
  for (i = 0; i < n; i++) {
    to_community = links[i].community;
    if (to_community != last) {
      igraph_vector_push_back(links_community, to_community);
      igraph_vector_push_back(links_weight, links[i].weight);
      last = to_community;
      c++;
    } else {    
      VECTOR(*links_weight)[c] += links[i].weight;
    }
  }

  igraph_free(links);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

igraph_real_t igraph_i_multilevel_community_modularity_gain(
  const igraph_i_multilevel_community_list *communities,
  igraph_integer_t community, igraph_integer_t vertex,
  igraph_real_t weight_all, igraph_real_t weight_inside) {
  IGRAPH_UNUSED(vertex);
  return weight_inside -
      communities->item[(long int)community].weight_all*weight_all/communities->weight_sum;
}

/* Shrinks communities into single vertices, keeping all the edges.
 * This method is internal because it destroys the graph in-place and
 * creates a new one -- this is fine for the multilevel community
 * detection where a copy of the original graph is used anyway.
 * The membership vector will also be rewritten by the underlying
 * igraph_membership_reindex call */
int igraph_i_multilevel_shrink(igraph_t *graph, igraph_vector_t *membership) {
  igraph_vector_t edges;
  long int no_of_nodes = igraph_vcount(graph);
  long int no_of_edges = igraph_ecount(graph);
  igraph_bool_t directed = igraph_is_directed(graph);

  long int i;
  igraph_eit_t eit;

  if (no_of_nodes == 0)
    return 0;

  if (igraph_vector_size(membership) < no_of_nodes) {
    IGRAPH_ERROR("cannot shrink graph, membership vector too short",
        IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);

  IGRAPH_CHECK(igraph_reindex_membership(membership, 0));

  /* Create the new edgelist */
  igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &eit);
  IGRAPH_FINALLY(igraph_eit_destroy, &eit);
  i = 0;
  while (!IGRAPH_EIT_END(eit)) {
    igraph_integer_t from, to;
    IGRAPH_CHECK(igraph_edge(graph, IGRAPH_EIT_GET(eit), &from, &to));
    VECTOR(edges)[i++] = VECTOR(*membership)[(long int) from];
    VECTOR(edges)[i++] = VECTOR(*membership)[(long int) to];
    IGRAPH_EIT_NEXT(eit);
  }
  igraph_eit_destroy(&eit);
  IGRAPH_FINALLY_CLEAN(1);

  /* Create the new graph */
  igraph_destroy(graph);
  no_of_nodes = (long int) igraph_vector_max(membership)+1;
  IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) no_of_nodes,
			     directed));

  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

/**
 * \ingroup communities
 * \function igraph_i_community_multilevel_step
 * \brief Performs a single step of the multi-level modularity optimization method
 *
 * This function implements a single step of the multi-level modularity optimization
 * algorithm for finding community structure, see VD Blondel, J-L Guillaume,
 * R Lambiotte and E Lefebvre: Fast unfolding of community hierarchies in large
 * networks, http://arxiv.org/abs/0803.0476 for the details.
 *
 * This function was contributed by Tom Gregorovic.
 *
 * \param graph   The input graph. It must be an undirected graph.
 * \param weights Numeric vector containing edge weights. If \c NULL, every edge
 *     has equal weight. The weights are expected to be non-negative.
 * \param membership The membership vector, the result is returned here.
 *     For each vertex it gives the ID of its community.
 * \param modularity The modularity of the partition is returned here.
 *     \c NULL means that the modularity is not needed.
 * \return Error code.
 *
 * Time complexity: in average near linear on sparse graphs.
 */
int igraph_i_community_multilevel_step(igraph_t *graph,
  igraph_vector_t *weights, igraph_vector_t *membership,
  igraph_real_t *modularity) {
  
  long int i, j;
  long int vcount = igraph_vcount(graph);
  long int ecount = igraph_ecount(graph);
  igraph_integer_t ffrom, fto;
  igraph_real_t q, pass_q;
  int pass;
  igraph_bool_t changed = 0;
  igraph_vector_t links_community;
  igraph_vector_t links_weight;
  igraph_vector_t edges;
  igraph_vector_t temp_membership;
  igraph_i_multilevel_community_list communities;

  /* Initial sanity checks on the input parameters */
  if (igraph_is_directed(graph)) {
    IGRAPH_ERROR("multi-level community detection works for undirected graphs only",
        IGRAPH_UNIMPLEMENTED);
  }
  if (igraph_vector_size(weights) < igraph_ecount(graph))
    IGRAPH_ERROR("multi-level community detection: weight vector too short", IGRAPH_EINVAL);
  if (igraph_vector_any_smaller(weights, 0))
    IGRAPH_ERROR("weights must be positive", IGRAPH_EINVAL);

  /* Initialize data structures */
  IGRAPH_VECTOR_INIT_FINALLY(&links_community, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&links_weight, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&temp_membership, vcount);
  IGRAPH_CHECK(igraph_vector_resize(membership, vcount));
 
  /* Initialize list of communities from graph vertices */
  communities.vertices_no = vcount;
  communities.communities_no = vcount;
  communities.weights = weights;
  communities.weight_sum = 2 * igraph_vector_sum(weights);
  communities.membership = membership;
  communities.item = igraph_Calloc(vcount, igraph_i_multilevel_community);
  if (communities.item == 0) {
    IGRAPH_ERROR("multi-level community structure detection failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, communities.item);

  /* Still initializing the communities data structure */
  for (i=0; i < vcount; i++) {
    VECTOR(*communities.membership)[i] = i;
    communities.item[i].size = 1;
    communities.item[i].weight_inside = 0;
    communities.item[i].weight_all = 0;
  }

  /* Some more initialization :) */
  for (i = 0; i < ecount; i++) {
    igraph_real_t weight = 1;
    igraph_edge(graph, (igraph_integer_t) i, &ffrom, &fto);

    weight = VECTOR(*weights)[i];
    communities.item[(long int) ffrom].weight_all += weight;
    communities.item[(long int) fto].weight_all += weight;
    if (ffrom == fto)
      communities.item[(long int) ffrom].weight_inside += 2*weight;
  }

  q = igraph_i_multilevel_community_modularity(&communities);
  pass = 1;

  do { /* Pass begin */
    long int temp_communities_no = communities.communities_no;

    pass_q = q;
    changed = 0;
    
    /* Save the current membership, it will be restored in case of worse result */
    IGRAPH_CHECK(igraph_vector_update(&temp_membership, communities.membership));

    for (i = 0; i < vcount; i++) {
      /* Exclude vertex from its current community */
      igraph_real_t weight_all = 0;
      igraph_real_t weight_inside = 0;
      igraph_real_t weight_loop = 0;
      igraph_real_t max_q_gain = 0;
      igraph_real_t max_weight;
      long int old_id, new_id, n;

      igraph_i_multilevel_community_links(graph, &communities, 
					  (igraph_integer_t) i, &edges,
					  &weight_all, &weight_inside, 
					  &weight_loop, &links_community,
					  &links_weight);

      old_id = (long int)VECTOR(*(communities.membership))[i];
      new_id = old_id;

      /* Update old community */
      igraph_vector_set(communities.membership, i, -1);
      communities.item[old_id].size--;
      if (communities.item[old_id].size == 0) {communities.communities_no--;}
      communities.item[old_id].weight_all -= weight_all;
      communities.item[old_id].weight_inside -= 2*weight_inside + weight_loop;

      /* debug("Remove %ld all: %lf Inside: %lf\n", i, -weight_all, -2*weight_inside + weight_loop); */

      /* Find new community to join with the best modification gain */
      max_q_gain = 0;
      max_weight = weight_inside;
      n = igraph_vector_size(&links_community);
      igraph_vector_shuffle(&links_community);
      for (j = 0; j < n; j++) {
        long int c = (long int) VECTOR(links_community)[j];
        igraph_real_t w = VECTOR(links_weight)[j];

        igraph_real_t q_gain = 
	  igraph_i_multilevel_community_modularity_gain(&communities, 
							(igraph_integer_t) c, 
							(igraph_integer_t) i,
							weight_all, w);
        /* debug("Link %ld -> %ld weight: %lf gain: %lf\n", i, c, (double) w, (double) q_gain); */
        if (q_gain > max_q_gain) {
          new_id = c;
          max_q_gain = q_gain;
          max_weight = w;
        }
      }

      /* debug("Added vertex %ld to community %ld (gain %lf).\n", i, new_id, (double) max_q_gain); */

      /* Add vertex to "new" community and update it */
      igraph_vector_set(communities.membership, i, new_id);
      if (communities.item[new_id].size == 0) {communities.communities_no++;}
      communities.item[new_id].size++;
      communities.item[new_id].weight_all += weight_all;
      communities.item[new_id].weight_inside += 2*max_weight + weight_loop;

      if (new_id != old_id) {
        changed++;
      }
    }

    q = igraph_i_multilevel_community_modularity(&communities);

    if (changed && (q > pass_q)) { 
      /* debug("Pass %d (changed: %d) Communities: %ld Modularity from %lf to %lf\n",
        pass, changed, communities.communities_no, (double) pass_q, (double) q); */
      pass++;
    } else {
      /* No changes or the modularity became worse, restore last membership */
      IGRAPH_CHECK(igraph_vector_update(communities.membership, &temp_membership));
      communities.communities_no = temp_communities_no;
      break;
    }

    IGRAPH_ALLOW_INTERRUPTION();
  } while (changed && (q > pass_q)); /* Pass end */

  if (modularity) {
    *modularity = q;
  }

  /* debug("Result Communities: %ld Modularity: %lf\n",
    communities.communities_no, (double) q); */

  IGRAPH_CHECK(igraph_reindex_membership(membership, 0));

  /* Shrink the nodes of the graph according to the present community structure
   * and simplify the resulting graph */

  /* TODO: check if we really need to copy temp_membership */
  IGRAPH_CHECK(igraph_vector_update(&temp_membership, membership));
  IGRAPH_CHECK(igraph_i_multilevel_shrink(graph, &temp_membership));
  igraph_vector_destroy(&temp_membership);
  IGRAPH_FINALLY_CLEAN(1);  
  
  /* Update edge weights after shrinking and simplification */
  /* Here we reuse the edges vector as we don't need the previous contents anymore */
  /* TODO: can we use igraph_simplify here? */
  IGRAPH_CHECK(igraph_i_multilevel_simplify_multiple(graph, &edges));

  /* We reuse the links_weight vector to store the old edge weights */
  IGRAPH_CHECK(igraph_vector_update(&links_weight, weights));
  igraph_vector_fill(weights, 0);
   
  for (i = 0; i < ecount; i++) {
    VECTOR(*weights)[(long int)VECTOR(edges)[i]] += VECTOR(links_weight)[i];
  }

  igraph_free(communities.item);
  igraph_vector_destroy(&links_community);
  igraph_vector_destroy(&links_weight);
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(4);
  
  return 0;
}

/**
 * \ingroup communities
 * \function igraph_community_multilevel
 * \brief Finding community structure by multi-level optimization of modularity
 * 
 * This function implements the multi-level modularity optimization
 * algorithm for finding community structure, see 
 * VD Blondel, J-L Guillaume, R Lambiotte and E Lefebvre: Fast unfolding of
 * community hierarchies in large networks, J Stat Mech P10008 (2008)
 * for the details (preprint: http://arxiv.org/abs/arXiv:0803.0476).
 *
 * It is based on the modularity measure and a hierarchical approach.
 * Initially, each vertex is assigned to a community on its own. In every step,
 * vertices are re-assigned to communities in a local, greedy way: each vertex
 * is moved to the community with which it achieves the highest contribution to
 * modularity. When no vertices can be reassigned, each community is considered
 * a vertex on its own, and the process starts again with the merged communities.
 * The process stops when there is only a single vertex left or when the modularity
 * cannot be increased any more in a step.
 *
 * This function was contributed by Tom Gregorovic.
 *
 * \param graph The input graph. It must be an undirected graph.
 * \param weights Numeric vector containing edge weights. If \c NULL, every edge
 *    has equal weight. The weights are expected to be non-negative.
 * \param membership The membership vector, the result is returned here.
 *    For each vertex it gives the ID of its community. The vector
 *    must be initialized and it will be resized accordingly.
 * \param memberships Numeric matrix that will contain the membership
 *     vector after each level, if not \c NULL. It must be initialized and
 *     it will be resized accordingly.
 * \param modularity Numeric vector that will contain the modularity score
 *     after each level, if not \c NULL. It must be initialized and it
 *     will be resized accordingly.
 * \return Error code.
 *
 * Time complexity: in average near linear on sparse graphs.
 * 
 * \example examples/simple/igraph_community_multilevel.c
 */

int igraph_community_multilevel(const igraph_t *graph,
  const igraph_vector_t *weights, igraph_vector_t *membership,
  igraph_matrix_t *memberships, igraph_vector_t *modularity) {
 
  igraph_t g;
  igraph_vector_t w, m, level_membership;
  igraph_real_t prev_q = -1, q = -1;
  int i, level = 1;
  long int vcount = igraph_vcount(graph);

  /* Make a copy of the original graph, we will do the merges on the copy */
  IGRAPH_CHECK(igraph_copy(&g, graph));
  IGRAPH_FINALLY(igraph_destroy, &g);

  if (weights) {
    IGRAPH_CHECK(igraph_vector_copy(&w, weights));   
    IGRAPH_FINALLY(igraph_vector_destroy, &w);  
  } else {
    IGRAPH_VECTOR_INIT_FINALLY(&w, igraph_ecount(&g));
    igraph_vector_fill(&w, 1);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&m, vcount);
  IGRAPH_VECTOR_INIT_FINALLY(&level_membership, vcount);

  if (memberships || membership) {
    /* Put each vertex in its own community */
    for (i = 0; i < vcount; i++) {
      VECTOR(level_membership)[i] = i;
    }
  }
  if (memberships) {
    /* Resize the membership matrix to have vcount columns and no rows */
    IGRAPH_CHECK(igraph_matrix_resize(memberships, 0, vcount));
  }
  if (modularity) {
    /* Clear the modularity vector */
    igraph_vector_clear(modularity);
  }
  
  while (1) {
    /* Remember the previous modularity and vertex count, do a single step */
    igraph_integer_t step_vcount = igraph_vcount(&g);

    prev_q = q;
    IGRAPH_CHECK(igraph_i_community_multilevel_step(&g, &w, &m, &q));

    /* Were there any merges? If not, we have to stop the process */
    if (igraph_vcount(&g) == step_vcount || q < prev_q)
      break;

    if (memberships || membership) {
      for (i = 0; i < vcount; i++) {
        /* Readjust the membership vector */
        VECTOR(level_membership)[i] = VECTOR(m)[(long int) VECTOR(level_membership)[i]];
      }
    }

    if (modularity) {
      /* If we have to return the modularity scores, add it to the modularity vector */
      IGRAPH_CHECK(igraph_vector_push_back(modularity, q));
    }

    if (memberships) {
      /* If we have to return the membership vectors at each level, store the new
       * membership vector */
      IGRAPH_CHECK(igraph_matrix_add_rows(memberships, 1));
      IGRAPH_CHECK(igraph_matrix_set_row(memberships, &level_membership, level - 1));
    }

    /* debug("Level: %d Communities: %ld Modularity: %f\n", level, (long int) igraph_vcount(&g),
      (double) q); */

    /* Increase the level counter */
    level++;
  }

  /* It might happen that there are no merges, so every vertex is in its 
     own community. We still might want the modularity score for that. */
  if (modularity && igraph_vector_size(modularity) == 0) {
    igraph_vector_t tmp;
    igraph_real_t mod;
    int i;
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, vcount);
    for (i=0; i<vcount; i++) { VECTOR(tmp)[i]=i; }
    IGRAPH_CHECK(igraph_modularity(graph, &tmp, &mod, weights));
    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_CHECK(igraph_vector_resize(modularity, 1));
    VECTOR(*modularity)[0]=mod;
  }

  /* If we need the final membership vector, copy it to the output */
  if (membership) {
    IGRAPH_CHECK(igraph_vector_resize(membership, vcount));   
    for (i = 0; i < vcount; i++) {
      VECTOR(*membership)[i] = VECTOR(level_membership)[i];
    }
  }

  /* Destroy the copy of the graph */
  igraph_destroy(&g);

  /* Destroy the temporary vectors */
  igraph_vector_destroy(&m);
  igraph_vector_destroy(&w);
  igraph_vector_destroy(&level_membership);
  IGRAPH_FINALLY_CLEAN(4);

  return 0;
}


void show_results(igraph_t *g, igraph_vector_t *membership, igraph_matrix_t *memberships, igraph_vector_t *modularity, FILE* f) {
  long int i, j, no_of_nodes = igraph_vcount(g);

  j=igraph_vector_which_max(modularity);
  for (i=0; i<igraph_vector_size(membership); i++) {
    if (VECTOR(*membership)[i] != MATRIX(*memberships, j, i)) {
      fprintf(f, "WARNING: best membership vector element %li does not match the best one in the membership matrix\n", i);
    }
  }

  fprintf(f, "Modularities:\n");
  igraph_vector_print(modularity);

  for (i=0; i < igraph_matrix_nrow(memberships); i++) {
    for (j=0; j < no_of_nodes; j++) {
      fprintf(f, "%ld ", (long int)MATRIX(*memberships, i, j));
    }
    fprintf(f, "\n");
  }

  fprintf(f, "\n");
}



int main() {
  igraph_t g;
  igraph_vector_t modularity, membership, edges;
  igraph_matrix_t memberships;
  int i, j, k;

  igraph_vector_init(&modularity,0);
  igraph_vector_init(&membership,0);
  igraph_matrix_init(&memberships,0,0);

  /* Unweighted test graph from the paper of Blondel et al */
  igraph_small(&g, 16, IGRAPH_UNDIRECTED,
      0, 2, 0, 3, 0, 4, 0, 5,
      1, 2, 1, 4, 1, 7,
      2, 4, 2, 5, 2, 6,
      3, 7,
      4, 10,
      5, 7, 5, 11,
      6, 7, 6, 11,
      8, 9, 8, 10, 8, 11, 8, 14, 8, 15,
      9, 12, 9, 14,
      10, 11, 10, 12, 10, 13, 10, 14,
      11, 13,
      -1);
  igraph_community_multilevel(&g, 0, &membership, &memberships, &modularity);
  show_results(&g, &membership, &memberships, &modularity, stdout);
  igraph_destroy(&g);

  /* Ring of 30 cliques */
 /* igraph_vector_init(&edges,0);
  for (i = 0; i < 30; i++) {
    for (j = 0; j < 5; j++) {
      for (k = j+1; k < 5; k++) {
        igraph_vector_push_back(&edges, i*5+j);
        igraph_vector_push_back(&edges, i*5+k);
      }
    }
  }
  for (i = 0; i < 30; i++) {
    igraph_vector_push_back(&edges, i*5 % 150);
    igraph_vector_push_back(&edges, (i*5+6) % 150);
  }
  igraph_create(&g, &edges, 150, 0);
  igraph_community_multilevel(&g, 0, &membership, &memberships, &modularity);
  show_results(&g, &membership, &memberships, &modularity, stdout);
  igraph_destroy(&g);*/

  igraph_vector_destroy(&modularity);
  igraph_vector_destroy(&membership);
  igraph_vector_destroy(&edges);
  igraph_matrix_destroy(&memberships);

#ifdef __APPLE__
  return 0;
#else
  return 77;
#endif
}
