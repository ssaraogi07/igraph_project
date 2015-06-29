#include <igraph.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "config.h"
#include "igraph_adjlist.h"
#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_dqueue.h"
#include "igraph_flow.h"
#include "igraph_interface.h"
#include "igraph_matching.h"
#include "igraph_structural.h"

/* #define MATCHING_DEBUG */

#ifdef _MSC_VER
/* MSVC does not support variadic macros */
#include <stdarg.h>
static void debug(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
#ifdef MATCHING_DEBUG
  vfprintf(stderr, fmt, args);
#endif
  va_end(args);
}
#else
#  ifdef MATCHING_DEBUG
#    define debug(...) fprintf(stderr, __VA_ARGS__)
#  else
#    define debug(...) 
#  endif
#endif

#define MATCHED(v) (VECTOR(match)[v] != -1)
#define UNMATCHED(v) (!MATCHED(v))

int igraph_i_maximum_bipartite_matching_weighted(const igraph_t* graph,
    const igraph_vector_bool_t* types, igraph_integer_t* matching_size,
    igraph_real_t* matching_weight, igraph_vector_long_t* matching,
    const igraph_vector_t* weights, igraph_real_t eps) {
  long int i, j, k, n, no_of_nodes, no_of_edges;
  igraph_integer_t u, v, w, msize;
  igraph_t newgraph;
  igraph_vector_long_t match;       /* will store the matching */
  igraph_vector_t slack;            /* will store the slack on each edge */
  igraph_vector_t parent;           /* parent vertices during a BFS */
  igraph_vector_t vec1, vec2;       /* general temporary vectors */
  igraph_vector_t labels;           /* will store the labels */
  igraph_dqueue_long_t q;           /* a FIFO for BST */
  igraph_bool_t smaller_set;        /* denotes which part of the bipartite graph is smaller */
  long int smaller_set_size;        /* size of the smaller set */
  igraph_real_t dual;               /* solution of the dual problem */
  igraph_adjlist_t tight_phantom_edges; /* adjacency list to manage tight phantom edges */
  igraph_integer_t alternating_path_endpoint;
  igraph_vector_t* neis;
  igraph_vector_int_t *neis2;
  igraph_inclist_t inclist;         /* incidence list of the original graph */ 

  /* The Hungarian algorithm is originally for complete bipartite graphs.
   * For non-complete bipartite graphs, a phantom edge of weight zero must be
   * added between every pair of non-connected vertices. We don't do this
   * explicitly of course. See the comments below about how phantom edges
   * are taken into account. */

  no_of_nodes = igraph_vcount(graph);
  no_of_edges = igraph_ecount(graph);
  if (eps < 0) {
    IGRAPH_WARNING("negative epsilon given, clamping to zero");
    eps = 0;
  }

  /* (1) Initialize data structures */
  IGRAPH_CHECK(igraph_vector_long_init(&match, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &match);
  IGRAPH_CHECK(igraph_vector_init(&slack, no_of_edges));
  IGRAPH_FINALLY(igraph_vector_destroy, &slack);
  IGRAPH_VECTOR_INIT_FINALLY(&vec1, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&vec2, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&labels, no_of_nodes);
  IGRAPH_CHECK(igraph_dqueue_long_init(&q, 0));
  IGRAPH_FINALLY(igraph_dqueue_long_destroy, &q);
  IGRAPH_VECTOR_INIT_FINALLY(&parent, no_of_nodes);
  IGRAPH_CHECK(igraph_adjlist_init_empty(&tight_phantom_edges, 
					 (igraph_integer_t) no_of_nodes));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &tight_phantom_edges);
  IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);

  /* (2) Find which set is the smaller one */
  j = 0;
  for (i = 0; i < no_of_nodes; i++) {
    if (VECTOR(*types)[i] == 0)
      j++;
  }
  smaller_set = (j > no_of_nodes / 2);
  smaller_set_size = smaller_set ? (no_of_nodes - j) : j;

  /* (3) Calculate the initial labeling and the set of tight edges. Use the
   *     smaller set only. Here we can assume that there are no phantom edges
   *     among the tight ones. */
  dual = 0;
  for (i = 0; i < no_of_nodes; i++) {
    igraph_real_t max_weight = 0;

    if (VECTOR(*types)[i] != smaller_set) {
      VECTOR(labels)[i] = 0;
      continue;
    }

    neis = igraph_inclist_get(&inclist, i);
    n = igraph_vector_size(neis);
    for (j = 0, k = 0; j < n; j++) {
      if (VECTOR(*weights)[(long int)VECTOR(*neis)[j]] > max_weight) {
        k = (long int) VECTOR(*neis)[j];
        max_weight = VECTOR(*weights)[k];
      }
    }

    VECTOR(labels)[i] = max_weight;
    dual += max_weight;
  }

  igraph_vector_clear(&vec1);
  IGRAPH_CHECK(igraph_get_edgelist(graph, &vec2, 0));
#define IS_TIGHT(i) (VECTOR(slack)[i] <= eps)
  for (i = 0, j = 0; i < no_of_edges; i++, j+=2) {
    u = (igraph_integer_t) VECTOR(vec2)[j];
    v = (igraph_integer_t) VECTOR(vec2)[j+1];
    VECTOR(slack)[i] = VECTOR(labels)[u] + VECTOR(labels)[v] - VECTOR(*weights)[i];
    if (IS_TIGHT(i)) {
      IGRAPH_CHECK(igraph_vector_push_back(&vec1, u));
      IGRAPH_CHECK(igraph_vector_push_back(&vec1, v));
    }
  }
  igraph_vector_clear(&vec2);

  /* (4) Construct a temporary graph on which the initial maximum matching
   *     will be calculated (only on the subset of tight edges) */
  IGRAPH_CHECK(igraph_create(&newgraph, &vec1,
			     (igraph_integer_t) no_of_nodes, 0));
  IGRAPH_FINALLY(igraph_destroy, &newgraph);
  IGRAPH_CHECK(igraph_maximum_bipartite_matching(&newgraph, types, &msize, 0, &match, 0, 0));
  igraph_destroy(&newgraph);
  IGRAPH_FINALLY_CLEAN(1);

  /* (5) Main loop until the matching becomes maximal */
  while (msize < smaller_set_size) {
    igraph_real_t min_slack, min_slack_2;
    igraph_integer_t min_slack_u, min_slack_v;

    /* (7) Fill the push queue with the unmatched nodes from the smaller set. */
    igraph_vector_clear(&vec1);
    igraph_vector_clear(&vec2);
    igraph_vector_fill(&parent, -1);
    for (i = 0; i < no_of_nodes; i++) {
      if (UNMATCHED(i) && VECTOR(*types)[i] == smaller_set) {
        IGRAPH_CHECK(igraph_dqueue_long_push(&q, i));
        VECTOR(parent)[i] = i;
        IGRAPH_CHECK(igraph_vector_push_back(&vec1, i));
      }
    }

#ifdef MATCHING_DEBUG
    debug("Matching:");
    igraph_vector_long_print(&match);
    debug("Unmatched vertices are marked by non-negative numbers:\n");
    igraph_vector_print(&parent);
    debug("Labeling:");
    igraph_vector_print(&labels);
    debug("Slacks:");
    igraph_vector_print(&slack);
#endif

    /* (8) Run the BFS */
    alternating_path_endpoint = -1;
    while (!igraph_dqueue_long_empty(&q)) {
      v = (int) igraph_dqueue_long_pop(&q);

      debug("Considering vertex %ld\n", (long int)v);

      /* v is always in the smaller set. Find the neighbors of v, which
       * are all in the larger set. Find the pairs of these nodes in
       * the smaller set and push them to the queue. Mark the traversed
       * nodes as seen.
       *
       * Here we have to be careful as there are two types of incident
       * edges on v: real edges and phantom ones. Real edges are
       * given by igraph_inclist_get. Phantom edges are not given so we
       * (ab)use an adjacency list data structure that lists the
       * vertices connected to v by phantom edges only. */
      neis = igraph_inclist_get(&inclist, v);
      n = igraph_vector_size(neis);
      for (i = 0; i < n; i++) {
        j = (long int) VECTOR(*neis)[i];
        /* We only care about tight edges */
        if (!IS_TIGHT(j))
          continue;
        /* Have we seen the other endpoint already? */
        u = IGRAPH_OTHER(graph, j, v);
        if (VECTOR(parent)[u] >= 0)
          continue;
        debug("  Reached vertex %ld via edge %ld\n", (long)u, (long)j);
        VECTOR(parent)[u] = v;
        IGRAPH_CHECK(igraph_vector_push_back(&vec2, u));
        w = (int) VECTOR(match)[u];
        if (w == -1) {
          /* u is unmatched and it is in the larger set. Therefore, we
           * could improve the matching by following the parents back
           * from u to the root.
           */
          alternating_path_endpoint = u;
          break;  /* since we don't need any more endpoints that come from v */
        } else {
          IGRAPH_CHECK(igraph_dqueue_long_push(&q, w));
          VECTOR(parent)[w] = u;
        }
        IGRAPH_CHECK(igraph_vector_push_back(&vec1, w));
      }

      /* Now do the same with the phantom edges */
      neis2 = igraph_adjlist_get(&tight_phantom_edges, v);
      n = igraph_vector_int_size(neis2);
      for (i = 0; i < n; i++) {
        u = (igraph_integer_t) VECTOR(*neis2)[i];
        /* Have we seen u already? */
        if (VECTOR(parent)[u] >= 0)
          continue;
        /* Check if the edge is really tight; it might have happened that the
         * edge became non-tight in the meanwhile. We do not remove these from
         * tight_phantom_edges at the moment, so we check them once again here.
         */
        if (fabs(VECTOR(labels)[(long int)v] + VECTOR(labels)[(long int)u]) > eps)
          continue;
        debug("  Reached vertex %ld via tight phantom edge\n", (long)u);
        VECTOR(parent)[u] = v;
        IGRAPH_CHECK(igraph_vector_push_back(&vec2, u));
        w = (int) VECTOR(match)[u];
        if (w == -1) {
          /* u is unmatched and it is in the larger set. Therefore, we
           * could improve the matching by following the parents back
           * from u to the root.
           */
          alternating_path_endpoint = u;
          break;  /* since we don't need any more endpoints that come from v */
        } else {
          IGRAPH_CHECK(igraph_dqueue_long_push(&q, w));
          VECTOR(parent)[w] = u;
        }
        IGRAPH_CHECK(igraph_vector_push_back(&vec1, w));
      }
    }

    /* Okay; did we have an alternating path? */
    if (alternating_path_endpoint != -1) {
#ifdef MATCHING_DEBUG
      debug("BFS parent tree:");
      igraph_vector_print(&parent);
#endif
      /* Increase the size of the matching with the alternating path. */
      v = alternating_path_endpoint;
      u = (igraph_integer_t) VECTOR(parent)[v];
      debug("Extending matching with alternating path ending in %ld.\n", (long int)v);

      while (u != v) {
        w = (int) VECTOR(match)[v];
        if (w != -1)
          VECTOR(match)[w] = -1;
        VECTOR(match)[v] = u;

        VECTOR(match)[v] = u;
        w = (int) VECTOR(match)[u];
        if (w != -1)
          VECTOR(match)[w] = -1;
        VECTOR(match)[u] = v;

        v = (igraph_integer_t) VECTOR(parent)[u];
	u = (igraph_integer_t) VECTOR(parent)[v];
      }

      msize++;

#ifdef MATCHING_DEBUG
      debug("New matching after update:");
      igraph_vector_long_print(&match);
      debug("Matching size is now: %ld\n", (long)msize);
#endif
      continue;
    }

#ifdef MATCHING_DEBUG
    debug("Vertices reachable from unmatched ones via tight edges:\n");
    igraph_vector_print(&vec1);
    igraph_vector_print(&vec2);
#endif

    /* At this point, vec1 contains the nodes in the smaller set (A)
     * reachable from unmatched nodes in A via tight edges only, while vec2
     * contains the nodes in the larger set (B) reachable from unmatched
     * nodes in A via tight edges only. Also, parent[i] >= 0 if node i
     * is reachable */

    /* Check the edges between reachable nodes in A and unreachable
     * nodes in B, and find the minimum slack on them.
     *
     * Since the weights are positive, we do no harm if we first
     * assume that there are no "real" edges between the two sets
     * mentioned above and determine an upper bound for min_slack
     * based on this. */
    min_slack = IGRAPH_INFINITY;
    min_slack_u = min_slack_v = 0;
    n = igraph_vector_size(&vec1);
    for (i = 0; i < no_of_nodes; i++) {
      if (VECTOR(*types)[i] == smaller_set)
        continue;
      if (VECTOR(labels)[i] < min_slack) {
        min_slack = VECTOR(labels)[i];
        min_slack_v = (igraph_integer_t) i;
      }
    }
    min_slack_2 = IGRAPH_INFINITY;
    for (i = 0; i < n; i++) {
      u = (igraph_integer_t) VECTOR(vec1)[i];
      /* u is surely from the smaller set, but we are interested in it
       * only if it is reachable from an unmatched vertex */
      if (VECTOR(parent)[u] < 0)
        continue;
      if (VECTOR(labels)[u] < min_slack_2) {
        min_slack_2 = VECTOR(labels)[u];
        min_slack_u = u;
      }
    }
    min_slack += min_slack_2;
    debug("Starting approximation for min_slack = %.4f (based on vertex pair %ld--%ld)\n",
        min_slack, (long int)min_slack_u, (long int)min_slack_v);

    n = igraph_vector_size(&vec1);
    for (i = 0; i < n; i++) {
      u = (igraph_integer_t) VECTOR(vec1)[i];
      /* u is a reachable node in A; get its incident edges.
       *
       * There are two types of incident edges: 1) real edges,
       * 2) phantom edges. Phantom edges were treated earlier
       * when we determined the initial value for min_slack. */
      debug("Trying to expand along vertex %ld\n", (long int)u);
      neis = igraph_inclist_get(&inclist, u);
      k = igraph_vector_size(neis);
      for (j = 0; j < k; j++) {
        /* v is the vertex sitting at the other end of an edge incident
         * on u; check whether it was reached */
        v = IGRAPH_OTHER(graph, VECTOR(*neis)[j], u);
        debug("  Edge %ld -- %ld (ID=%ld)\n", (long int)u, (long int)v, (long int)VECTOR(*neis)[j]);
        if (VECTOR(parent)[v] >= 0) {
          /* v was reached, so we are not interested in it */
          debug("    %ld was reached, so we are not interested in it\n", (long int)v);
          continue;
        }
        /* v is the ID of the edge from now on */
        v = (igraph_integer_t) VECTOR(*neis)[j];
        if (VECTOR(slack)[v] < min_slack) {
          min_slack = VECTOR(slack)[v];
          min_slack_u = u;
          min_slack_v = IGRAPH_OTHER(graph, v, u);
        }
        debug("    Slack of this edge: %.4f, min slack is now: %.4f\n",
            VECTOR(slack)[v], min_slack);
      }
    }
    debug("Minimum slack: %.4f on edge %d--%d\n", min_slack, (int)min_slack_u, (int)min_slack_v);

    if (min_slack > 0) {
      /* Decrease the label of reachable nodes in A by min_slack.
       * Also update the dual solution */
      n = igraph_vector_size(&vec1);
      for (i = 0; i < n; i++) {
        u = (igraph_integer_t) VECTOR(vec1)[i];
        VECTOR(labels)[u] -= min_slack;
        neis = igraph_inclist_get(&inclist, u);
        k = igraph_vector_size(neis);
        for (j = 0; j < k; j++) {
          debug("  Decreasing slack of edge %ld (%ld--%ld) by %.4f\n",
              (long)VECTOR(*neis)[j], (long)u,
              (long)IGRAPH_OTHER(graph, VECTOR(*neis)[j], u), min_slack);
          VECTOR(slack)[(long int)VECTOR(*neis)[j]] -= min_slack;
        }
        dual -= min_slack;
      }

      /* Increase the label of reachable nodes in B by min_slack.
       * Also update the dual solution */
      n = igraph_vector_size(&vec2);
      for (i = 0; i < n; i++) {
        u = (igraph_integer_t) VECTOR(vec2)[i];
        VECTOR(labels)[u] += min_slack;
        neis = igraph_inclist_get(&inclist, u);
        k = igraph_vector_size(neis);
        for (j = 0; j < k; j++) {
          debug("  Increasing slack of edge %ld (%ld--%ld) by %.4f\n",
              (long)VECTOR(*neis)[j], (long)u,
              (long)IGRAPH_OTHER(graph, (long)VECTOR(*neis)[j], u), min_slack);
          VECTOR(slack)[(long int)VECTOR(*neis)[j]] += min_slack;
        }
        dual += min_slack;
      }
    }

    /* Update the set of tight phantom edges.
     * Note that we must do it even if min_slack is zero; the reason is that
     * it can happen that min_slack is zero in the first step if there are
     * isolated nodes in the input graph.
     *
     * TODO: this is O(n^2) here. Can we do it faster? */
    for (u = 0; u < no_of_nodes; u++) {
      if (VECTOR(*types)[u] != smaller_set)
        continue;

      for (v = 0; v < no_of_nodes; v++) {
        if (VECTOR(*types)[v] == smaller_set)
          continue;

        if (VECTOR(labels)[(long int)u] + VECTOR(labels)[(long int)v] <= eps) {
          /* Tight phantom edge found. Note that we don't have to check whether
           * u and v are connected; if they were, then the slack of this edge
           * would be negative. */
          neis2 = igraph_adjlist_get(&tight_phantom_edges, u);
          if (!igraph_vector_int_binsearch(neis2, v, &i)) {
            debug("New tight phantom edge: %ld -- %ld\n", (long)u, (long)v);
            IGRAPH_CHECK(igraph_vector_int_insert(neis2, i, v));
          }
        }
      }
    }

#ifdef MATCHING_DEBUG
    debug("New labels:");
    igraph_vector_print(&labels);
    debug("Slacks after updating with min_slack:");
    igraph_vector_print(&slack);
#endif
  }

  /* Cleanup: remove phantom edges from the matching */
  for (i = 0; i < no_of_nodes; i++) {
    if (VECTOR(*types)[i] != smaller_set)
      continue;

    if (VECTOR(match)[i] != -1) {
      j = VECTOR(match)[i];
      neis2 = igraph_adjlist_get(&tight_phantom_edges, i);
      if (igraph_vector_int_binsearch(neis2, j, 0)) {
        VECTOR(match)[i] = VECTOR(match)[j] = -1;
        msize--;
      }
    }
  }

  /* Fill the output parameters */
  if (matching != 0) {
    IGRAPH_CHECK(igraph_vector_long_update(matching, &match));
  }
  if (matching_size != 0) {
    *matching_size = msize;
  }
  if (matching_weight != 0) {
    *matching_weight = 0;
    for (i = 0; i < no_of_edges; i++) {
      if (IS_TIGHT(i)) {
        IGRAPH_CHECK(igraph_edge(graph, (igraph_integer_t) i, &u, &v));
        if (VECTOR(match)[u] == v)
          *matching_weight += VECTOR(*weights)[i];
      }
    }
  }

  /* Release everything */
#undef IS_TIGHT
  igraph_inclist_destroy(&inclist);
  igraph_adjlist_destroy(&tight_phantom_edges);
  igraph_vector_destroy(&parent);
  igraph_dqueue_long_destroy(&q);
  igraph_vector_destroy(&labels);
  igraph_vector_destroy(&vec1);
  igraph_vector_destroy(&vec2);
  igraph_vector_destroy(&slack);
  igraph_vector_long_destroy(&match);
  IGRAPH_FINALLY_CLEAN(9);

  return IGRAPH_SUCCESS;
}

int test_weighted_graph_from_mit_notes() {
  /* Test graph from the following lecture notes:
   * http://math.mit.edu/~goemans/18433S07/matching-notes.pdf
   */
  igraph_t graph;
  igraph_vector_bool_t types;
  igraph_vector_long_t matching;
  igraph_vector_t weights;
  igraph_integer_t matching_size;
  igraph_real_t matching_weight;
  igraph_bool_t is_matching;
  igraph_real_t weight_array[] = { 2, 7, 2, 3,
      1, 3, 9, 3, 3,
      1, 3, 3, 1, 2,
      4, 1, 2,
      3 };
  int i;

  igraph_small(&graph, 0, 0,
      0, 6, 0, 7, 0, 8, 0, 9,
      1, 5, 1, 6, 1, 7, 1, 8, 1, 9,
      2, 5, 2, 6, 2, 7, 2, 8, 2, 9,
      3, 5, 3, 7, 3, 9,
      4, 7, -1);
  igraph_vector_bool_init(&types, 10);
  for (i = 0; i < 10; i++)
    VECTOR(types)[i] = (i >= 5);
  igraph_vector_long_init(&matching, 0);
  igraph_vector_init_copy(&weights, weight_array,
      sizeof(weight_array) / sizeof(weight_array[0]));

  igraph_i_maximum_bipartite_matching_weighted(&graph, &types, &matching_size,
      &matching_weight, &matching, &weights, 0);
  if (matching_size != 4) {
    printf("matching_size is %ld, expected: 4\n", (long)matching_size);
    return 1;
  }
  if (matching_weight != 19) {
    printf("matching_weight is %ld, expected: 19\n", (long)matching_weight);
    return 2;
  }
  igraph_is_maximal_matching(&graph, &types, &matching, &is_matching);
  if (!is_matching) {
    printf("not a matching: ");
    igraph_vector_long_print(&matching);
    return 3;
  }
  igraph_vector_long_print(&matching);
  igraph_vector_destroy(&weights);
  igraph_vector_long_destroy(&matching);
  igraph_vector_bool_destroy(&types);
  igraph_destroy(&graph);

  return 0;
}

int test_weighted_graph_generated() {
  /* Several randomly generated small test graphs */
  igraph_t graph;
  igraph_vector_bool_t types;
  igraph_vector_long_t matching;
  igraph_vector_t weights;
  igraph_integer_t matching_size;
  igraph_real_t matching_weight;
  igraph_real_t weight_array_1[] = { 8, 5, 9, 18, 20, 13 };
  igraph_real_t weight_array_2[] = { 20, 4, 20, 3, 13, 1 };
  int i;

  igraph_vector_bool_init(&types, 10);
  for (i = 0; i < 10; i++)
    VECTOR(types)[i] = (i >= 5);
  igraph_vector_long_init(&matching, 0);

  /* Case 1 */

  igraph_small(&graph, 0, 0, 0, 8, 2, 7, 3, 7, 3, 8, 4, 5, 4, 9, -1);
  igraph_vector_init_copy(&weights, weight_array_1,
      sizeof(weight_array_1) / sizeof(weight_array_1[0]));
  igraph_i_maximum_bipartite_matching_weighted(&graph, &types, &matching_size,
      &matching_weight, &matching, &weights, 0);
  if (matching_weight != 43) {
    printf("matching_weight is %ld, expected: 43\n", (long)matching_weight);
    return 2;
  }
  igraph_vector_long_print(&matching);
  igraph_vector_destroy(&weights);
  igraph_destroy(&graph);

  /* Case 2 */

  igraph_small(&graph, 0, 0, 0, 5, 0, 6, 1, 7, 2, 5, 3, 5, 3, 9, -1);
  igraph_vector_init_copy(&weights, weight_array_2,
      sizeof(weight_array_2) / sizeof(weight_array_2[0]));
  igraph_i_maximum_bipartite_matching_weighted(&graph, &types, &matching_size,
      &matching_weight, &matching, &weights, 0);
  if (matching_weight != 41) {
    printf("matching_weight is %ld, expected: 41\n", (long)matching_weight);
    return 2;
  }
  igraph_vector_destroy(&weights);
  igraph_destroy(&graph);
  igraph_vector_long_print(&matching);
  igraph_vector_long_destroy(&matching);
  igraph_vector_bool_destroy(&types);

  return 0;
}

int main() {
  igraph_i_set_attribute_table(&igraph_cattribute_table);


  if (test_weighted_graph_from_mit_notes())
    return 2;

  if (test_weighted_graph_generated())
    return 3;

  if (!IGRAPH_FINALLY_STACK_EMPTY) {
    printf("Finally stack still has %d elements.\n", IGRAPH_FINALLY_STACK_SIZE());
    return 4;
  }

  return 0;
}
