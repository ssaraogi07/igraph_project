#include "igraph_math.h"
#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "igraph_interrupt_internal.h"
#include "igraph_topology.h"
#include "igraph_types_internal.h"
#include "igraph_stack.h"
#include "igraph_dqueue.h"
#include "config.h"

#include "bigint.h"
#include "prpack.h"


void igraph_i_destroy_biguints(igraph_biguint_t *p) {
  igraph_biguint_t *p2 = p;
  while ( *((long int*)(p)) ) {
    igraph_biguint_destroy(p);
    p++;
  }
  igraph_Free(p2);
}
int igraph_bet_est(const igraph_t *graph, igraph_vector_t *res, 
				const igraph_vs_t vids, igraph_bool_t directed,
				igraph_real_t cutoff, 
				const igraph_vector_t *weights, 
				igraph_bool_t nobigint) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  long int *distance;
  unsigned long long int *nrgeo=0;  /* must be long long; consider grid
				       graphs for example */
  igraph_biguint_t *big_nrgeo=0;
  double *tmpscore;
  igraph_stack_t stack=IGRAPH_STACK_NULL;
  long int source;
  long int j, k, nneis;
  igraph_vector_int_t *neis;
  igraph_vector_t v_tmpres, *tmpres=&v_tmpres;
  igraph_vit_t vit;

  igraph_adjlist_t adjlist_out, adjlist_in;
  igraph_adjlist_t *adjlist_out_p, *adjlist_in_p;

  igraph_biguint_t D, R, T;

  if (weights) { 
    return igraph_i_betweenness_estimate_weighted(graph, res, vids, directed,
						cutoff, weights, nobigint);
  }

  if (!igraph_vs_is_all(&vids)) {
    /* subset */
    IGRAPH_VECTOR_INIT_FINALLY(tmpres, no_of_nodes);
  } else {
    /* only  */
    IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
    igraph_vector_null(res);
    tmpres=res;
  }

  directed=directed && igraph_is_directed(graph);
  if (directed) {
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist_out, IGRAPH_OUT));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist_out);
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist_in, IGRAPH_IN));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist_in);
    adjlist_out_p=&adjlist_out;
    adjlist_in_p=&adjlist_in;
  } else {
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist_out, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist_out);
    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist_in, IGRAPH_ALL));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist_in);
    adjlist_out_p=&adjlist_out;
    adjlist_in_p=&adjlist_in;
  }
  for (j=0; j<no_of_nodes; j++) {
    igraph_vector_int_clear(igraph_adjlist_get(adjlist_in_p, j));
  }
  
  distance=igraph_Calloc(no_of_nodes, long int);
  if (distance==0) {
    IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, distance);
  if (nobigint) {
    nrgeo=igraph_Calloc(no_of_nodes, unsigned long long int);
    if (nrgeo==0) {
      IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, nrgeo);
  } else {
    /* +1 is to have one containing zeros, when we free it, we stop
       at the zero */
    big_nrgeo=igraph_Calloc(no_of_nodes+1, igraph_biguint_t);
    if (!big_nrgeo) {
      IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_i_destroy_biguints, big_nrgeo);
    for (j=0; j<no_of_nodes; j++) {
      IGRAPH_CHECK(igraph_biguint_init(&big_nrgeo[j]));
    }
    IGRAPH_CHECK(igraph_biguint_init(&D));
    IGRAPH_FINALLY(igraph_biguint_destroy, &D);
    IGRAPH_CHECK(igraph_biguint_init(&R));
    IGRAPH_FINALLY(igraph_biguint_destroy, &R);
    IGRAPH_CHECK(igraph_biguint_init(&T));
    IGRAPH_FINALLY(igraph_biguint_destroy, &T);
  }
  tmpscore=igraph_Calloc(no_of_nodes, double);
  if (tmpscore==0) {
    IGRAPH_ERROR("betweenness failed", IGRAPH_ENOMEM);
  }
 IGRAPH_FINALLY(igraph_free, tmpscore);

  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  igraph_stack_init(&stack, no_of_nodes);
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);
    
  /* here we go */
  
  for (source=0; source<no_of_nodes; source++) {
    IGRAPH_PROGRESS("Betweenness centrality: ", 100.0*source/no_of_nodes, 0);
    IGRAPH_ALLOW_INTERRUPTION();

    IGRAPH_CHECK(igraph_dqueue_push(&q, source));
    if (nobigint) { 
      nrgeo[source]=1;
    } else {
      igraph_biguint_set_limb(&big_nrgeo[source], 1);
    }
    distance[source]=1;
    
    while (!igraph_dqueue_empty(&q)) {
      long int actnode=(long int) igraph_dqueue_pop(&q);
      IGRAPH_CHECK(igraph_stack_push(&stack, actnode));

      if (cutoff >= 0 && distance[actnode] >= cutoff+1) { continue; }
      
      neis = igraph_adjlist_get(adjlist_out_p, actnode);
      nneis = igraph_vector_int_size(neis);
      for (j=0; j<nneis; j++) {
        long int neighbor=(long int) VECTOR(*neis)[j];
        if (distance[neighbor]==0) {
	  distance[neighbor]=distance[actnode]+1;
	  IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	} 
	if (distance[neighbor]==distance[actnode]+1) {
	  igraph_vector_int_t *v=igraph_adjlist_get(adjlist_in_p, 
						    neighbor);
	  igraph_vector_int_push_back(v, actnode);
	  if (nobigint) { 
	    nrgeo[neighbor]+=nrgeo[actnode];
	  } else {
	    IGRAPH_CHECK(igraph_biguint_add(&big_nrgeo[neighbor],
					    &big_nrgeo[neighbor], 
					    &big_nrgeo[actnode]));
	  }
	}
      }
    } /* while !igraph_dqueue_empty */
    
    /* Ok, we've the distance of each node and also the number of
       shortest paths to them. Now we do an inverse search, starting
       with the farthest nodes. */
    while (!igraph_stack_empty(&stack)) {
      long int actnode=(long int) igraph_stack_pop(&stack);
      neis = igraph_adjlist_get(adjlist_in_p, actnode);
      nneis = igraph_vector_int_size(neis);
      for (j=0; j<nneis; j++) {
        long int neighbor=(long int) VECTOR(*neis)[j];
	if (nobigint) {
	  tmpscore[neighbor] +=  (tmpscore[actnode]+1)*
	    ((double)(nrgeo[neighbor]))/nrgeo[actnode];
	} else {
	  if (!igraph_biguint_compare_limb(&big_nrgeo[actnode], 0)) {
	    tmpscore[neighbor] = IGRAPH_INFINITY;
	  } else {
	    double div;
	    limb_t shift=1000000000L;
	    IGRAPH_CHECK(igraph_biguint_mul_limb(&T, &big_nrgeo[neighbor], 
						 shift));	  
	    igraph_biguint_div(&D, &R, &T, &big_nrgeo[actnode]);
	    div=igraph_biguint_get(&D) / shift;
	    tmpscore[neighbor] += (tmpscore[actnode]+1) * div;
	  }
	}
      }
      
      if (actnode != source) { VECTOR(*tmpres)[actnode] += tmpscore[actnode]; }

      distance[actnode]=0;
      if (nobigint) { 
	nrgeo[actnode]=0;
      } else {
	igraph_biguint_set_limb(&big_nrgeo[actnode], 0);
      }
      tmpscore[actnode]=0;
      igraph_vector_int_clear(igraph_adjlist_get(adjlist_in_p, actnode));
    }

  } /* for source < no_of_nodes */

  IGRAPH_PROGRESS("Betweenness centrality: ", 100.0, 0);

  /* clean  */
  igraph_Free(distance);
  if (nobigint) {
    igraph_Free(nrgeo); 
  } else {
    igraph_biguint_destroy(&T);
    igraph_biguint_destroy(&R);
    igraph_biguint_destroy(&D);
    IGRAPH_FINALLY_CLEAN(3);
    igraph_i_destroy_biguints(big_nrgeo);
  }
  igraph_Free(tmpscore);
  
  igraph_dqueue_destroy(&q);
  igraph_stack_destroy(&stack);
  IGRAPH_FINALLY_CLEAN(5);

  /* Keep only the requested vertices */
  if (!igraph_vs_is_all(&vids)) { 
    IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
    IGRAPH_FINALLY(igraph_vit_destroy, &vit);
    IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));

    for (k=0, IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit);
	 IGRAPH_VIT_NEXT(vit), k++) {
      long int node=IGRAPH_VIT_GET(vit);
      VECTOR(*res)[k] = VECTOR(*tmpres)[node];
    }

    igraph_vit_destroy(&vit);
    igraph_vector_destroy(tmpres);
    IGRAPH_FINALLY_CLEAN(2);
  }     

  /* divide by 2 for undirected graph */
  if (!directed) {
    nneis=igraph_vector_size(res);
    for (j=0; j<nneis; j++) {
      VECTOR(*res)[j] /= 2.0;
    }
  }
  
  igraph_adjlist_destroy(&adjlist_out);
  igraph_adjlist_destroy(&adjlist_in);
  IGRAPH_FINALLY_CLEAN(2);

  return 0;
}



int main() {

  igraph_t g;
   long int i;
   igraph_vs_t vs;
  igraph_vector_t res;
  igraph_vector_init(&res,0);
  long int number;
  srand(time(NULL));
  igraph_ring(&g, 10, IGRAPH_UNDIRECTED, 0, 1);
  number=rand()%10;
  igraph_vs_vector_small(&vs, 1, 3,5,4, 6, -1);
  number=0;
  igraph_bet_est(&g,&res,vs,IGRAPH_OUT, (igraph_real_t)number,
				     NULL,0);
  igraph_vector_print(&res);

  igraph_destroy(&g);
  number=rand()%10;
  igraph_small(&g, 10, IGRAPH_DIRECTED, 
	       0,1, 0,2, 0,3,    1,2, 1,4, 1,5,
	       2,3, 2,6,         3,2, 6,3,
	       4,5, 4,7,         5,6, 5,8, 5,9,
	       7,5, 7,8,         8,9,
	       5,2,
	       2,1,
	       -1);
  number = 0;
  igraph_bet_est(&g,&res,vs,IGRAPH_OUT,(igraph_real_t)number,
				     NULL,1);
  igraph_vector_print(&res);
 
  igraph_vs_destroy(&vs);
  igraph_destroy(&g);
  igraph_vector_destroy(&res);

  return 0;
}


