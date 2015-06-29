#include <igraph.h>
#include <igraph_interrupt_internal.h>
int igraph_close_est( const igraph_t *graph, igraph_vector_t *res, 
		              const igraph_vs_t vids, igraph_neimode_t mode,
                              igraph_real_t cutoff,
			      const igraph_vector_t *weights,
			      igraph_bool_t normalized) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t already_counted;
  igraph_vector_int_t *neis;
  long int i, j;
  long int nodes_reached;
  igraph_adjlist_t allneis;

  igraph_dqueue_t q;
  
  long int nodes_to_calc;
  igraph_vit_t vit;

  if (weights) { 
    return igraph_i_closeness_estimate_weighted(graph, res, vids, mode, cutoff,
						weights, normalized);
  }

  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);

  nodes_to_calc=IGRAPH_VIT_SIZE(vit);
  
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("calculating closeness", IGRAPH_EINVMODE);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&already_counted, no_of_nodes);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);

  IGRAPH_CHECK(igraph_adjlist_init(graph, &allneis, mode));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &allneis);

  IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
  igraph_vector_null(res);
  
  for (IGRAPH_VIT_RESET(vit), i=0; 
       !IGRAPH_VIT_END(vit); 
       IGRAPH_VIT_NEXT(vit), i++) {
    igraph_dqueue_clear(&q);
    IGRAPH_CHECK(igraph_dqueue_push(&q, IGRAPH_VIT_GET(vit)));
    IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
    nodes_reached=1;
    VECTOR(already_counted)[(long int)IGRAPH_VIT_GET(vit)]=i+1;

    IGRAPH_PROGRESS("Closeness: ", 100.0*i/no_of_nodes, NULL);
    IGRAPH_ALLOW_INTERRUPTION();
    
    while (!igraph_dqueue_empty(&q)) {
      long int act=(long int) igraph_dqueue_pop(&q);
      long int actdist=(long int) igraph_dqueue_pop(&q);
      VECTOR(*res)[i] += actdist;

      if (cutoff>0 && actdist>=cutoff) continue;   /* NOT break!!! */

      neis=igraph_adjlist_get(&allneis, act);
      for (j=0; j<igraph_vector_int_size(neis); j++) {
        long int neighbor=(long int) VECTOR(*neis)[j];
        if (VECTOR(already_counted)[neighbor] == i+1) { continue; }
        VECTOR(already_counted)[neighbor] = i+1;
        nodes_reached++;
        IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
        IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
      }
    }
    /* using igraph_real_t here instead of igraph_integer_t to avoid overflow */
    printf("%li",nodes_reached);
    VECTOR(*res)[i] += ((igraph_real_t)no_of_nodes * (no_of_nodes-nodes_reached));
    VECTOR(*res)[i] = (no_of_nodes-1) / VECTOR(*res)[i];
  }

  if (!normalized) {
    for (i=0; i<nodes_to_calc; i++) {
      VECTOR(*res)[i] /= (no_of_nodes-1);
    }
  }

  IGRAPH_PROGRESS("Closeness: ", 100.0, NULL);

  /* Clean */
  igraph_dqueue_destroy(&q);
  igraph_vector_destroy(&already_counted);
  igraph_vit_destroy(&vit);
  igraph_adjlist_destroy(&allneis);
  IGRAPH_FINALLY_CLEAN(4);
  
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
  igraph_close_est(&g,&res,vs,IGRAPH_OUT, (igraph_real_t)number,
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
  igraph_close_est(&g,&res,vs,IGRAPH_OUT,(igraph_real_t)number,
				     NULL,1);
  igraph_vector_print(&res);
 
  igraph_vs_destroy(&vs);
  igraph_destroy(&g);
  igraph_vector_destroy(&res);

  return 0;
}


