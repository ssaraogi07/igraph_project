#include <igraph.h>
#include <igraph_types_internal.h>
#include "igraph_interrupt_internal.h"
#include <stdlib.h>

int igraph_get_shortest_paths_dijkstra(const igraph_t *graph,
                                       igraph_vector_ptr_t *vertices,
				       igraph_vector_ptr_t *edges,
				       igraph_integer_t from,
				       igraph_vs_t to,
				       const igraph_vector_t *weights,
				       igraph_neimode_t mode,
                                       igraph_vector_long_t *predecessors,
                                       igraph_vector_long_t *inbound_edges) {
  /* Implementation details. This is the basic Dijkstra algorithm, 
     with a binary heap. The heap is indexed, i.e. it stores not only
     the distances, but also which vertex they belong to. The other
     mapping, i.e. getting the distance for a vertex is not in the
     heap (that would by the double-indexed heap), but in the result
     matrix.

     Dirty tricks:
     - the opposite of the distance is stored in the heap, as it is a
       maximum heap and we need a minimum heap.
     - we don't use IGRAPH_INFINITY in the distance vector during the
       computation, as IGRAPH_FINITE() might involve a function call 
       and we want to spare that. So we store distance+1.0 instead of 
       distance, and zero denotes infinity.
     - `parents' assigns the inbound edge IDs of all vertices in the
       shortest path tree to the vertices. In this implementation, the
       edge ID + 1 is stored, zero means unreachable vertices.
  */
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_vit_t vit;
  igraph_2wheap_t Q;
  igraph_lazy_inclist_t inclist;
  igraph_vector_t dists;
  long int *parents;
  igraph_bool_t *is_target;
  long int i,to_reach;

  if (!weights) {
    return igraph_get_shortest_paths(graph, vertices, edges, from, to, mode,
        predecessors, inbound_edges);
  }
  
  if (igraph_vector_size(weights) != no_of_edges) {
    IGRAPH_ERROR("Weight vector length does not match", IGRAPH_EINVAL);
  }
  if (igraph_vector_min(weights) < 0) {
    IGRAPH_ERROR("Weight vector must be non-negative", IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_vit_create(graph, to, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);

  if (vertices && IGRAPH_VIT_SIZE(vit) != igraph_vector_ptr_size(vertices)) {
    IGRAPH_ERROR("Size of `vertices' and `to' should match", IGRAPH_EINVAL);
  }
  if (edges && IGRAPH_VIT_SIZE(vit) != igraph_vector_ptr_size(edges)) {
    IGRAPH_ERROR("Size of `edges' and `to' should match", IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_2wheap_init(&Q, no_of_nodes));
  IGRAPH_FINALLY(igraph_2wheap_destroy, &Q);
  IGRAPH_CHECK(igraph_lazy_inclist_init(graph, &inclist, mode));
  IGRAPH_FINALLY(igraph_lazy_inclist_destroy, &inclist);

  IGRAPH_VECTOR_INIT_FINALLY(&dists, no_of_nodes);
  igraph_vector_fill(&dists, -1.0);

  parents = igraph_Calloc(no_of_nodes, long int);
  if (parents == 0) IGRAPH_ERROR("Can't calculate shortest paths", IGRAPH_ENOMEM);
  IGRAPH_FINALLY(igraph_free, parents);
  is_target = igraph_Calloc(no_of_nodes, igraph_bool_t);
  if (is_target == 0) IGRAPH_ERROR("Can't calculate shortest paths", IGRAPH_ENOMEM);
  IGRAPH_FINALLY(igraph_free, is_target);

  /* Mark the vertices we need to reach */
  to_reach=IGRAPH_VIT_SIZE(vit);
  for (IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit)) {
    if (!is_target[ (long int) IGRAPH_VIT_GET(vit) ]) {
      is_target[ (long int) IGRAPH_VIT_GET(vit) ] = 1;
    } else {
      to_reach--;		/* this node was given multiple times */
    }
  }

  VECTOR(dists)[(long int)from] = 0.0;	/* zero distance */
  parents[(long int)from] = 0;
  igraph_2wheap_push_with_index(&Q, from, 0);
 
    int ar[10],a;
   srand(time(NULL));
    for(a=0;a<10;a++)
    ar[a]=rand()%100+1;
    int x=0;    

  while (!igraph_2wheap_empty(&Q) && to_reach > 0) {
    long int nlen, minnei=igraph_2wheap_max_index(&Q);
    igraph_real_t mindist=-igraph_2wheap_delete_max(&Q);
    igraph_vector_t *neis;
  

    IGRAPH_ALLOW_INTERRUPTION();

    if (is_target[minnei]) {
      is_target[minnei] = 0;
	  to_reach--;
	}

    /* Now check all neighbors of 'minnei' for a shorter path */
    neis=igraph_lazy_inclist_get(&inclist , (igraph_integer_t) minnei);
    nlen=igraph_vector_size(neis); 
   
    for (i=0; i<nlen; i++) {
      long int edge=(long int) VECTOR(*neis)[i];
      long int tto=IGRAPH_OTHER(graph, edge, minnei);
      igraph_real_t altdist=mindist + VECTOR(*weights)[edge];
      igraph_real_t curdist=VECTOR(dists)[tto];
      if (curdist < 0) {
        /* This is the first finite distance */
        VECTOR(dists)[tto] = altdist;
        parents[tto] = edge+1;
        IGRAPH_CHECK(igraph_2wheap_push_with_index(&Q, tto, -altdist));
      } else if (altdist < curdist) {
	      /* This is a shorter path */
        VECTOR(dists)[tto] = altdist;
        parents[tto] = edge+1;
        IGRAPH_CHECK(igraph_2wheap_modify(&Q, tto, -altdist));
      }
      else if(altdist==curdist)
      {
        
        int n=ar[x++];
        if (n>50) 
        parents[tto] = edge+1;
        if(x==10)
         x=0;
      }
    }
     
  } /* !igraph_2wheap_empty(&Q) */

  if (to_reach > 0) IGRAPH_WARNING("Couldn't reach some vertices");

  /* Create `predecessors' if needed */
  if (predecessors) {
    IGRAPH_CHECK(igraph_vector_long_resize(predecessors, no_of_nodes));

    for (i = 0; i < no_of_nodes; i++) {
      if (i == from) {
        /* i is the start vertex */
        VECTOR(*predecessors)[i] = i;
      } else if (parents[i] <= 0) {
        /* i was not reached */
        VECTOR(*predecessors)[i] = -1;
      } else {
        /* i was reached via the edge with ID = parents[i] - 1 */
        VECTOR(*predecessors)[i] = IGRAPH_OTHER(graph, parents[i]-1, i);
      }
    }
  }
  
  /* Create `inbound_edges' if needed */
  if (inbound_edges) {
    IGRAPH_CHECK(igraph_vector_long_resize(inbound_edges, no_of_nodes));

    for (i = 0; i < no_of_nodes; i++) {
      if (parents[i] <= 0) {
        /* i was not reached */
        VECTOR(*inbound_edges)[i] = -1;
      } else {
        /* i was reached via the edge with ID = parents[i] - 1 */
        VECTOR(*inbound_edges)[i] = parents[i]-1;
      }
    }
  }
  
  /* Reconstruct the shortest paths based on vertex and/or edge IDs */
  if (vertices || edges) {
    for (IGRAPH_VIT_RESET(vit), i=0; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit), i++) {
      long int node=IGRAPH_VIT_GET(vit);
      igraph_vector_t *vvec=0, *evec=0;
      if (vertices) {
	vvec=VECTOR(*vertices)[i];
	igraph_vector_clear(vvec);
      }
      if (edges) {
	evec=VECTOR(*edges)[i];
	igraph_vector_clear(evec);
      } 
  
      IGRAPH_ALLOW_INTERRUPTION();
      
      if (parents[node]>0) {
	long int size=0;
	long int act=node;
	long int edge;
	while (parents[act]) {
	  size++;
	  edge=parents[act]-1;
	  act=IGRAPH_OTHER(graph, edge, act);
	}
	if (vvec) { 
	  IGRAPH_CHECK(igraph_vector_resize(vvec, size+1)); 
	  VECTOR(*vvec)[size]=node;
	}
	if (evec) {
	  IGRAPH_CHECK(igraph_vector_resize(evec, size));
	}
	act=node;
	while (parents[act]) {
	  edge=parents[act]-1;
	  act=IGRAPH_OTHER(graph, edge, act);
	  size--;
	  if (vvec) { VECTOR(*vvec)[size]=act; }
	  if (evec) { VECTOR(*evec)[size]=edge; }
	}
      }
    }
  }
  
  igraph_lazy_inclist_destroy(&inclist);
  igraph_2wheap_destroy(&Q);
  igraph_vector_destroy(&dists);
  igraph_Free(is_target);
  igraph_Free(parents);
  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(6);
  
  return 0;
}

void print_vector(igraph_vector_t *v) {
  long int i, l=igraph_vector_size(v);
  for (i=0; i<l; i++) {
    printf(" %li", (long int) VECTOR(*v)[i]);
  }
  printf("\n");
}

int check_evecs(const igraph_t *graph, const igraph_vector_ptr_t *vecs,
		const igraph_vector_ptr_t *evecs, int error_code) {

  igraph_bool_t directed=igraph_is_directed(graph);
  long int i, n=igraph_vector_ptr_size(vecs);
  if (igraph_vector_ptr_size(evecs) != n) { exit(error_code+1); }
  
  for (i=0; i<n; i++) {
    igraph_vector_t *vvec=VECTOR(*vecs)[i];
    igraph_vector_t *evec=VECTOR(*evecs)[i];
    long int j, n2=igraph_vector_size(evec);
    if (igraph_vector_size(vvec) == 0 && n2==0) { continue; }
    if (igraph_vector_size(vvec) != n2+1) { exit(error_code+2); }
    for (j=0; j<n2; j++) {
      long int edge=VECTOR(*evec)[j];
      long int from=VECTOR(*vvec)[j];
      long int to=VECTOR(*vvec)[j+1];
      if (directed) {
	if (from != IGRAPH_FROM(graph, edge) ||
	    to   != IGRAPH_TO  (graph, edge)) {
	  exit(error_code);
	}
      } else {
	long int from2=IGRAPH_FROM(graph, edge);
	long int to2=IGRAPH_TO(graph, edge);
	long int min1= from < to ? from : to;
	long int max1= from < to ? to : from;
	long int min2= from2 < to2 ? from2 : to2;
	long int max2= from2 < to2 ? to2 : from2;
	if (min1 != min2 || max1 != max2) { exit(error_code+3); }
      }
    }
  }

  return 0;
}

int check_pred_inbound(const igraph_t* graph, const igraph_vector_long_t* pred,
        const igraph_vector_long_t* inbound, int start, int error_code) {
  long int i, n = igraph_vcount(graph);

  if (igraph_vector_long_size(pred) != n ||
      igraph_vector_long_size(inbound) != n) {
    exit(error_code);
  }

  if (VECTOR(*pred)[start] != start || VECTOR(*inbound)[start] != -1)
    exit(error_code+1);

  for (i = 0; i < n; i++) {
    if (VECTOR(*pred)[i] == -1) {
      if (VECTOR(*inbound)[i] != -1) {
        exit(error_code+2);
      }
    } else if (VECTOR(*pred)[i] == i) {
      if (i != start) {
        exit(error_code+3);
      }
      if (VECTOR(*inbound)[i] != -1) {
        exit(error_code+4);
      }
    } else {
      long int eid = VECTOR(*inbound)[i];
      long int u = IGRAPH_FROM(graph, eid), v = IGRAPH_TO(graph, eid);
      if (v != i && !igraph_is_directed(graph)) {
        long int dummy = u;
        u = v; v = dummy;
      }
      if (v != i) {
        exit(error_code+5);
      } else if (u != VECTOR(*pred)[i]) {
        exit(error_code+6);
      }
    }
  }

  return 0;
}

int main() {

  igraph_t g;
  igraph_vector_ptr_t vecs, evecs;
  igraph_vector_long_t pred, inbound;
  long int i;
  igraph_real_t weights[] = { 1, 2, 3, 4, 5, 1, 1, 1, 1, 1 }; 
  igraph_real_t weights2[] = { 0,2,1, 0,5,2, 1,1,0, 0,2,8, 1,1,3, 1,1,4, 2,1 };
  igraph_vector_t weights_vec;
  igraph_vs_t vs;

 
  igraph_ring(&g, 10, IGRAPH_UNDIRECTED, 0, 1);
  
  igraph_vector_ptr_init(&vecs, 5);
  igraph_vector_ptr_init(&evecs, 5);
  igraph_vector_long_init(&pred, 0);
  igraph_vector_long_init(&inbound, 0);

  for (i=0; i<igraph_vector_ptr_size(&vecs); i++) {
    VECTOR(vecs)[i] = calloc(1, sizeof(igraph_vector_t));
    igraph_vector_init(VECTOR(vecs)[i], 0);
    VECTOR(evecs)[i] = calloc(1, sizeof(igraph_vector_t));
    igraph_vector_init(VECTOR(evecs)[i], 0);
  }
  igraph_vs_vector_small(&vs, 1, 3,5,4, 6, -1);
  printf("The follwing examples show the shortest path for the vertices 1,3,5,4,6 respectively \n\n"); 
 

  printf(" A ring,  with weights \n");

  igraph_vector_view(&weights_vec, weights, sizeof(weights)/sizeof(igraph_real_t));
  igraph_get_shortest_paths_dijkstra(&g, /*vertices=*/ &vecs, 
				     /*edges=*/ &evecs, /*from=*/ 0, /*to=*/ vs, 
				     &weights_vec, IGRAPH_OUT,
				     /*predecessors=*/ &pred,
				     /*inbound_edges=*/ &inbound);
  
  check_evecs(&g, &vecs, &evecs, 20);
  check_pred_inbound(&g, &pred, &inbound, /* from= */ 0, 50);

  for (i=0; i<igraph_vector_ptr_size(&vecs); i++) {
    print_vector(VECTOR(vecs)[i]);
  }

  igraph_destroy(&g);

  printf("A complicated graph \n");

  igraph_small(&g, 10, IGRAPH_DIRECTED, 
	       0,1, 0,2, 0,3,    1,2, 1,4, 1,5,
	       2,3, 2,6,         3,2, 6,3,
	       4,5, 4,7,         5,6, 5,8, 5,9,
	       7,5, 7,8,         8,9,
	       5,2,
	       2,1,
	       -1);
  
  igraph_vector_view(&weights_vec, weights2, sizeof(weights2)/sizeof(igraph_real_t));
  igraph_get_shortest_paths_dijkstra(&g, /*vertices=*/ &vecs, 
				     /*edges=*/ &evecs, /*from=*/ 0, /*to=*/ vs, 
				     &weights_vec, IGRAPH_OUT,
				     /*predecessors=*/ &pred,
				     /*inbound_edges=*/ &inbound);

  check_evecs(&g, &vecs, &evecs, 30);
  check_pred_inbound(&g, &pred, &inbound, /* from= */ 0, 60);
  
  for (i=0; i<igraph_vector_ptr_size(&vecs); i++) {
    print_vector(VECTOR(vecs)[i]);
    igraph_vector_destroy(VECTOR(vecs)[i]);
    free(VECTOR(vecs)[i]);
    igraph_vector_destroy(VECTOR(evecs)[i]);
    free(VECTOR(evecs)[i]);
  }

  igraph_vector_ptr_destroy(&vecs);
  igraph_vector_ptr_destroy(&evecs);
  igraph_vector_long_destroy(&pred);
  igraph_vector_long_destroy(&inbound);
  
  igraph_vs_destroy(&vs);
  igraph_destroy(&g);

  if (!IGRAPH_FINALLY_STACK_EMPTY) return 1;

  return 0;
}


