#include <igraph.h>
#include "igraph_structural.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_interrupt_internal.h"
#include "igraph_memory.h"
#include "igraph_progress.h"
#include "igraph_types_internal.h"

int igraph_i_minimum_spanning_tree_prim(const igraph_t* graph,
    igraph_vector_t* res, const igraph_vector_t *weights) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  char *already_added;
  char *al_ad;
  char *added_edges;

  igraph_d_indheap_t heap=IGRAPH_D_INDHEAP_NULL;
  igraph_integer_t mode=IGRAPH_ALL;
  igraph_integer_t ed,ef,et,fe;
  
  igraph_vector_t adj;
  igraph_vector_t temp;

  long int i, j, b, arr[50],edno;

  //igraph_vector_clear(res);
  igraph_vector_init(&temp,no_of_edges);
  if (weights == 0)
    return igraph_i_minimum_spanning_tree_unweighted(graph, res);

  if (igraph_vector_size(weights) != igraph_ecount(graph)) {
    IGRAPH_ERROR("Invalid weights length", IGRAPH_EINVAL);
  }

  added_edges=igraph_Calloc(no_of_edges, char);
  if (added_edges==0) {
    IGRAPH_ERROR("prim spanning tree failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, added_edges);
  already_added=igraph_Calloc(no_of_nodes, char);
  if (already_added == 0) {
    IGRAPH_ERROR("prim spanning tree failed", IGRAPH_ENOMEM);
  }
  al_ad=igraph_Calloc(no_of_nodes, char);
  if (al_ad == 0) {
    IGRAPH_ERROR("prim spanning tree failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, already_added);
  IGRAPH_FINALLY(igraph_free, al_ad);
  IGRAPH_CHECK(igraph_d_indheap_init(&heap, 0));
  IGRAPH_FINALLY(igraph_d_indheap_destroy, &heap);
  IGRAPH_VECTOR_INIT_FINALLY(&adj, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&temp, 0);
  for(b=0;b<no_of_edges;b++)
  {
	arr[b]=-1;
  }

  for (i=0; i<no_of_nodes; i++) {
    if (already_added[i]>0) { continue; }
    IGRAPH_ALLOW_INTERRUPTION();
    if (al_ad[i]>0) { continue; }
    IGRAPH_ALLOW_INTERRUPTION();
    already_added[i]=1;
    /* add all edges of the first vertex */
    igraph_incident(graph, &adj, (igraph_integer_t) i, (igraph_neimode_t) mode);
    for (j=0; j<igraph_vector_size(&adj); j++) {
      long int edgeno=(long int) VECTOR(adj)[j];
      igraph_integer_t edgefrom, edgeto;
      long int neighbor;
      igraph_edge(graph, (igraph_integer_t) edgeno, &edgefrom, &edgeto);
      neighbor= edgefrom != i ? edgefrom : edgeto;
      if (already_added[neighbor] == 0) {
	IGRAPH_CHECK(igraph_d_indheap_push(&heap, -VECTOR(*weights)[edgeno], i,
					   edgeno));
      }
    }
 }
  for(b=0;b<no_of_nodes;b++)
  {
	already_added[b]=0;
	al_ad[b]=-1;
  }
    while(! igraph_d_indheap_empty(&heap)) {
      /* Get minimal edge */
      long int from, edge,edge_w,a,x=0;
	//igraph_vector_init(&temp,no_of_edges);
      igraph_integer_t tmp, to;
      igraph_d_indheap_max_index(&heap, &from, &edge);
      igraph_edge(graph, (igraph_integer_t) edge, &tmp, &to);
	edge_w=VECTOR(*weights)[edge];
	//printf("Edge : %li Edge_w : %li ",edge,edge_w);
	x=0;

	while(1)
	{
		if(VECTOR(*weights)[edge]>edge_w)
		{ //printf(" now it breaks ");
		break;
		}
		else
		{
			if(igraph_d_indheap_empty(&heap)) break;
			igraph_d_indheap_max_index(&heap, &from, &edge);
			arr[edge]=from;		
			igraph_vector_insert(&temp,x,(igraph_real_t)edge);
			//igraph_vector_resize_min(&temp);
			//printf("x : %li ",x);
			//printf(" size : %li ",igraph_vector_size(&temp));
			x++;	
			if(igraph_d_indheap_empty(&heap)) break;
			igraph_d_indheap_delete_max(&heap);
			if(igraph_d_indheap_empty(&heap)) break;
			igraph_d_indheap_max_index(&heap, &from, &edge);
			
		}
	}
	igraph_vector_shuffle(&temp);
	printf("TEMP : ");
        igraph_vector_print(&temp);
      /* Erase it */
      //igraph_d_indheap_delete_max(&heap);

      /* Is this edge already included? */
	long int size=igraph_vector_size(&temp);
	long int nei;
     for(a=0;a<size;a++)
     {
	edge=VECTOR(temp)[a];
	from=arr[edge];
	igraph_edge(graph, (igraph_integer_t) edge, &tmp, &to);
	//printf("EDGE: %li TO: %li",edge,(long int) to);
      if (added_edges[edge]==0) 
      {
        if (from==to) { to=tmp;}
        /* Does it point to a visited node? */ 
	printf("EDGE: %li FROM: %li TO: %li ",edge,from,(long int)to);
        if (already_added[(long int)to]==0) 
	{
	  	if(al_ad[from]==1 && al_ad[(long int)to]==1)
	   	continue;
	  	//if(already_added[from]==1 && al_ad[(long int)to]==1)
	   	//continue;
          	already_added[(long int)to]=1;
          	added_edges[edge]=1;
		if(already_added[from]==0)
	  	al_ad[from]=1;
		//al_ad[(long int)to]=edge;
	  	IGRAPH_CHECK(igraph_vector_push_back(res, edge));
	 
         }
	/*else 
	{
		if(al_ad[from]==no_of_edges && al_ad[(long int)to]==no_of_edges)
	   	continue;
	  	if(already_added[from]==1 && al_ad[(long int)to]==no_of_edges)
	   	continue;
		if(already_added[from]==1)
		{
			igraph_edge(graph, (igraph_integer_t) al_ad[from], &ef, &et);
			//printf("for EDGE: %li ef: %li et: %li",&edge,(long int)ef,(long int)et);
			igraph_edge(graph, (igraph_integer_t) al_ad[(long int)to], &fe, &et);
			//printf("for EDGE: %li fe: %li et: %li",&edge,(long int)fe,(long int)et);
			if((fe==ef)||(already_added[(long int)ef]==1 && already_added[(long int)fe]==1))
			continue;
		}
          	already_added[(long int)to]=1;
          	added_edges[edge]=1;
		if(already_added[from]==0)
	  	al_ad[from]=no_of_edges;
		al_ad[(long int)to]=edge;
	  	IGRAPH_CHECK(igraph_vector_push_back(res, edge));
	}*/
	else if(already_added[from]==0 && al_ad[from]==0)
	{
		already_added[(long int)from]=1;
          	added_edges[edge]=1;
		if(already_added[from]==0)
	  	al_ad[from]=no_of_edges;
	  	IGRAPH_CHECK(igraph_vector_push_back(res, edge));
	 }
		
       } /* if !already_added */
      }
      igraph_vector_clear(&temp);
    } /* while in the same component */
   
   
   /* for all nodes */

  igraph_d_indheap_destroy(&heap);
  igraph_Free(already_added);
  igraph_Free(al_ad);
  igraph_vector_destroy(&adj);
  igraph_vector_destroy(&temp);
  igraph_Free(added_edges);
  IGRAPH_FINALLY_CLEAN(4);

  return IGRAPH_SUCCESS;
}

int main() {
  
  igraph_t g;
  igraph_vector_t  res;
  long int i;
  igraph_vector_init(&res, 0);
  igraph_real_t weights2[] = { 0,2,1, 0,5,2, 1,1,0, 0,2,8, 1,1,3, 1,1,4, 2,1 };
  igraph_vector_t weights_vec;
  igraph_small(&g, 10, IGRAPH_DIRECTED, 
	       0,1, 0,2, 0,3,    1,2, 1,4, 1,5,
	       2,3, 2,6,         3,2, 3,6,
	       4,5, 4,7,         5,6, 5,8, 5,9,
	       7,5, 7,8,         8,9,
	       5,2,
	       2,1,
	       -1);
  
   igraph_vector_view(&weights_vec, weights2, sizeof(weights2)/sizeof(igraph_real_t));
  igraph_i_minimum_spanning_tree_prim(&g, &res, &weights_vec);
  igraph_vector_print(&res);
  igraph_vector_destroy(&res);
  igraph_destroy(&g);
  
  return 0;
}
