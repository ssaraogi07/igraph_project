#include <igraph.h>
#include <time.h>
#include <math.h>
#include <string.h>
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
  char *added_edges;

  igraph_d_indheap_t heap=IGRAPH_D_INDHEAP_NULL;
  igraph_integer_t mode=IGRAPH_ALL;
  igraph_integer_t ed,ef,et;
  
  igraph_vector_t adj;
  igraph_vector_t temp;

  long int i, j, b, arr[5000],edno,cnt,nei,val,shf,max;

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
  IGRAPH_FINALLY(igraph_free, already_added);
  IGRAPH_CHECK(igraph_d_indheap_init(&heap, 0));
  IGRAPH_FINALLY(igraph_d_indheap_destroy, &heap);
  IGRAPH_VECTOR_INIT_FINALLY(&adj, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&temp, 0);
  memset(arr,0,sizeof(arr));
  for(b=0;b<no_of_edges;b++)
  {
	arr[(long int)VECTOR(*weights)[b]]++;
  }
  int ar[10],a;
   srand(time(NULL));
    for(a=0;a<10;a++)
    ar[a]=rand()%100+1;
    int x=0;    
  max=0;
  for(b=0;b<5000;b++)
  {
	if(arr[b]>max)
		max=arr[b];
  }
  shf=(long int)ceil(log2(max));
  srand(time(NULL));
  val=0;
  for(i=0;i<shf;i++)
	{
		val+=pow(2,i);
	}
  for(b=0;b<no_of_edges;b++)
  {
	VECTOR(*weights)[b]=(long int)VECTOR(*weights)[b]<<shf;
	if(val)
	VECTOR(*weights)[b]+=ar[x]%val +1;
 	x++;
	if(x==9)
	x=0;
  }
  cnt=0;
  i=rand()%no_of_nodes;
  for (; cnt<no_of_nodes; i++) {
    
    cnt++;
    if(i==no_of_nodes-1)
	i=0;
    if (already_added[i]>0) { continue; }
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

    while(! igraph_d_indheap_empty(&heap)) {
      /* Get minimal edge */
      long int from, edge,edge_w,a,x=0;
      igraph_integer_t tmp, to;
      igraph_d_indheap_max_index(&heap, &from, &edge);
      igraph_edge(graph, (igraph_integer_t) edge, &tmp, &to);
	
      /* Erase it */
      igraph_d_indheap_delete_max(&heap);

      /* Is this edge already included? */
	
      if (added_edges[edge]==0) {
        if (from==to) { to=tmp; }
        /* Does it point to a visited node? */ 
	//printf("EDGE: %li FROM: %li TO: %li",edge,from,(long int)to);     
        if (already_added[(long int)to]==0) {
          already_added[(long int)to]=1;
          added_edges[edge]=1;
          IGRAPH_CHECK(igraph_vector_push_back(res, edge));
          /* add all outgoing edges */
          igraph_incident(graph, &adj, to, (igraph_neimode_t) mode);
          for (j=0; j<igraph_vector_size(&adj); j++) {
            long int edgeno=(long int) VECTOR(adj)[j];
            igraph_integer_t edgefrom, edgeto;
            long int neighbor;
            igraph_edge(graph, (igraph_integer_t) edgeno, &edgefrom, &edgeto);
            neighbor= edgefrom != to ? edgefrom : edgeto;
	    nei= edgefrom == to ? edgefrom : edgeto;
            if (already_added[neighbor] == 0) {
              IGRAPH_CHECK(igraph_d_indheap_push(&heap, -VECTOR(*weights)[edgeno], to,
                   edgeno));
		
            }
          }
        } /* for */
       } /* if !already_added */
    } /* while in the same component */
   }
  /* for all nodes */

  igraph_d_indheap_destroy(&heap);
  igraph_Free(already_added);
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
