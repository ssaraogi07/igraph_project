#include <igraph.h>
#include <igraph_structural.h>
#include <igraph_dqueue.h>
#include <igraph_interface.h>
#include <igraph_interrupt_internal.h>
#include <igraph_memory.h>
#include <igraph_progress.h>
#include <igraph_types_internal.h>

int igraph_i_minimum_spanning_tree_unweighted(const igraph_t* graph,
    igraph_vector_t* res) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  char *already_added;
  char *added_edges;
  
  igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
  igraph_vector_t tmp=IGRAPH_VECTOR_NULL;
  long int i, j;

  igraph_vector_clear(res);

  added_edges=igraph_Calloc(no_of_edges, char);
  if (added_edges==0) {
    IGRAPH_ERROR("unweighted spanning tree failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, added_edges);
  already_added=igraph_Calloc(no_of_nodes, char);
  if (already_added==0) {
    IGRAPH_ERROR("unweighted spanning tree failed", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, already_added);
  IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
  IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
  
  for (i=0; i<no_of_nodes; i++) {
    if (already_added[i]>0) { continue; }

    IGRAPH_ALLOW_INTERRUPTION();

    already_added[i]=1;
    IGRAPH_CHECK(igraph_dqueue_push(&q, i));
    while (! igraph_dqueue_empty(&q)) {
      long int act_node=(long int) igraph_dqueue_pop(&q);
      IGRAPH_CHECK(igraph_incident(graph, &tmp, (igraph_integer_t) act_node,
				   IGRAPH_ALL));
	igraph_vector_sort(&tmp);
      for (j=0; j<igraph_vector_size(&tmp); j++) {
        long int edge=(long int) VECTOR(tmp)[j];
        if (added_edges[edge]==0) {
          igraph_integer_t from, to;
          igraph_edge(graph, (igraph_integer_t) edge, &from, &to);
          if (act_node==to) { to=from; }
          if (already_added[(long int) to]==0) {
            already_added[(long int) to]=1;
            added_edges[edge]=1;
            IGRAPH_CHECK(igraph_vector_push_back(res, edge));
            IGRAPH_CHECK(igraph_dqueue_push(&q, to));
          }
        }
      }
    }
  }
  
  igraph_dqueue_destroy(&q);
  igraph_Free(already_added);
  igraph_vector_destroy(&tmp);
  igraph_Free(added_edges);
  IGRAPH_FINALLY_CLEAN(4);

  return IGRAPH_SUCCESS;
}
int main() {
  
  igraph_t g;
  igraph_vector_t  res;
  long int i;
  
  igraph_vector_init(&res,0); 
  igraph_small(&g, 10, IGRAPH_DIRECTED, 
	       0,1, 0,2, 0,3,    1,2, 1,4, 1,5,
	       2,3, 2,6,         3,2, 6,3,
	       4,5, 4,7,         5,6, 5,8, 5,9,
	       7,5, 7,8,         8,9,
	       5,2,
	       2,1,
	       -1);
 
  igraph_i_minimum_spanning_tree_unweighted(&g,&res);
  printf("The required mst is: ");
  igraph_vector_print(&res);
  igraph_vector_destroy(&res);
  igraph_destroy(&g);
  
  return 0;
}
