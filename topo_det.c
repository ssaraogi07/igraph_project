#include <igraph.h>
#include "igraph_structural.h"
#include "igraph_transitivity.h"
#include "igraph_paths.h"
#include "igraph_math.h"
#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_adjlist.h"
#include "igraph_interface.h"
#include "igraph_progress.h"
#include "igraph_interrupt_internal.h"
#include "igraph_centrality.h"
#include "igraph_components.h"
#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_types_internal.h"
#include "igraph_dqueue.h"
#include "igraph_attributes.h"
#include "igraph_neighborhood.h"
#include "igraph_topology.h"
#include "igraph_qsort.h"
#include "config.h"



int igraph_topological_sorting(const igraph_t* graph, igraph_vector_t *res,
			       igraph_neimode_t mode) {
  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t degrees, neis,temp;
  igraph_dqueue_t sources;
  igraph_neimode_t deg_mode;
  long int node, i, j, b, x;

  if (mode == IGRAPH_ALL || !igraph_is_directed(graph)) {
    IGRAPH_ERROR("topological sorting does not make sense for undirected graphs", IGRAPH_EINVAL);
  } else if (mode == IGRAPH_OUT) {
    deg_mode = IGRAPH_IN;
  } else if (mode == IGRAPH_IN) {
    deg_mode = IGRAPH_OUT;
  } else {
    IGRAPH_ERROR("invalid mode", IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&degrees, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&temp, 0);
  IGRAPH_CHECK(igraph_dqueue_init(&sources, 0));
  IGRAPH_FINALLY(igraph_dqueue_destroy, &sources);
  IGRAPH_CHECK(igraph_degree(graph, &degrees, igraph_vss_all(), deg_mode, 0));

  igraph_vector_clear(res);
  x=0;
  /* Do we have nodes with no incoming vertices? */
  for (i=0; i<no_of_nodes; i++) {
    if (VECTOR(degrees)[i] == 0)
      /*IGRAPH_CHECK(igraph_dqueue_push(&sources, i));*/
	{
	igraph_vector_insert(&temp,x,(igraph_real_t)i);
	x++;
	}
  }
 long int sz=igraph_vector_size(&temp);
 igraph_vector_sort(&temp);
 for(b=0;b<sz;b++)
 {
	IGRAPH_CHECK(igraph_dqueue_push(&sources, (long int)VECTOR(temp)[b]));
 }
  /* Take all nodes with no incoming vertices and remove them */
  while (!igraph_dqueue_empty(&sources)) {
    igraph_real_t tmp=igraph_dqueue_pop(&sources); node=(long) tmp;
    /* Add the node to the result vector */
    igraph_vector_push_back(res, node);
    /* Exclude the node from further source searches */
    VECTOR(degrees)[node]=-1;
    /* Get the neighbors and decrease their degrees by one */
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) node, mode));
    j=igraph_vector_size(&neis);
    igraph_vector_sort(&neis);
    for (i=0; i<j; i++) {
      VECTOR(degrees)[(long)VECTOR(neis)[i]]--;
      if (VECTOR(degrees)[(long)VECTOR(neis)[i]] == 0)
	IGRAPH_CHECK(igraph_dqueue_push(&sources, VECTOR(neis)[i]));
    }
  }

  if (igraph_vector_size(res)<no_of_nodes)
    IGRAPH_WARNING("graph contains a cycle, partial result is returned");
  igraph_vector_destroy(&degrees);
  igraph_vector_destroy(&neis);
  igraph_vector_destroy(&temp);
  igraph_dqueue_destroy(&sources);
  IGRAPH_FINALLY_CLEAN(4);

  return 0;
}

void print_vector(igraph_vector_t *v, FILE *f) {
  long int i;
  for (i=0; i<igraph_vector_size(v); i++) {
    fprintf(f, " %d", (int)VECTOR(*v)[i]);
  }
  fprintf(f, "\n");
}

void igraph_warning_handler_print_stdout(const char *reason, const char *file,
					 int line, int igraph_errno) {
  fprintf(stdout, "Warning: %s\n", reason);
}


int main() {

  igraph_t g;
  igraph_vector_t v, res;
  igraph_bool_t is_dag;
  int ret;

  igraph_set_warning_handler(igraph_warning_handler_print_stdout);

  /* Test graph taken from http://en.wikipedia.org/wiki/Topological_sorting
   * @ 05.03.2006 */
  igraph_small(&g, 8, 1, 0, 3, 0, 4, 1, 3, 2, 4, 2, 7, \
	       3, 5, 3, 6, 3, 7, 4, 6, -1);

  igraph_vector_init(&res, 0);

  igraph_is_dag(&g, &is_dag);
  if (!is_dag)
    return 2;

  igraph_topological_sorting(&g, &res, IGRAPH_OUT);
  print_vector(&res, stdout);
  igraph_topological_sorting(&g, &res, IGRAPH_IN);
  print_vector(&res, stdout);

  /* Add a circle: 5 -> 0 */
  igraph_vector_init_int(&v, 2, 5, 0);
  igraph_add_edges(&g, &v, 0);
  igraph_is_dag(&g, &is_dag);
  if (is_dag)
    return 3;
  igraph_topological_sorting(&g, &res, IGRAPH_OUT);
  print_vector(&res, stdout);
  igraph_vector_destroy(&v);
  igraph_destroy(&g);

  /* Error handling */
  igraph_set_error_handler(igraph_error_handler_ignore);

  /* This graph is the same but undirected */
  igraph_small(&g, 8, 0, 0, 3, 0, 4, 1, 3, 2, 4, 2, 7, \
	       3, 5, 3, 6, 3, 7, 4, 6, -1);
  igraph_is_dag(&g, &is_dag);
  if (is_dag)
    return 4;
  ret=igraph_topological_sorting(&g, &res, IGRAPH_ALL);
  if (ret != IGRAPH_EINVAL) {
    return 1;
  }
  ret=igraph_topological_sorting(&g, &res, IGRAPH_OUT);
  if (ret != IGRAPH_EINVAL) {
    return 1;
  }  
  igraph_destroy(&g);

  igraph_vector_destroy(&res);  
  return 0;
}
