#include<igraph.h>

int dfs_d(const igraph_t *graph, igraph_integer_t root,
	       igraph_neimode_t mode, igraph_bool_t unreachable, 
	       igraph_vector_t *order,
	       igraph_vector_t *order_out, igraph_vector_t *father,
	       igraph_vector_t *dist, igraph_dfshandler_t *in_callback,
	       igraph_dfshandler_t *out_callback,
	       void *extra) {
  
  long int no_of_nodes=igraph_vcount(graph);
  igraph_lazy_adjlist_t adjlist;
  igraph_stack_t stack;
  igraph_vector_char_t added;
  igraph_vector_long_t nptr;
  long int actroot;
  long int act_rank=0;
  long int rank_out=0;
  long int act_dist=0;

  if (root < 0 || root >= no_of_nodes) { 
    IGRAPH_ERROR("Invalid root vertex for DFS", IGRAPH_EINVAL);
  }

  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
  }
  
  if (!igraph_is_directed(graph)) { mode=IGRAPH_ALL; }

  IGRAPH_CHECK(igraph_vector_char_init(&added, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_char_destroy, &added);
  IGRAPH_CHECK(igraph_stack_init(&stack, 100));
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);
  IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, mode, /*simplify=*/ 0));  
  IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);
  IGRAPH_CHECK(igraph_vector_long_init(&nptr, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &nptr);

# define FREE_ALL() do { 			\
  igraph_vector_long_destroy(&nptr);            \
  igraph_lazy_adjlist_destroy(&adjlist);        \
  igraph_stack_destroy(&stack);                 \
  igraph_vector_char_destroy(&added);           \
  IGRAPH_FINALLY_CLEAN(4); } while (0)

  /* Resize result vectors and fill them with IGRAPH_NAN */
  
# define VINIT(v) if (v) {                      \
    igraph_vector_resize(v, no_of_nodes);       \
    igraph_vector_fill(v, IGRAPH_NAN); }
  
  VINIT(order);
  VINIT(order_out);
  VINIT(father);
  VINIT(dist);

# undef VINIT

  IGRAPH_CHECK(igraph_stack_push(&stack, root));
  VECTOR(added)[(long int)root] = 1;
  if (father) { VECTOR(*father)[(long int)root] = -1; }
  if (order) { VECTOR(*order)[act_rank++] = root; }
  if (dist) { VECTOR(*dist)[(long int)root] = 0; }
  if (in_callback) {
    igraph_bool_t terminate=in_callback(graph, root, 0, extra);
    if (terminate) { FREE_ALL(); return 0; }
  }

  for (actroot=0; actroot<no_of_nodes; actroot++) {

    /* 'root' first, then all other vertices */
    if (igraph_stack_empty(&stack)) {
      if (!unreachable) { break; }
      if (VECTOR(added)[actroot]) { continue; }
      IGRAPH_CHECK(igraph_stack_push(&stack, actroot));
      VECTOR(added)[actroot] = 1;
      if (father) { VECTOR(*father)[actroot] = -1; }
      if (order) { VECTOR(*order)[act_rank++] = actroot; }
      if (dist) { VECTOR(*dist)[actroot] = 0; }

      if (in_callback) {
	igraph_bool_t terminate=in_callback(graph, (igraph_integer_t) actroot,
					    0, extra);
	if (terminate) { FREE_ALL(); return 0; }
      }
    }
    
    while (!igraph_stack_empty(&stack)) {
      long int actvect=(long int) igraph_stack_top(&stack);
      igraph_vector_t *neis=igraph_lazy_adjlist_get(&adjlist, 
						    (igraph_integer_t) actvect);
      long int n=igraph_vector_size(neis);
      long int *ptr=igraph_vector_long_e_ptr(&nptr, actvect);

      igraph_vector_sort(neis);
      /* Search for a neighbor that was not yet visited */
      igraph_bool_t any=0;
      long int nei;
      while (!any && (*ptr) <n) {
	nei=(long int) VECTOR(*neis)[(*ptr)];
	any=!VECTOR(added)[nei];
	(*ptr) ++;
      }
      if (any) {
	/* There is such a neighbor, add it */
	IGRAPH_CHECK(igraph_stack_push(&stack, nei));
	VECTOR(added)[nei] = 1;
	if (father) { VECTOR(*father)[ nei ] = actvect; }
	if (order) { VECTOR(*order)[act_rank++] = nei; }
	act_dist++;
	if (dist) { VECTOR(*dist)[nei] = act_dist; }

	if (in_callback) {
	  igraph_bool_t terminate=in_callback(graph, (igraph_integer_t) nei,
					      (igraph_integer_t) act_dist, 
					      extra);
	  if (terminate) { FREE_ALL(); return 0; }
	}

      } else {
	/* There is no such neighbor, finished with the subtree */
	igraph_stack_pop(&stack);
	if (order_out) { VECTOR(*order_out)[rank_out++] = actvect; }
	act_dist--;

	if (out_callback) {
	  igraph_bool_t terminate=out_callback(graph, (igraph_integer_t) 
					       actvect, (igraph_integer_t) 
					       act_dist, extra);
	  if (terminate) { FREE_ALL(); return 0; }
	}
      }
    }      
  }

  FREE_ALL();
# undef FREE_ALL

  return 0;
}

igraph_bool_t dfs_callback(const igraph_t *graph,
			   igraph_integer_t vid, 
			   igraph_integer_t dist,
			   void *extra) {
  printf(" %li", (long int) vid);
  return 0;
}		   

int main() {
  
  igraph_t graph, ring;
  igraph_vector_t order, order_out, father, dist;
  long int i;
  
  igraph_ring(&ring, 10, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/ 1);
  igraph_disjoint_union(&graph, &ring, &ring);
  igraph_destroy(&ring);
  
  igraph_vector_init(&order, 0);
  igraph_vector_init(&order_out, 0);
  igraph_vector_init(&father, 0);
  //igraph_vector_init(&pred, 0);
  //igraph_vector_init(&succ, 0);
  igraph_vector_init(&dist, 0);
  
  dfs_d(&graph, /*root=*/0,/*neimode=*/ IGRAPH_OUT, 
	     /*unreachable=*/ 1,&order,&order_out, &father, &dist, 
	     /*in_callback=*/ 0, /*out_callback=*/ 0, /*extra=*/ 0);
  
  igraph_vector_print(&order);
  igraph_vector_print(&order_out);
  igraph_vector_print(&father);
  igraph_vector_print(&dist);

  igraph_vector_destroy(&order);
  igraph_vector_destroy(&order_out);
  igraph_vector_destroy(&father);
  igraph_vector_destroy(&dist);

  /* Test the callback */

  dfs_d(&graph, /*root=*/ 0, /*neimode=*/ IGRAPH_OUT, 
	     /*unreachable=*/ 1,
	     0, 0, 0, 0, &dfs_callback, &dfs_callback, 0);
  printf("\n");
  
  /* Test different roots */

  dfs_d(&graph, /*root=*/ 2, /*neimode=*/ IGRAPH_OUT, 
	     /*unreachable=*/ 1,
	     0, 0, 0, 0, &dfs_callback, 0, 0);
  printf("\n");

  dfs_d(&graph, /*root=*/ 2, /*neimode=*/ IGRAPH_OUT, 
	     /*unreachable=*/ 1,
	     0, 0, 0, 0, 0, &dfs_callback, 0);
  printf("\n");
  igraph_destroy(&graph);
  
  return 0;
}
