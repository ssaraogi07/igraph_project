#include<igraph.h>
#include<stdio.h>

long int flag=0;
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
  igraph_vector_char_t visited;
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
  IGRAPH_CHECK(igraph_vector_char_init(&visited, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_char_destroy, &visited);
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
  igraph_vector_char_destroy(&visited);           \
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
	if(flag==1 && VECTOR(visited)[actvect]!=1)
      {
	VECTOR(visited)[actvect]=1;
	igraph_vector_shuffle(neis);
	//printf("neis was shuffled ");
      }
      //printf("printing neis contents : ");
	//igraph_vector_print(neis);
      /* Search for a neighbor that was not yet visited */
      igraph_bool_t any=0;
      long int nei;
      while (!any && (*ptr) <n) {
	nei=(long int) VECTOR(*neis)[(*ptr)];
	any=!VECTOR(added)[nei];
	//printf("printing ptr contents : %li",(*ptr));
	(*ptr) ++;
      }
      //printf("The while loop ends here ");
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
	//printf("Now the control comes here ");
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
  printf("%li",(long int) vid);
  return 0;
}	 

int main() {
  
  igraph_t graph;
  igraph_vector_t order, order_out, father, dist;
  long int i,j,n,c,p,x,y,d[105];
  scanf("%li %li",&n,&c);
  long int mat[n][n];
  igraph_tree(&graph,(igraph_integer_t) n,(igraph_integer_t) c, IGRAPH_TREE_UNDIRECTED);
  printf("Enter the node ");
  scanf("%li",&x);
  igraph_vector_init(&order, 0);
  igraph_vector_init(&order_out, 0);
  igraph_vector_init(&father, 0);
  igraph_vector_init(&dist, 0);

  flag=0;
  for(i=0;i<n;i++)
  {
	dfs_d(&graph, /*root=*/ (igraph_integer_t)i, /*neimode=*/ IGRAPH_ALL, 
	     /*unreachable=*/ 1,
	     &order, 0, 0, &dist, 0, 0, 0);
	for(j=0;j<n;j++)
	{
		mat[i][j]=(long int)VECTOR(dist)[j];
	}
   }
  flag=1;
  y=0;
  for(i=0;i<n;i++)
  {
	for(j=0;j<n;j++)
	 printf("%li",mat[i][j]);
	printf("\n");
  }
  for(i=0;i<10;i++)
  {
  	dfs_d(&graph, /*root=*/ 0, /*neimode=*/ IGRAPH_ALL, /*unreachable=*/ 1,
	     &order, 0, 0, &dist, 0, 0, 0);
	j=0;p=0;
	igraph_vector_print(&order);
	while(j<n && x!=(long int)VECTOR(order)[j])
	{
		p=p+mat[(long int)VECTOR(order)[j]][(long int)VECTOR(order)[j+1]];
		j++;
	}
	d[i]=p;
	y=y+p;
  }
  y=y/10;
  printf("Printing the distance array\n");
  for(i=0;i<10;i++)
  printf("%li ",d[i]);
  printf("\n");
  for(i=0;i<n;i++)
  {
	for(j=0;j<n;j++)
	 printf("%li",mat[i][j]);
	printf("\n");
  }
  printf("\nAverage path length %li",y);
  
  //printf("\n");
  igraph_destroy(&graph);
   igraph_vector_destroy(&order);
  igraph_vector_destroy(&order_out);
  igraph_vector_destroy(&father);
  igraph_vector_destroy(&dist);
  return 0;
}
