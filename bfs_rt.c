#include <igraph.h>
#include<stdio.h>
void vector_print(igraph_vector_t *v) {
  long int i;
  for (i=0; i<igraph_vector_size(v); i++) {
    printf(" %li", (long int) VECTOR(*v)[i]);
  }
  printf("\n");
}

long int flag=0;
int bfsd(const igraph_t *graph, 
	       igraph_integer_t root,igraph_vector_t *roots,
	       igraph_neimode_t mode, igraph_bool_t unreachable,
	       const igraph_vector_t *restricted,
	       igraph_vector_t *order, igraph_vector_t *rank,
	       igraph_vector_t *father,
	       igraph_vector_t *pred, igraph_vector_t *succ,
	       igraph_vector_t *dist, igraph_bfshandler_t *callback,
	       void *extra) {
  
  igraph_dqueue_t Q;
  long int no_of_nodes=igraph_vcount(graph);
  long int actroot=0;
  igraph_vector_char_t added;

  igraph_lazy_adjlist_t adjlist;
  
  long int act_rank=0;
  long int pred_vec=-1;
  
  long int rootpos=0;
  long int noroots= roots ? igraph_vector_size(roots) : 1;

  if (!roots && (root < 0 || root >= no_of_nodes)) {
    IGRAPH_ERROR("Invalid root vertex in BFS", IGRAPH_EINVAL);
  }
  
  if (roots) {
    igraph_real_t min, max;
    igraph_vector_minmax(roots, &min, &max);
    if (min < 0 || max >= no_of_nodes) {
      IGRAPH_ERROR("Invalid root vertex in BFS", IGRAPH_EINVAL);
    }
  }

  if (restricted) {
    igraph_real_t min, max;
    igraph_vector_minmax(restricted, &min, &max);
    if (min < 0 || max >= no_of_nodes) {
      IGRAPH_ERROR("Invalid vertex id in restricted set", IGRAPH_EINVAL);
    }
  }

  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
  }
  
  if (!igraph_is_directed(graph)) { mode=IGRAPH_ALL; }

  IGRAPH_CHECK(igraph_vector_char_init(&added, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_char_destroy, &added);
  IGRAPH_CHECK(igraph_dqueue_init(&Q, 100));
  IGRAPH_FINALLY(igraph_dqueue_destroy, &Q);

  IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, mode, /*simplify=*/ 0));
  IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);

  /* Mark the vertices that are not in the restricted set, as already
     found. Special care must be taken for vertices that are not in
     the restricted set, but are to be used as 'root' vertices. */
  if (restricted) {
    long int i, n=igraph_vector_size(restricted);
    igraph_vector_char_fill(&added, 1);
    for (i=0; i<n; i++) {
      long int v=(long int) VECTOR(*restricted)[i];
      VECTOR(added)[v]=0;
    }
  }

  /* Resize result vectors, and fill them with IGRAPH_NAN */

# define VINIT(v) if (v) {                      \
  igraph_vector_resize((v), no_of_nodes);	\
  igraph_vector_fill((v), IGRAPH_NAN); }

  VINIT(order);
  VINIT(rank);
  VINIT(father);
  VINIT(pred);
  VINIT(succ);
  VINIT(dist);
# undef VINIT

  //vector_print(roots);

  while (1) { 

    /* Get the next root vertex, if any */

    if (roots && rootpos < noroots) {
      /* We are still going through the 'roots' vector */
      actroot=(long int) VECTOR(*roots)[rootpos++];
    } else if (!roots && rootpos==0) {
      /* We have a single root vertex given, and start now */
      actroot=root;
      rootpos++;
    } else if (rootpos==noroots && unreachable) {
      /* We finished the given root(s), but other vertices are also
	 tried as root */
      actroot=0;
      rootpos++;
    } else if (unreachable && actroot+1 < no_of_nodes) {
      /* We are already doing the other vertices, take the next one */
      actroot++;
    } else {
      /* No more root nodes to do */
      break;
    }

    /* OK, we have a new root, start BFS */
    if (VECTOR(added)[actroot]) { continue; }
    IGRAPH_CHECK(igraph_dqueue_push(&Q, actroot));
    IGRAPH_CHECK(igraph_dqueue_push(&Q, 0));
    VECTOR(added)[actroot] = 1;
    if (father) { VECTOR(*father)[actroot] = -1; }
      
    pred_vec=-1;

    while (!igraph_dqueue_empty(&Q)) {
      long int actvect=(long int) igraph_dqueue_pop(&Q);
      long int actdist=(long int) igraph_dqueue_pop(&Q);
      long int succ_vec;
      igraph_vector_t *neis=igraph_lazy_adjlist_get(&adjlist, 
						    (igraph_integer_t) actvect);
      long int i, n=igraph_vector_size(neis);  
      if(flag==1)  
       igraph_vector_shuffle(neis);
      if (pred) { VECTOR(*pred)[actvect] = pred_vec; }
      if (rank) { VECTOR(*rank) [actvect] = act_rank; }
      if (order) { VECTOR(*order)[act_rank++] = actvect; }
      if (dist) { VECTOR(*dist)[actvect] = actdist; }   
      //printf("printing neis...");  
      igraph_vector_shuffle(neis);
      for (i=0; i<n; i++) {
	long int nei=(long int) VECTOR(*neis)[i];
	if (! VECTOR(added)[nei]) {
	  VECTOR(added)[nei] = 1;
	  IGRAPH_CHECK(igraph_dqueue_push(&Q, nei));
	  IGRAPH_CHECK(igraph_dqueue_push(&Q, actdist+1));
	  if (father) { VECTOR(*father)[nei] = actvect; }
	}
      }

      succ_vec = igraph_dqueue_empty(&Q) ? -1L : 
	(long int) igraph_dqueue_head(&Q);
      if (callback) {
	igraph_bool_t terminate=
	  callback(graph, (igraph_integer_t) actvect, (igraph_integer_t) 
		   pred_vec, (igraph_integer_t) succ_vec, 
		   (igraph_integer_t) act_rank-1, (igraph_integer_t) actdist, 
		   extra);
	if (terminate) {
	  igraph_lazy_adjlist_destroy(&adjlist);
	  igraph_dqueue_destroy(&Q);
	  igraph_vector_char_destroy(&added);
	  IGRAPH_FINALLY_CLEAN(3);
	  return 0;
	}
      }

      if (succ) { VECTOR(*succ)[actvect] = succ_vec; }
      pred_vec=actvect;

    } /* while Q !empty */

  } /* for actroot < no_of_nodes */

  igraph_lazy_adjlist_destroy(&adjlist);
  igraph_dqueue_destroy(&Q);
  igraph_vector_char_destroy(&added);
  IGRAPH_FINALLY_CLEAN(3);
  
  return 0;
}


igraph_bool_t bfs_callback(const igraph_t *graph,
			   igraph_integer_t vid, 
			   igraph_integer_t pred, 
			   igraph_integer_t succ,
			   igraph_integer_t rank,
			   igraph_integer_t dist,
			   void *extra) {
  printf(" %li", (long int) vid);
  return 0;
}		   

int main() {
  
  igraph_t graph;
  igraph_vector_t order,dist;
  long int i,j,n,c,p,x,y,d[105];
  scanf("%li %li",&n,&c);
  long int mat[n][n];
  igraph_tree(&graph,(igraph_integer_t) n,(igraph_integer_t) c, IGRAPH_TREE_UNDIRECTED);
  printf("Enter the node ");
  scanf("%li",&x);
  igraph_vector_init(&order, 0);
  igraph_vector_init(&dist, 0);
  FILE *fp=fopen("qwe.gml","w");
  
  bfsd(&graph, /*root=*/0, /*roots=*/ 0, /*neimode=*/ IGRAPH_ALL, 
	     /*unreachable=*/ 1, /*restricted=*/ 0,
	     &order, 0, 0, 0, 0, &dist, 
	     /*callback=*/ 0, /*extra=*/ 0);
  
  vector_print(&order);
  igraph_vector_print(&dist);
  flag=0;
  for(i=0;i<n;i++)
  {
	bfsd(&graph, /*root=*/i, /*roots=*/ 0, /*neimode=*/ IGRAPH_ALL, 
	     /*unreachable=*/ 1, /*restricted=*/ 0,
	     &order, 0, 0, 0, 0, &dist, 
	     /*callback=*/ 0, /*extra=*/ 0);
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
  igraph_write_graph_gml(&graph,fp,NULL,"test suite");
  for(i=0;i<10;i++)
  {
  	bfsd(&graph, /*root=*/0, /*roots=*/ 0, /*neimode=*/ IGRAPH_ALL, 
	     /*unreachable=*/ 1, /*restricted=*/ 0,
	     &order, 0, 0, 0, 0, &dist, 
	     /*callback=*/ 0, /*extra=*/ 0);
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
  
  fclose(fp);
  igraph_vector_destroy(&order);
  igraph_vector_destroy(&dist);
  igraph_destroy(&graph);
  
  return 0;
}
