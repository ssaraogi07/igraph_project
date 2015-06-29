#include <igraph.h>
#include <stdio.h>
int main() {
  
  igraph_t g, tree;
  igraph_vector_t  edges;
  long int i;
  //FILE *fp;
  //fp=fopen("qwe.gml","w");
  /*igraph_small(&g, 0, IGRAPH_UNDIRECTED, 
	       0,  1,  0,  2,  0,  3,  0,  4,  0,  5,
	       0,  6,  0,  7,  0,  8,  0, 10,  0, 11,
	       0, 12,  0, 13,  0, 17,  0, 19,  0, 21,
	       0, 31,  1,  2,  1,  3,  1,  7,  1, 13,
	       1, 17,  1, 19,  1, 21,  1, 30,  2,  3,
	       2,  7,  2,  8,  2,  9,  2, 13,  2, 27,
	       2, 28,  2, 32,  3,  7,  3, 12,  3, 13,
	       4,  6,  4, 10,  5,  6,  5, 10,  5, 16,
	       6, 16,  8, 30,  8, 32,  8, 33,  9, 33,
	       13, 33, 14, 32, 14, 33, 15, 32, 15, 33,
	       18, 32, 18, 33, 19, 33, 20, 32, 20, 33,
	       22, 32, 22, 33, 23, 25, 23, 27, 23, 29,
	       23, 32, 23, 33, 24, 25, 24, 27, 24, 31,
	       25, 31, 26, 29, 26, 33, 27, 33, 28, 31,
	       28, 33, 29, 32, 29, 33, 30, 32, 30, 33,
	       31, 32, 31, 33, 32, 33,
	       -1);*/
  igraph_real_t weights2[] = { 0,2,1, 0,5,2, 1,1,0, 0,2,8, 1,1,3, 1,1,4, 2,1 };
  igraph_vector_t weights_vec;
  igraph_small(&g, 10, IGRAPH_DIRECTED, 
	       0,1, 0,2, 0,3,    1,2, 1,4, 1,5,
	       2,3, 2,6,         3,2, 6,3,
	       4,5, 4,7,         5,6, 5,8, 5,9,
	       7,5, 7,8,         8,9,
	       5,2,
	       2,1,
	       -1);
  //igraph_write_graph_gml(&g, fp, 
	//		   NULL, 0);
   igraph_vector_view(&weights_vec, weights2, sizeof(weights2)/sizeof(igraph_real_t));
  igraph_minimum_spanning_tree_prim(&g, &tree, &weights_vec);
  igraph_write_graph_edgelist(&tree, stdout);

  igraph_vector_init(&edges, 0);
  igraph_minimum_spanning_tree(&g, &edges, &weights_vec);
  igraph_vector_print(&edges);
  igraph_vector_destroy(&edges);
  //fclose(fp);
  igraph_destroy(&tree);
  igraph_destroy(&g);
  
  return 0;
}
