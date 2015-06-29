
#include <igraph.h>
#include <cmath>
#include "igraph_interface.h"
#include "igraph_community.h"
#include "igraph_interrupt_internal.h"


#include "infomap_Node.h"
#include "infomap_Greedy.h"
void swap (int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}
 
// A utility function to print an array
void printArray (int arr[], int n)
{
    for (int i = 0; i < n; i++)
        printf("%d ", arr[i]);
    printf("\n");
}
 
// A function to generate a random permutation of arr[]
void randomize ( int arr[], int n )
{
    // Use a different seed value so that we don't get same
    // result each time we run this program
    srand ( time(NULL) );
 
    // Start from the last element and swap one by one. We don't
    // need to run for the first element that's why i > 0
    for (int i = n-1; i > 0; i--)
    {
        // Pick a random index from 0 to i
        int j = rand() % (i+1);
 
        // Swap arr[i] with the element at random index
        swap(&arr[i], &arr[j]);
    }
}
/****************************************************************************/
int infomap_partition(FlowGraph * fgraph, bool rcall) {
  Greedy * greedy;

  // save the original graph
  FlowGraph * cpy_fgraph = new FlowGraph(fgraph);
  IGRAPH_FINALLY(delete_FlowGraph, cpy_fgraph);
  
  int Nnode = cpy_fgraph->Nnode; 
  // "real" number of vertex, ie. number of vertex of the graph	

  int iteration = 0;
  double outer_oldCodeLength, newCodeLength;
  
  int *initial_move = NULL;
  bool initial_move_done = true;
  
  do { // Main loop
    outer_oldCodeLength = fgraph->codeLength;
    
    if (iteration > 0) {
      /**********************************************************************/
      //  FIRST PART: re-split the network (if need)
      // ===========================================
      
      // intial_move indicate current clustering
      initial_move = new int[Nnode];    
      // new_cluster_id --> old_cluster_id (save curent clustering state)
      
      IGRAPH_FINALLY(operator delete [], initial_move);
      initial_move_done = false;
      
      int *subMoveTo = NULL; // enventual new partitionment of original graph
			
      if ((iteration % 2 == 0) && (fgraph->Nnode > 1)) { 
	// 0/ Submodule movements : partition each module of the 
	// current partition (rec. call)
	
	subMoveTo = new int[Nnode];       
	// vid_cpy_fgraph  --> new_cluster_id (new partition)

	IGRAPH_FINALLY(operator delete [], subMoveTo);
				
	int subModIndex = 0;

	for (int i=0 ; i < fgraph->Nnode ; i++) {
	  // partition each non trivial module
	  int sub_Nnode = fgraph->node[i]->members.size();
	  if (sub_Nnode > 1) { // If the module is not trivial
	    int *sub_members  = new int[sub_Nnode];      // id_sub --> id
	    IGRAPH_FINALLY(operator delete [], sub_members); 
	printArray(sub_members,sub_Nnode);    
	printf("sub_members after shuffling array :\n ");
	randomize(sub_members,sub_Nnode);
	printArray(sub_members,sub_Nnode);
	    for (int j=0 ; j < sub_Nnode ; j++)
	      sub_members[j] = fgraph->node[i]->members[j];
	    
	    // extraction of the subgraph
	    FlowGraph *sub_fgraph = new FlowGraph(cpy_fgraph, sub_Nnode, 
						  sub_members);
	    IGRAPH_FINALLY(delete_FlowGraph, sub_fgraph);
	    sub_fgraph->initiate();
	    
	    // recursif call of partitionment on the subgraph
	    infomap_partition(sub_fgraph, true);
	    
	    // Record membership changes
	    for (int j=0; j < sub_fgraph->Nnode; j++) {
	      int Nmembers = sub_fgraph->node[j]->members.size();
	      for (int k=0; k<Nmembers; k++) {
		subMoveTo[sub_members[sub_fgraph->node[j]->members[k]]] = 
		  subModIndex;
	      }
	      initial_move[subModIndex] = i;
	      subModIndex++;
	    }
	    
	    delete sub_fgraph;
	    IGRAPH_FINALLY_CLEAN(1);
	    delete [] sub_members;
	    IGRAPH_FINALLY_CLEAN(1);
	  } else{
	    subMoveTo[fgraph->node[i]->members[0]] = subModIndex;
	    initial_move[subModIndex] = i;
	    subModIndex++;
	  }
	}
      } else {
	// 1/ Single-node movements : allows each node to move (again)
	// save current modules
	for (int i=0; i < fgraph->Nnode; i++) { // for each module
	
	  int Nmembers = fgraph->node[i]->members.size(); // Module size
	  for (int j=0;j<Nmembers;j++) { // for each vertex (of the module)
	    initial_move[fgraph->node[i]->members[j]] = i;
	  }
	}
      }
      
      fgraph->back_to(cpy_fgraph);
      if (subMoveTo) {
	Greedy *cpy_greedy = new Greedy(fgraph);
	IGRAPH_FINALLY(delete_Greedy, cpy_greedy);
	
	cpy_greedy->setMove(subMoveTo);
	cpy_greedy->apply(false);
	
	delete_Greedy(cpy_greedy);
	IGRAPH_FINALLY_CLEAN(1);
	delete [] subMoveTo;
	IGRAPH_FINALLY_CLEAN(1);
      }
    }
    /**********************************************************************/
    //  SECOND PART: greedy optimizing it self
    // ===========================================
    double oldCodeLength;
    
    do {
      // greedy optimizing object creation
      greedy = new Greedy(fgraph);
      IGRAPH_FINALLY(delete_Greedy, greedy);
      
      // Initial move to apply ?
      if (!initial_move_done && initial_move) {
	initial_move_done = true;
	greedy->setMove(initial_move);
      }
      
      oldCodeLength = greedy->codeLength;
      bool moved = true;
      int Nloops = 0;
      //int count = 0;
      double inner_oldCodeLength = 1000;
      
      while (moved) { // main greedy optimizing loop
	inner_oldCodeLength = greedy->codeLength;
	moved = greedy->optimize();

	Nloops++;
	//count++;
	//printf("iteration : %d greedy: %lf inner: %lf",iteration,greedy->codeLength,inner_oldCodeLength);
	if (fabs(greedy->codeLength - inner_oldCodeLength) < 1.0e-10) 
	  // if the move does'n reduce the codelenght -> exit !
	  moved = false;
	
	//if (count == 10) {
	//	greedy->tune();
	//	count = 0;
	//}
      }
      
      // transform the network to network of modules:
      greedy->apply(true);
      newCodeLength = greedy->codeLength;
      
      // destroy greedy object
      delete greedy;
      IGRAPH_FINALLY_CLEAN(1);
      
    } while (oldCodeLength - newCodeLength >  1.0e-10); 
    // while there is some improvement
		
    if (iteration > 0) {
      delete [] initial_move;
      IGRAPH_FINALLY_CLEAN(1);
    }
    
    iteration++;
    if (!rcall) IGRAPH_ALLOW_INTERRUPTION();
  } while (outer_oldCodeLength - newCodeLength > 1.0e-10);
  
  delete cpy_fgraph;
  IGRAPH_FINALLY_CLEAN(1);
  return IGRAPH_SUCCESS;
}


int igraph_community_infomap(const igraph_t * graph,
                             const igraph_vector_t *e_weights,
                             const igraph_vector_t *v_weights,
                             int nb_trials,
                             igraph_vector_t *membership,
                             igraph_real_t *codelength) {

  FlowGraph * fgraph = new FlowGraph(graph, e_weights, v_weights);
  IGRAPH_FINALLY(delete_FlowGraph, fgraph);
	
  // compute stationary distribution
  fgraph->initiate();
	
  FlowGraph * cpy_fgraph ;
  double shortestCodeLength = 1000.0;
  
  // create membership vector
  int Nnode = fgraph->Nnode;
  IGRAPH_CHECK(igraph_vector_resize(membership, Nnode));
  
  for (int trial = 0; trial < nb_trials; trial++) {
    cpy_fgraph = new FlowGraph(fgraph);
    IGRAPH_FINALLY(delete_FlowGraph, cpy_fgraph);
    
    //partition the network
    IGRAPH_CHECK(infomap_partition(cpy_fgraph, false));
    
    // if better than the better...
    if (cpy_fgraph->codeLength < shortestCodeLength) {
      shortestCodeLength = cpy_fgraph->codeLength;
      // ... store the partition
      for (int i=0 ; i < cpy_fgraph->Nnode ; i++) {
	int Nmembers = cpy_fgraph->node[i]->members.size();
	for (int k=0; k < Nmembers; k++) {
	  //cluster[ cpy_fgraph->node[i]->members[k] ] = i;
	  VECTOR(*membership)[cpy_fgraph->node[i]->members[k]] = i;
	}
      }
    }
    
    delete_FlowGraph(cpy_fgraph);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  *codelength = (igraph_real_t) shortestCodeLength/log(2.0);
  
  delete fgraph;
  IGRAPH_FINALLY_CLEAN(1);
  return IGRAPH_SUCCESS;
}

int main() {
  igraph_t g;
  igraph_vector_t membership;
  igraph_real_t codelength;
  int i, j, k;

  igraph_vector_init(&membership,0);

  /* Unweighted test graph from the paper of Blondel et al */
  igraph_small(&g, 16, IGRAPH_UNDIRECTED,
      0, 2, 0, 3, 0, 4, 0, 5,
      1, 2, 1, 4, 1, 7,
      2, 4, 2, 5, 2, 6,
      3, 7,
      4, 10,
      5, 7, 5, 11,
      6, 7, 6, 11,
      8, 9, 8, 10, 8, 11, 8, 14, 8, 15,
      9, 12, 9, 14,
      10, 11, 10, 12, 10, 13, 10, 14, 11, 13, -1);
  FILE *fp;
  fp=fopen("abc.gml","w");
  igraph_write_graph_gml(&g,fp,NULL,0);
  igraph_community_infomap(&g, NULL, NULL, 1,&membership, &codelength);
  igraph_vector_print(&membership);
  printf("Coodelength : %li ",(long int)codelength);
  igraph_destroy(&g);
  igraph_vector_destroy(&membership);


#ifdef __APPLE__
  return 0;
#else
  return 77;
#endif
}
