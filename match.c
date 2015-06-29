#include <igraph.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "config.h"
#include "igraph_adjlist.h"
#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_dqueue.h"
#include "igraph_flow.h"
#include "igraph_interface.h"
#include "igraph_matching.h"
#include "igraph_structural.h"

/* #define MATCHING_DEBUG */

#ifdef _MSC_VER
/* MSVC does not support variadic macros */
#include <stdarg.h>
static void debug(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
#ifdef MATCHING_DEBUG
  vfprintf(stderr, fmt, args);
#endif
  va_end(args);
}
#else
#  ifdef MATCHING_DEBUG
#    define debug(...) fprintf(stderr, __VA_ARGS__)
#  else
#    define debug(...) 
#  endif
#endif

#define MATCHED(v) (VECTOR(match)[v] != -1)
#define UNMATCHED(v) (!MATCHED(v))

int igraph_i_maximum_bipartite_matching_unweighted(const igraph_t* graph,
    const igraph_vector_bool_t* types, igraph_integer_t* matching_size,
    igraph_vector_long_t* matching) {
  long int i, j, k, n, no_of_nodes = igraph_vcount(graph);
  long int num_matched;             /* number of matched vertex pairs */
  igraph_vector_long_t match;       /* will store the matching */
  igraph_vector_t labels;           /* will store the labels */
  igraph_vector_t neis;             /* used to retrieve the neighbors of a node */
  igraph_dqueue_long_t q;           /* a FIFO for push ordering */
  igraph_bool_t smaller_set;        /* denotes which part of the bipartite graph is smaller */
  long int label_changed = 0;       /* Counter to decide when to run a global relabeling */
  long int relabeling_freq = no_of_nodes / 2;

  /* We will use:
   * - FIFO push ordering
   * - global relabeling frequency: n/2 steps where n is the number of nodes
   * - simple greedy matching for initialization
   */

  /* (1) Initialize data structures */
  IGRAPH_CHECK(igraph_vector_long_init(&match, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &match);
  IGRAPH_VECTOR_INIT_FINALLY(&labels, no_of_nodes);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_dqueue_long_init(&q, 0));
  IGRAPH_FINALLY(igraph_dqueue_long_destroy, &q);

  /* (2) Initially, every node is unmatched */
  igraph_vector_long_fill(&match, -1);

  /* (3) Find an initial matching in a greedy manner.
   *     At the same time, find which side of the graph is smaller. */
  num_matched = 0; j = 0;
  for (i = 0; i < no_of_nodes; i++) {
    if (VECTOR(*types)[i])
      j++;
    if (MATCHED(i))
      continue;
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) i,
				  IGRAPH_ALL));
    n = igraph_vector_size(&neis);
    for (j = 0; j < n; j++) {
      k = (long int) VECTOR(neis)[j];
      if (UNMATCHED(k)) {
        /* We match vertex i to vertex VECTOR(neis)[j] */
        VECTOR(match)[k] = i;
        VECTOR(match)[i] = k;
        num_matched++;
        break;
      }
    }
  }
  smaller_set = (j <= no_of_nodes/2);

  /* (4) Set the initial labeling -- lines 1 and 2 in the tech report */
  IGRAPH_CHECK(igraph_i_maximum_bipartite_matching_unweighted_relabel(
      graph, types, &labels, &match, smaller_set));

  /* (5) Fill the push queue with the unmatched nodes from the smaller set. */
  for (i = 0; i < no_of_nodes; i++) {
    if (UNMATCHED(i) && VECTOR(*types)[i] == smaller_set)
      IGRAPH_CHECK(igraph_dqueue_long_push(&q, i));
  }

  /* (6) Main loop from the referenced tech report -- lines 4--13 */
  label_changed = 0;
  while (!igraph_dqueue_long_empty(&q)) {
    long int v = igraph_dqueue_long_pop(&q);             /* Line 13 */
    long int u = -1, label_u = 2 * no_of_nodes;
    long int w;

    if (label_changed >= relabeling_freq) {
      /* Run global relabeling */
      IGRAPH_CHECK(igraph_i_maximum_bipartite_matching_unweighted_relabel(
            graph, types, &labels, &match, smaller_set));
      label_changed = 0;
    }

    debug("Considering vertex %ld\n", v);

    /* Line 5: find row u among the neighbors of v s.t. label(u) is minimal */
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) v, 
				  IGRAPH_ALL));
    n = igraph_vector_size(&neis);
    for (i = 0; i < n; i++) {
      if (VECTOR(labels)[(long int)VECTOR(neis)[i]] < label_u) {
        u = (long int) VECTOR(neis)[i];
        label_u = (long int) VECTOR(labels)[u];
        label_changed++;
      }
    }

    debug("  Neighbor with smallest label: %ld (label=%ld)\n", u, label_u);

    if (label_u < no_of_nodes) {                         /* Line 6 */
      VECTOR(labels)[v] = VECTOR(labels)[u] + 1;         /* Line 7 */
      if (MATCHED(u)) {                                  /* Line 8 */
        w = VECTOR(match)[u];
        debug("  Vertex %ld is matched to %ld, performing a double push\n", u, w);
        if (w != v) {
          VECTOR(match)[u] = -1; VECTOR(match)[w] = -1;  /* Line 9 */
          IGRAPH_CHECK(igraph_dqueue_long_push(&q, w));  /* Line 10 */
          debug("  Unmatching & activating vertex %ld\n", w);
          num_matched--;
        }
      }
      VECTOR(match)[u] = v; VECTOR(match)[v] = u;      /* Line 11 */
      num_matched++;
      VECTOR(labels)[u] += 2;                          /* Line 12 */
      label_changed++;
    }
	printf("MATCH: ");
	igraph_vector_long_print(&match);
	printf("LABELS ");
	igraph_vector_print(&labels);
  }

  /* Fill the output parameters */
  if (matching != 0) {
    IGRAPH_CHECK(igraph_vector_long_update(matching, &match));
  }
  if (matching_size != 0) {
    *matching_size = (igraph_integer_t) num_matched;
  }

  /* Release everything */
  igraph_dqueue_long_destroy(&q);
  igraph_vector_destroy(&neis);
  igraph_vector_destroy(&labels);
  igraph_vector_long_destroy(&match);
  IGRAPH_FINALLY_CLEAN(4);

  return IGRAPH_SUCCESS;
}

int igraph_i_maximum_bipartite_matching_unweighted_relabel(const igraph_t* graph,
    const igraph_vector_bool_t* types, igraph_vector_t* labels,
    igraph_vector_long_t* match, igraph_bool_t smaller_set) {
  long int i, j, n, no_of_nodes = igraph_vcount(graph), matched_to;
  igraph_dqueue_long_t q;
  igraph_vector_t neis;

  debug("Running global relabeling.\n");

  /* Set all the labels to no_of_nodes first */
  igraph_vector_fill(labels, no_of_nodes);

  /* Allocate vector for neighbors */
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

  /* Create a FIFO for the BFS and initialize it with the unmatched rows
   * (i.e. members of the larger set) */
  IGRAPH_CHECK(igraph_dqueue_long_init(&q, 0));
  IGRAPH_FINALLY(igraph_dqueue_long_destroy, &q);
  for (i = 0; i < no_of_nodes; i++) {
    if (VECTOR(*types)[i] != smaller_set && VECTOR(*match)[i] == -1) {
      IGRAPH_CHECK(igraph_dqueue_long_push(&q, i));
      VECTOR(*labels)[i] = 0;
    }
  }

  /* Run the BFS */
  while (!igraph_dqueue_long_empty(&q)) {
    long int v = igraph_dqueue_long_pop(&q);
    long int w;

    IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) v,
				  IGRAPH_ALL));

    n = igraph_vector_size(&neis);
    //igraph_vector_shuffle(&neis);
    for (j = 0; j < n; j++) {
      w = (long int) VECTOR(neis)[j];
      if (VECTOR(*labels)[w] == no_of_nodes) {
        VECTOR(*labels)[w] = VECTOR(*labels)[v] + 1;
        matched_to = VECTOR(*match)[w];
        if (matched_to != -1 && VECTOR(*labels)[matched_to] == no_of_nodes) {
          IGRAPH_CHECK(igraph_dqueue_long_push(&q, matched_to));
          VECTOR(*labels)[matched_to] = VECTOR(*labels)[w] + 1;
        }
      }
    }
  }
  printf("Inside relabel : ");
  igraph_vector_print(labels);
  igraph_dqueue_long_destroy(&q);
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(2);

  return IGRAPH_SUCCESS;
}
int test_graph_from_leda_tutorial() {
  /* Test graph from the LEDA tutorial:
   * http://www.leda-tutorial.org/en/unofficial/ch05s03s05.html
   */
  igraph_t graph;
  igraph_vector_bool_t types;
  igraph_vector_long_t matching;
  igraph_integer_t matching_size;
  igraph_bool_t is_matching;
  int i;

  igraph_small(&graph, 0, 0,
      0, 8, 0, 12, 0, 14,
      1, 9, 1, 10, 1, 13,
      2, 8, 2, 9,
      3, 10, 3, 11, 3, 13,
      4, 9, 4, 14,
      5, 14,
      6, 9, 6, 14,
      7, 8, 7, 12, 7, 14
      , -1);
  igraph_vector_bool_init(&types, 15);
  for (i = 0; i < 15; i++)
    VECTOR(types)[i] = (i >= 8);
  igraph_vector_long_init(&matching, 0);

  igraph_i_maximum_bipartite_matching_unweighted(&graph, &types, &matching_size, &matching);
  if (matching_size != 6) {
    printf("matching_size is %ld, expected: 6\n", (long)matching_size);
    return 1;
  }
  igraph_is_maximal_matching(&graph, &types, &matching, &is_matching);
  if (!is_matching) {
    printf("not a matching: ");
    igraph_vector_long_print(&matching);
    return 3;
  }
  igraph_vector_long_print(&matching);
  igraph_vector_long_destroy(&matching);
  igraph_vector_bool_destroy(&types);
  igraph_destroy(&graph);

  return 0;
}
int main() {
  igraph_i_set_attribute_table(&igraph_cattribute_table);

  if (test_graph_from_leda_tutorial())
    return 1;
  if (!IGRAPH_FINALLY_STACK_EMPTY) {
    printf("Finally stack still has %d elements.\n", IGRAPH_FINALLY_STACK_SIZE());
    return 4;
  }

  return 0;
}
