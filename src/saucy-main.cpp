#include <Rcpp.h>

using namespace Rcpp;

/*
 * main.c
 * Program entry and option parsing
 *
 * by Paul T. Darga <pdarga@umich.edu>
 * and Mark Liffiton <liffiton@umich.edu>
 * and Hadi Katebi <hadik@eecs.umich.edu>
 *
 * Copyright (C) 2004, The Regents of the University of Michigan
 * See the LICENSE file for details.
 */

#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include "saucy.h"
#include "amorph.h"
#include "util.h"
#include "platform.h"

// static char *filename;   /* Graph file we're reading */
// static int timeout = 0;      /* Seconds before quitting after refinement */
static sig_atomic_t timeout_flag = 0; /* Has the alarm gone off yet? */
static int stats_mode;  /* Print out stats when we're done */
static int quiet_mode;  /* Don't output automorphisms */
static int gap_mode;    /* Do GAP I/O (interface with Shatter) */
static int cnf_mode;   /* Read CNF instead of graphs */
static int digraph_mode; /* Read digraphs; order matters in input */
// static int repeat = 1; /* Repeat count, for benchmarking */
static int first;      /* Have we seen the first automorphism? (for gap) */
static char *marks;    /* "Bit" vector for printing */

/* Stats are global so we can print them from the signal handler */
struct saucy_stats sstats;

/* Do nothing but set a flag to be tested during the search */
static void timeout_handler(void) {
  timeout_flag = 1;
}

static int on_automorphism(int n, const int *gamma, int k, int *support, void *arg) {
  struct amorph_graph *g = (struct amorph_graph *)arg;
  if (!quiet_mode) {
    qsort_integers(support, k);
    if (gap_mode) {
      putchar(!first ? '[' : ',');
      putchar('\n');
      first = 1;
    }
    g->consumer(n, gamma, k, support, g, marks);
  }
  return !timeout_flag;
}

// ----- HACK -----------
// --cnf, --digraph, and --shatter are mutually exclusive"

// [[Rcpp::export]]
Rcpp::List saucy_int(std::string filename, std::string mode, int timeout, int rpt) {

  // Rcout << "filename " << filename << " mode " << mode << " os " << output_stats <<
  //   " q " << quiet << " to " << timeout << " rpt " << rpt << std::endl;

  int repeat = rpt;

  quiet_mode = 0;
  stats_mode = 0;

  cnf_mode = mode == "cnf" ? 1 : 0;
  digraph_mode = mode == "digraph" ? 1 : 0;
  gap_mode = mode == "shatter" ? 1 : 0;

  struct saucy *s;
  struct amorph_graph *g = NULL;
  long cpu_time;
  int i, n;

  /* Read the input file */
  if (gap_mode) {
    g = amorph_read_gap(filename.c_str());
  } else if (cnf_mode) {
    g = amorph_read_dimacs(filename.c_str());
  } else {
    g = amorph_read(filename.c_str(), digraph_mode);
  }

  if (!g) {
    Rf_warning("unable to read input file");
    return(NULL);
  }

  n = g->sg.n;

  /* Allocate some memory to facilitate printing */
  marks = (char *)calloc(n, sizeof(char));
  if (!marks) {
    Rf_warning("out of memory");
    return(NULL);
  }

  /* Allocate saucy space */
  s = saucy_alloc(n);
  if (s == NULL) {
    Rf_warning("saucy initialization failed");
    return(NULL);
  }

  /* Set up the alarm for timeouts */
  if (timeout > 0) platform_set_timer(timeout, timeout_handler);

  /* Print statistics when signaled */
  // platform_set_user_signal(stats_handler);

  /* Start timing */
  cpu_time = platform_clock();

  /* Run the search */
  for (i = 0; i < repeat; ++i) {
    saucy_search(s, &g->sg, digraph_mode, g->colors, on_automorphism, g, &sstats);
  }

  /* Finish timing */
  cpu_time = platform_clock() - cpu_time;

  /* Warn if timeout */
  if (timeout_flag) Rf_warning("search timed out");

  Rcpp::List ret = Rcpp::List::create(
    _["input_file"] = filename,
    _["vertices"] = n,
    _["edges"] = g->sg.e,
    _["group_size_base"] =  sstats.grpsize_base,
    _["group_size_exp"] = sstats.grpsize_exp,
    _["levels"] = sstats.levels,
    _["nodes"] = sstats.nodes,
    _["generators"] = sstats.gens,
    _["total_support"] = sstats.support,
    _["average_support"] = divide(sstats.support, sstats.gens),
    _["nodes_per_generator"] = divide(sstats.nodes, sstats.gens),
    _["bad_nodes"] = sstats.bads
  );

  /* Cleanup */
  saucy_free(s);
  g->free(g);
  free(marks);

  return(ret);

}
