/*
 * shatter.c
 * A reimplementation of Shatter by Fadi A. Aloul
 *
 * by Paul T. Darga
 *
 * This implementation encodes CNF formulas as graphs, where binary
 * clauses are represented as a single edge between literals.  Only
 * the symmetries directly found by saucy are used to form
 * symmetry-breaking predicates (that is, no powers or other
 * enumeration).
 */
#include <Rcpp.h>

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <errno.h>
#include "saucy.h"
#include "amorph.h"
#include "util.h"
#include "platform.h"

using namespace Rcpp;

static const char *sbpfile;
static int *p;
static int *supp;
static FILE *sbp;
static int orig_vars;
static int vars;
static int clauses;
static int literals;
static char *marks;
static int violations;
static long shatter_time;

static int name(int k) {
	return k >= orig_vars ? k - orig_vars : k;
}

static int negate(int k) {
	return k >= orig_vars ? k - orig_vars : k + orig_vars;
}

static void clause(int x, ...) {
	va_list args;
	va_start(args, x);

	while (x != 0) {
		++literals;
		fprintf(sbp, "%d ", x);
		x = va_arg(args, int);
	}
	fputs("0\n", sbp);

	va_end(args);
	++clauses;
}

static int shatter(int n, const int *perm, int nsupp, int *support, void *arg) {
	int i, ns, j, k, x, z, big;

	/* Boolean consistency check */
	for (i = 0; i < nsupp; ++i) {
		k = support[i];
		if (k >= 2 * orig_vars) continue;
		if (negate(perm[k]) != perm[negate(k)]) {
			++violations;
			return 1;
		}
	}

	/*
	 * We only care about positive literals.  If a positive
	 * literal is mapped, so is the negative one.  And we don't
	 * care at all about clauses.  So put just the positive
	 * literals in the "support" before we start.
	 *
	 * Additionally, we don't want to print clauses for the
	 * largest (lex) variable in each orbit.  So we walk each orbit,
	 * and don't include elements that are largest.
	 */
	ns = 0;
	for (i = 0; i < nsupp; ++i) {
		if (support[i] >= 2 * orig_vars) continue;

		k = name(support[i]);
		if (marks[k]) continue;
		marks[k] = 1;

		if (k == name(perm[k])) {
			supp[ns++] = k + 1;
			continue;
		}

		big = k;
		for (j = name(perm[k]); j != k; j = name(perm[j])) {
			marks[j] = 1;
			if (big < j) big = j;
		}

		k = name(support[i]);
		if (k != big) supp[ns++] = k + 1;
		for (j = name(perm[k]); j != k; j = name(perm[j])) {
			if (j != big) supp[ns++] = j + 1;
		}
	}

	/* Do nothing on clause-only symmetries */
	if (ns == 0) return 1;

	qsort_integers(supp, ns);

	for (i = 0; i < nsupp; ++i) {
		if (support[i] < 2 * orig_vars) {
			marks[name(support[i])] = 0;
		}
	}

	/*
	 * Build secondary mapping array in the domain of DIMACS
	 * literals rather than the graph vertices.  This simplifies
	 * the printing later on.
	 */
	for (i = 0; i < ns; ++i) {
		k = supp[i];
		x = perm[k-1];
		p[k] = x < orig_vars ? x + 1 : orig_vars - x - 1;
	}

	z = supp[0];

	/* short circuit simple phase shifts */
	if (p[z] == -z) {
		clause(-z, 0);
		return 1;
	}

	clause(-z, p[z], 0);

	++vars;
	clause(vars, 0);

	for (i = 1; i < ns; ++i) {
		x = supp[i];

		/* again, terminate at phase shift */
		if (p[x] == -x) {
			clause(-vars, -z, -x, 0);
			clause(-vars, p[z], -x, 0);
			break;
		}

		clause(-vars, -z, -x, p[x], 0);
		clause(-vars, -z, vars+1, 0);
		clause(-vars, p[z], -x, p[x], 0);
		clause(-vars, p[z], vars+1, 0);

		++vars;
		z = x;
	}

	return 1;
}

static int time_shatter(int n, const int *perm, int nsupp, int *support, void *arg) {
	long cpu_time = platform_clock();
	int ret = shatter(n, perm, nsupp, support, arg);
	shatter_time += platform_clock() - cpu_time;
	return ret;
}


// [[Rcpp::export]]
List shatter_int(std::string filename, std::string sbp_file, std::string gfile) {

	struct amorph_graph *g;
	struct saucy *s;
	struct saucy_stats stats;
	struct dimacs_info *info;
	FILE *f = NULL;
	int n, orig_clauses;
	long cpu_time;

	sbpfile = sbp_file.c_str();

	g = amorph_read_dimacs(filename.c_str());
	if (!g) {
	  Rf_warning("unable to read CNF input file");
	  return(NULL);
	}

	info = (struct dimacs_info *)g->data;

	n = g->sg.n;
	vars = orig_vars = info->vars;
	clauses = orig_clauses = info->orig_clauses;

	supp = (int *)malloc(vars * sizeof(int));
	p = (int *)malloc((vars+1) * sizeof(int));
	marks = (char *)calloc(vars, sizeof(char));
	if (!supp || !p || !marks) {
	  Rf_warning("can't allocate memory");
	  return(NULL);
	}

	f = fopen(gfile.c_str(), "w+");

	sbp = sbpfile ? fopen(sbpfile, "w+") : tmpfile();
	if (!sbp) {
	  Rf_warning("can't create SBP file");
	  return(NULL);
	}

	s = saucy_alloc(n);
	if (!s) {
	  Rf_warning("unable to initialize saucy");
	  return(NULL);
	}

	cpu_time = platform_clock();

	saucy_search(s, &g->sg, 0, g->colors, time_shatter, 0, &stats);
	cpu_time = platform_clock() - cpu_time;

	saucy_free(s);

	g->stats(g, f);
	fclose(f);

	Rcpp::List ret = Rcpp::List::create(
    _["input_file"] = filename,
    _["vertices"] = g->sg.n,
    _["edges"] = g->sg.e,
    _["group_size_base"] = stats.grpsize_base,
    _["group_size_exp"] =  stats.grpsize_exp,
    _["nodes"] = stats.nodes,
    _["generators"] = stats.gens,
    _["bad_nodes"] = stats.bads,
    _["discovery_time"] = divide(cpu_time - shatter_time, PLATFORM_CLOCKS_PER_SEC),
    _["symmetry_breaking_clauses"] = clauses - orig_clauses,
    _["additional_variables"] = vars - orig_vars,
    _["additional_literals"] = literals,
    _["consistency_violations"] = violations,
    _["sbp_generation_time"] = divide(shatter_time, PLATFORM_CLOCKS_PER_SEC),
    _["total_time"] = divide(cpu_time, PLATFORM_CLOCKS_PER_SEC)
	);

	g->free(g);

	return(ret);

}
