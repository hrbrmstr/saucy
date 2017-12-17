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

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <errno.h>
#include "saucy.h"
#include "amorph.h"
#include "util.h"
#include "platform.h"

static const char *sbpfile;
static int *p;
static int *supp;
static FILE *sbp;
static int orig_vars;
static int vars;
static int clauses;
static int literals;
static char *marks;
static int stats_mode;
static int quiet_mode;
static int violations;
static long shatter_time;

static void arg_sbpfile(char *arg) { sbpfile = arg; }
static void arg_stats(char *arg) { stats_mode = 1; }
static void arg_quiet(char *arg) { quiet_mode = 1; }

static void arg_help(char *arg);

static void arg_version(char *arg)
{
	fprintf(stderr, "shatter (saucy) %s\n", SAUCY_VERSION);
	exit(0);
}

static struct option options[] = {
	{ "sbpfile", 'o', "FILE", arg_sbpfile,
	"put symmetry breaking predicates in FILE" },
	{ "stats", 's', 0, arg_stats, "print statistics after execution" },
	{ "quiet", 'q', 0, arg_quiet,
	"don't output final CNF formula (for use with -s or -o)" },
	{ "help", 0, 0, arg_help, "print this help message" },
	{ "version", 0, 0, arg_version, "version information" },
	{ 0, 0, 0, 0, 0 }
};

static void
arg_help(char *arg)
{
	fprintf(stderr, "usage: shatter [OPTION...] FILE\n");
	print_options(options);
	exit(0);
}

static int
name(int k)
{
	return k >= orig_vars ? k - orig_vars : k;
}

static int
negate(int k)
{
	return k >= orig_vars ? k - orig_vars : k + orig_vars;
}

static void
clause(int x, ...)
{
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

static int
shatter(int n, const int *perm, int nsupp, int *support, void *arg)
{
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

static int
time_shatter(int n, const int *perm, int nsupp, int *support, void *arg)
{
	long cpu_time = platform_clock();
	int ret = shatter(n, perm, nsupp, support, arg);
	shatter_time += platform_clock() - cpu_time;
	return ret;
}

static void
eat_line(FILE *f)
{
	int c;
	while ((c = getc(f)) != '\n') {
		if (c == EOF) die("unexpected end of file");
	}
}

static void
print_file(FILE *f)
{
	char buf[4096];
	int c;
	while ((c = fread(buf, 1, sizeof(buf), f))) {
		fwrite(buf, 1, c, stdout);
	}
	if (ferror(f)) die("error reading file");
	fclose(f);
}

int entry_main(int argc, char **argv) {
	const char *filename;
	struct amorph_graph *g;
	struct saucy *s;
	struct saucy_stats stats;
	struct dimacs_info *info;
	FILE *f;
	int c, n, orig_clauses;
	long cpu_time;

	parse_arguments(&argc, &argv, options);
	if (argc > 1) die("trailing arguments");
	if (argc < 1) die("missing filename");
	filename = *argv;

	g = amorph_read_dimacs(filename);
	if (!g) die("unable to read CNF input file");
	info = g->data;

	n = g->sg.n;
	vars = orig_vars = info->vars;
	clauses = orig_clauses = info->orig_clauses;

	supp = malloc(vars * sizeof(int));
	p = malloc((vars+1) * sizeof(int));
	marks = calloc(vars, sizeof(char));
	if (!supp || !p || !marks) bang("can't allocate memory");

	sbp = sbpfile ? fopen(sbpfile, "w+") : tmpfile();
	if (!sbp) bang("can't create SBP file");

	s = saucy_alloc(n);
	if (!s) die("unable to initialize saucy");

	cpu_time = platform_clock();
	saucy_search(s, &g->sg, 0, g->colors, time_shatter, 0, &stats);
	cpu_time = platform_clock() - cpu_time;

	saucy_free(s);

	if (!quiet_mode) {
		f = fopen(filename, "r");
		if (!f) bang("unable to reopen CNF file");

		while ((c = getc(f)) == 'c') {
			eat_line(f);
		}
		if (c != 'p') die("can't read CNF header");
		eat_line(f);

		printf("p cnf %d %d\n", vars, clauses);
		print_file(f);

		errno = 0;
		rewind(sbp);
		if (errno != 0) bang("rewinding SBP file failed");

		print_file(sbp);
	}

	if (stats_mode) {
		f = quiet_mode ? stdout : stderr;
		fprintf(f, "----------- formula info ----------\n");
		fprintf(f, "input file = %s\n", filename);
		g->stats(g, f);
		fprintf(f, "-------- symmetry discovery -------\n");
		fprintf(f, "vertices = %d\n", g->sg.n);
		fprintf(f, "edges = %d\n", g->sg.e);
		fprintf(f, "group size = %fe%d\n",
			stats.grpsize_base, stats.grpsize_exp);
		fprintf(f, "nodes = %d\n", stats.nodes);
		fprintf(f, "generators = %d\n", stats.gens);
		fprintf(f, "bad nodes = %d\n", stats.bads);
		fprintf(f, "discovery time (s) = %.2f\n",
			divide(cpu_time - shatter_time,
				PLATFORM_CLOCKS_PER_SEC));
		fprintf(f, "----------- shatter info ----------\n");
		fprintf(f, "symmetry breaking clauses = %d\n",
			clauses - orig_clauses);
		fprintf(f, "additional variables = %d\n",
			vars - orig_vars);
		fprintf(f, "additional literals = %d\n", literals);
		fprintf(f, "consistency violations = %d\n", violations);
		fprintf(f, "SBP generation time (s) = %.2f\n",
			divide(shatter_time, PLATFORM_CLOCKS_PER_SEC));
		fprintf(f, "total time (s) = %.2f\n",
			divide(cpu_time, PLATFORM_CLOCKS_PER_SEC));
	}

	g->free(g);
	return 0;
}
