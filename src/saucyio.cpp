/*
 * saucyio.c
 * Input/output routines for saucy
 *
 * by Paul T. Darga <pdarga@umich.edu>
 * and Mark Liffiton <liffiton@umich.edu>
 *
 * Copyright (C) 2004, The Regents of the University of Michigan
 * See the LICENSE file for details.
 */

#include <Rcpp.h>

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "saucy.h"
#include "amorph.h"
#include "util.h"
#include "platform.h"

using namespace Rcpp;

static int
read_int(gzFile f, int *k)
{
	int c, r = 0, neg = 0;
	for (c = gzgetc(f); c != '-' && !isdigit(c); ) {
		for (; isspace(c); c = gzgetc(f));
		for (; c == 'c'; c = gzgetc(f)) {
			while ((c = gzgetc(f)) != '\n') {
				if (c == EOF) return 0;
			}
		}
	}
	if (c == '-') {
		neg = 1;
		c = gzgetc(f);
	}
	if (!isdigit(c)) return 0;
	for (; isdigit(c); c = gzgetc(f)) {
		r = r * 10 + c - '0';
	}
	*k = neg ? -r : r;
	return isspace(c);
}

void
amorph_graph_free(struct amorph_graph *g)
{
	free(g->sg.adj);
	free(g->sg.edg);
	free(g->colors);
	free(g);
}

static int
init_fixadj1(int n, int *adj)
{
	int val, sum, i;

	/* Translate adj values to real locations */
	val = adj[0]; sum = 0; adj[0] = 0;
	for (i = 1; i < n; ++i) {
		sum += val;
		val = adj[i];
		adj[i] = sum;
	}
	return sum + val;
}

static void
init_fixadj2(int n, int e, int *adj)
{
	int i;

	/* Translate again-broken sizes to adj values */
	for (i = n-1; i > 0; --i) {
		adj[i] = adj[i-1];
	}
	adj[0] = 0;
	adj[n] = e;
}

static void
add_edge(int a, int b, int *adj, int *edg)
{
	edg[adj[a]++] = b;
	edg[adj[b]++] = a;
}

static void
amorph_print_automorphism(
	int n, const int *gamma, int nsupp, const int *support,
	struct amorph_graph *g, char *marks)
{
	int i, j, k;

	/* We presume support is already sorted */
	for (i = 0; i < nsupp; ++i) {
		k = support[i];

		/* Skip elements already seen */
		if (marks[k]) continue;

		/* Start an orbit */
		marks[k] = 1;
		Rcout << "(" << k;
		// printf("(%d", k);

		/* Mark and notify elements in this orbit */
		for (j = gamma[k]; j != k; j = gamma[j]) {
			marks[j] = 1;
		  Rcout << " " << j;
			// printf(" %d", j);
		}

		/* Finish off the orbit */
    Rcout << ")" ;
		// putchar(')');
	}
	Rcout << std::endl;
	// putchar('\n');

	/* Clean up after ourselves */
	for (i = 0; i < nsupp; ++i) {
		marks[support[i]] = 0;
	}
}

/* return value >1 indicates error */
static int
dupe_check(int n, int *adj, int *edg)
{
	int i, j, self_loop_ctr;
	int *dupe_tmp = (int *)calloc(n, sizeof(int));
	if (!dupe_tmp) {
		warn("can't allocate memory");
		free(dupe_tmp);
		return 2;
	}

	/* Check outgoing edges of each vertex for duplicate endpoints */
	for (i = 0; i < n; ++i) {
		self_loop_ctr = 0;
		for (j = adj[i] ; j < adj[i+1] ; j++) {
			/* Self-loops lead to two entries of edg[j]==i,
			 * so we check for those and only worry if we see
			 * 3 hits of edg[j]==i (which means 2 self-loops).
			 */
			if (edg[j] == i) {
				++self_loop_ctr;
				if (self_loop_ctr > 2) {
					warn("duplicate edge in input");
					free(dupe_tmp);
					return 1;
				}
			}
			/* If we have recorded this vertex as connected to i,
			 * we have a dupe.
			 * Using i+1 because we used calloc above, and 0 is
			 * a valid vertex index.
			 */
			else if (dupe_tmp[edg[j]] == i+1) {
				warn("duplicate edge in input");
				free(dupe_tmp);
				return 1;
			}
			dupe_tmp[edg[j]] = i+1;
		}
	}

	free(dupe_tmp);
	return 0;
}

struct amorph_graph *
amorph_read(const char *filename, int digraph)
{
	int i, j, k, n, e, p, *aout, *eout, *ain, *ein, *colors;
	struct amorph_graph *g = NULL;
	gzFile f;

	/* Open the file */
	f = gzopen(filename, "r");
	if (f == NULL) goto out;

	/* Read the sizes */
	if (!read_int(f, &n) || !read_int(f, &e) || !read_int(f, &p)) {
		goto out_close;
	}

	/* Allocate everything */
	g = (struct amorph_graph *)malloc(sizeof(struct amorph_graph));
	aout = (int *)calloc(digraph ? (2*n+2) : (n+1), sizeof(int));
	eout = (int *)malloc(2 * e * sizeof(int));
	colors = (int *)malloc(n * sizeof(int));
	if (!g || !aout || !eout || !colors) goto out_free;

	g->sg.n = n;
	g->sg.e = e;
	g->sg.adj = aout;
	g->sg.edg = eout;
	g->colors = colors;

	ain = aout + (digraph ? n+1 : 0);
	ein = eout + (digraph ? e : 0);

	/* Initial coloring with provided splits */
	for (i = j = 0; i < p - 1; ++i) {
		if (!read_int(f, &k)) goto out_free;
		while (j < k) {
			colors[j++] = i;
		}
	}
	while (j < n) {
		colors[j++] = i;
	}

	/* Count the size of each adjacency list */
	for (i = 0; i < e; ++i) {
		if (!read_int(f, &j) || !read_int(f, &k)) goto out_free;
		++aout[j]; ++ain[k];
	}

	/* Fix that */
	init_fixadj1(n, aout);
	if (digraph) init_fixadj1(n, ain);

	/* Go back to the front of the edges */
	gzrewind(f);
	for (i = 0; i < p + 2; ++i) {
		if (!read_int(f, &k)) goto out_free;
	}

	/* Insert adjacencies */
	for (i = 0; i < e; ++i) {
		if (!read_int(f, &j) || !read_int(f, &k)) goto out_free;

		/* Simple input validation: check vertex values */
		if (j >= n || j < 0) {
			warn("invalid vertex in input: %d", j);
			goto out_free;
		}
		if (k >= n || k < 0) {
			warn("invalid vertex in input: %d", k);
			goto out_free;
		}

		eout[aout[j]++] = k;
		ein[ain[k]++] = j;
	}

	/* Fix that too */
	if (digraph) {
		init_fixadj2(n, e, aout);
		init_fixadj2(n, e, ain);
	}
	else {
		init_fixadj2(n, 2 * e, aout);
	}

	/* Check for duplicate edges */
	if (dupe_check(n, aout, eout)) goto out_free;

	/* Assign the functions */
	g->consumer = amorph_print_automorphism;
	g->free = amorph_graph_free;
	g->stats = NULL;
	goto out_close;

out_free:
	free(g);
	free(aout);
	free(eout);
	free(colors);
	g = NULL;
out_close:
	gzclose(f);
out:
	return g;
}

static void
gap_print_automorphism(
	int n, const int *gamma, int nsupp, const int *support,
	struct amorph_graph *g, char *marks)
{
	int i, j, k;

	/* We presume support is already sorted */
	for (i = 0; i < nsupp; ++i) {
		k = support[i];

		/* Skip elements already seen */
		if (marks[k]) continue;

		/* Start an orbit */
		marks[k] = 1;
		Rcout << "(" << k+1;
		// printf("(%d", k+1);

		/* Mark and notify elements in this orbit */
		for (j = gamma[k]; j != k; j = gamma[j]) {
			marks[j] = 1;
		  Rcout << "," << j+1;
			// printf(",%d", j+1);
		}

		/* Finish off the orbit */
		Rcout << ")"; //putchar(')');
	}

	/* Clean up after ourselves */
	for (i = 0; i < nsupp; ++i) {
		marks[support[i]] = 0;
	}
}

struct amorph_graph *
amorph_read_gap(const char *filename)
{
	int c, i, j, k, n, e, *adj, *edg, *colors;
	struct amorph_graph *g = NULL;
	FILE *f;
	fpos_t fpos;

	/* Open the file */
	f = fopen(filename, "r");
	if (f == NULL) goto out;

	/* Skip leading chaff */
	do {
		while ((c = getc(f)) != '[' && c != EOF);
	} while ((c = getc(f)) != '[' && c != EOF);
	if (c == EOF) goto out_close;
	ungetc('[', f);

	/* Remember this spot in the file */
	if (fgetpos(f, &fpos) == -1) goto out_close;

	/* First pass: count the edges */
	e = 0;
	do {
		++e;
		if (fscanf(f, "[%d,%d]", &j, &k) != 2) goto out_close;
	} while ((c = getc(f)) == ',');
	if (c == EOF) goto out_close;

	/* Read the number of vertices */
	if (fscanf(f, ", %d)), [", &n) != 1) goto out_close;

	/* Do the allocation */
	g = (struct amorph_graph *)malloc(sizeof(struct amorph_graph));
	adj = (int *)calloc(n+1, sizeof(int));
	edg = (int *)malloc(2 * e * sizeof(int));
	colors = (int *)malloc(n * sizeof(int));
	if (!g || !adj || !edg || !colors) goto out_free;

	g->sg.n = n;
	g->sg.e = e;
	g->sg.adj = adj;
	g->sg.edg = edg;
	g->colors = colors;

	/* Go back to the front of the edges */
	if (fsetpos(f, &fpos) == -1) goto out_free;

	/* Second pass: count adjacencies for each vertex */
	do {
		if (fscanf(f, "[%d,%d]", &j, &k) != 2) goto out_free;
		++adj[j-1]; ++adj[k-1];
	} while ((c = getc(f)) == ',');
	if (c == EOF) goto out_free;

	/* Fix that */
	init_fixadj1(n, adj);

	/* Go back again */
	if (fsetpos(f, &fpos) == -1) goto out_free;

	/* Third pass: actually insert the edges */
	do {
		if (fscanf(f, "[%d,%d]", &j, &k) != 2) goto out_free;
		add_edge(j-1, k-1, adj, edg);
	} while ((c = getc(f)) == ',');
	if (c == EOF) goto out_free;

	/* Fix that too */
	init_fixadj2(n, 2 * e, adj);

	/* Read the number of vertices (again) */
	if (fscanf(f, ", %d)), [", &n) != 1) goto out_free;

	/* We've read the edges; now read the coloring */
	i = 0; do {

		/* Eat a left bracket */
		if (getc(f) == EOF) goto out_free;

		/* Try eating a right bracket; otherwise keep going */
		if ((c = getc(f)) == ']') continue;
		ungetc(c, f);

		/* Read in the entire color */
		do {
			if (fscanf(f, "%d", &j) != 1) goto out_free;
			colors[j-1] = i;
		} while ((c = getc(f)) == ',');
		if (c == EOF) goto out_free;

		/* Mark the end of the cell */
		++i;
	} while ((c = getc(f)) == ',');
	if (c == EOF) goto out_free;

	/* Check for duplicate edges */
	if (dupe_check(n, adj, edg)) goto out_free;

	/* Assign the functions */
	g->consumer = gap_print_automorphism;
	g->free = amorph_graph_free;
	g->stats = NULL;
	goto out_close;

out_free:
	free(g);
	free(adj);
	free(edg);
	free(colors);
	g = NULL;
out_close:
	fclose(f);
out:
	return g;
}

static int verify(gzFile f, const char *s) {
	for (; *s; ++s) {
		if (gzgetc(f) != *s) return 0;
	}
	return 1;
}

static int dimacs_header(gzFile f, int *v, int *nc) {
	int c;
	gzrewind(f);
	while ((c = gzgetc(f)) == 'c') {
		while ((c = gzgetc(f)) != '\n') {
			if (c == EOF) break;
		}
	}
	if (c != 'p' || !verify(f, " cnf ")) {
		warn("invalid DIMACS header");
		return 0;
	}
	if (!read_int(f, v) || !read_int(f, nc)) {
		warn("invalid DIMACS header");
		return 0;
	}
	return 1;
}

static int
l2v(int k, int v)
{
	return (k > 0 ? k : v - k) - 1;
}

static int
v2l(int k, int v)
{
	return k < v ? k + 1 : v - k - 1;
}

static void
dimacs_graph_free(struct amorph_graph *g)
{
	free(g->data);
	amorph_graph_free(g);
}

static void
dimacs_stats(struct amorph_graph *g, FILE *f)
{
	struct dimacs_info *info = (struct dimacs_info *)g->data;
	fprintf(f, "variables = %d\n", info->vars);
	fprintf(f, "clauses = %d\n", info->orig_clauses);
	fprintf(f, "non-binary clauses = %d\n", info->clauses);
	fprintf(f, "literals = %d\n", info->literals);
}

static void
dimacs_print_automorphism(
	int n, const int *gamma, int nsupp, const int *support,
	struct amorph_graph *g, char *marks)
{
	struct dimacs_info *info = (struct dimacs_info *)g->data;
	int i, j, k, v = info->vars;
	int printed = 0;

	/* We presume support is already sorted */
	for (i = 0; i < nsupp; ++i) {
		k = support[i];
		if (k >= 2*v) break;
		if (marks[k]) continue;

		/* Start an orbit */
		printed = 1;
		marks[k] = 1;
		Rcout << "(" << v2l(k,v);
		// printf("(%d", v2l(k,v));

		/* Mark and notify elements in this orbit */
		for (j = gamma[k]; j != k; j = gamma[j]) {
			marks[j] = 1;
		  Rcout << " " << v2l(j,v);
			// printf(" %d", v2l(j,v));
		}

		/* Finish off the orbit */
		Rcout << ")";
		// putchar(')');
	}
	if (printed) Rcout << std::endl; // putchar('\n');

	/* Clean up after ourselves */
	for (i = 0; i < nsupp; ++i) {
		marks[support[i]] = 0;
	}
}

struct amorph_graph *
amorph_read_dimacs(const char *filename)
{
	struct amorph_graph *g = NULL;
	struct dimacs_info *info = NULL;
	gzFile f;
	int v, c, i, x, y, z, n, e, nc, *adj, *edg, *colors;
	int lits = 0;

	/* Set to NULL to simplify freeing in error path */
	adj = edg = colors = NULL;;

	f = gzopen(filename, "r");
	if (f == NULL) goto out;
	if (!dimacs_header(f, &v, &c)) goto out_close;

	/*
	 * We start out assuming that we will need vertices for every
	 * clause.  As we read the formula the first time, we keep
	 * track of the clause vertices we actually need, with "nc".
	 * We allocate n integers for the adj array no matter what,
	 * though, which may waste a bit of space.  Oh well.
	 */
	n = 2 * v + c;
	adj = (int *)calloc(n+1, sizeof(int));
	if (adj == NULL) goto out_close;

	/* Count boolean consistency edges */
	for (i = 0; i < 2*v; ++i) {
		++adj[i];
	}

	/* Pass 1: Count degrees */
	for (i = nc = 2*v; i < n; ++i) {

		/* Disallow empty clauses */
		if (!read_int(f, &x)) goto out_free;
		if (x == 0) goto out_free;
		++adj[l2v(x,v)];
		++lits;

		/* Treat unary clauses like other clauses */
		if (!read_int(f, &y)) goto out_free;
		if (y == 0) {
			++adj[nc++];
			continue;
		}
		++adj[l2v(y,v)];
		++lits;
		if (!read_int(f, &z)) goto out_free;
		if (z != 0) {
			/*
			 * The clause is not binary; add an adjacency
			 * per variable and the total to the clauses.
			 */
			adj[nc] = 2;
			do {
				++adj[l2v(z,v)];
				++lits;
				++adj[nc];
				if (!read_int(f, &z)) goto out_free;
			} while (z != 0);
			++nc;
		}
	}

	/* Fix up */
	e = init_fixadj1(nc, adj);

	/* Allocate for edges, now that we know how many we'll have */
	edg = (int *)malloc(e * sizeof(int));
	if (edg == NULL) goto out_free;

	/* Rewind */
	if (!dimacs_header(f, &v, &c)) goto out_free;

	/* Populate boolean consistency */
	for (i = 0; i < v; ++i) {
		add_edge(i, i+v, adj, edg);
	}

	/* Pass 2: populate edge array */
	for (i = nc = 2*v; i < n; ++i) {

		/* Disallow empty clauses */
		if (!read_int(f, &x)) goto out_free;
		if (x == 0) goto out_free;
		x = l2v(x,v);

		/* Treat unary clauses like other clauses */
		if (!read_int(f, &y)) goto out_free;
		if (y == 0) {
			add_edge(x, nc, adj, edg);
			nc++;
			continue;
		}
		y = l2v(y,v);
		if (!read_int(f, &z)) goto out_free;
		if (z == 0) {
			add_edge(x, y, adj, edg);
		}
		else {
			add_edge(x, nc, adj, edg);
			add_edge(y, nc, adj, edg);
			do {
				z = l2v(z,v);
				add_edge(z, nc, adj, edg);
				if (!read_int(f, &z)) goto out_free;
			} while (z != 0);
			++nc;
		}
	}

	/* Fix adj again */
	init_fixadj2(nc, e, adj);

	/* Check for duplicate edges */
	if (dupe_check(n, adj, edg)) goto out_free;

	/* Initialize colors; clauses are separate from variables */
	colors = (int *)malloc(nc * sizeof(int));
	if (colors == NULL) goto out_free;
	for (i = 0; i < 2*v; ++i) {
		colors[i] = 0;
	}
	for (i = 2*v; i < nc; ++i) {
		colors[i] = 1;
	}

	/* Allocate the other structures now */
	g = (struct amorph_graph *)malloc(sizeof(struct amorph_graph));
	info = (struct dimacs_info *)malloc(sizeof(struct dimacs_info));
	if (!g || !info) goto out_free;

	/* Graph data */
	g->sg.n = nc;
	g->sg.e = e/2;
	g->sg.adj = adj;
	g->sg.edg = edg;
	g->colors = colors;

	/* Dimacs data */
	info->vars = v;
	info->clauses = nc - 2*v;
	info->literals = lits;
	info->orig_clauses = c;

	/* Amorph functions */
	g->consumer = dimacs_print_automorphism;
	g->data = info;
	g->free = dimacs_graph_free;
	g->stats = dimacs_stats;
	goto out_close;

out_free:
	free(g);
	free(colors);
	free(adj);
	free(edg);
	free(info);
	g = NULL;
out_close:
	gzclose(f);
out:
	return g;
}
