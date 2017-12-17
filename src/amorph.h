#ifndef SAUCY_AMORPH_H
#define SAUCY_AMORPH_H

#include "saucy.h"

#ifdef __cplusplus
extern "C" {
#endif

struct amorph_graph {
	struct saucy_graph sg;
	int *colors;
	void *data;
	void (*consumer)(int, const int *, int, const int *,
		struct amorph_graph *g, char *);
	void (*free)(struct amorph_graph *);
	void (*stats)(struct amorph_graph *, FILE *f);
};

struct dimacs_info {
	int vars;
	int clauses;
	int literals;
	int orig_clauses;
};

struct amorph_graph *amorph_read(const char *filename, int digraph);
struct amorph_graph *amorph_read_gap(const char *filename);
struct amorph_graph *amorph_read_dimacs(const char *filename);

#ifdef __cplusplus
}
#endif

#endif
