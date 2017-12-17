#ifndef UTIL_H
#define UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef __GNUC__
#define __attribute__(x)
#else
#define _inline static __attribute__((always_inline,unused))
#endif

void warn(const char *fmt, ...) __attribute__((format(printf,1,2)));
void die(const char *fmt, ...) __attribute__((noreturn, format(printf,1,2)));
void bang(const char *fmt, ...) __attribute__((noreturn, format(printf,1,2)));

struct option {
	char *name;
	char letter;
	char *argname;
	void (*callback)(char *);
	char *description;
};

void print_options(const struct option *options);
void parse_arguments(int *argc, char ***argv, const struct option *options);

_inline int
integer_compare(const void *a, const void *b)
{
	const int *aa = (const int *)a, *bb = (const int *)b;
	return *aa < *bb ? -1 : *aa == *bb ? 0 : 1;
}

_inline void
qsort_integers(int *a, int n)
{
	qsort(a, n, sizeof(int), integer_compare);
}

_inline double
divide(int a, int b)
{
	return ((double) a) / ((double) b);
}

#ifdef __cplusplus
}
#endif

#undef _inline
#endif
