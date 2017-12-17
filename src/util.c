#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "util.h"

static const char *progname;

static void
vwarn(const char *fmt, va_list args)
{
	if (progname) fprintf(stderr, "%s: ", progname);
	vfprintf(stderr, fmt, args);
	putc('\n', stderr);
}

void
warn(const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	vwarn(fmt, args);
	va_end(args);
}

void
die(const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	vwarn(fmt, args);
	va_end(args);
	exit(EXIT_FAILURE);
}

void
bang(const char *fmt, ...)
{
	const char *err = strerror(errno);
	va_list args;
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, ": %s\n", err);
	exit(EXIT_FAILURE);
}

static void
space(int k)
{
	while (k--) putchar(' ');
}

void
print_options(const struct option *opt)
{
	int len, max;
	const struct option *p;
	const char *s;

	max = 0;
	for (p = opt; p->name; ++p) {
		if (*p->description == '*') continue;
		len = strlen(p->name);
		if (p->argname) len += strlen(p->argname) + 1;
		if (max < len) max = len;
	}

	for (p = opt; p->name; ++p) {
		if (p->letter) {
			printf(" -%c,", p->letter);
		}
		else {
			space(4);
		}
		if (p->argname) {
			len = printf(" --%s=%s", p->name, p->argname);
		}
		else {
			len = printf(" --%s", p->name);
		}
		s = p->description;
		if (*s == '*') {
			space(3);
			puts(s+1);
		}
		else {
			space(max - len + 6);
			puts(s);
		}
	}
}

void
parse_arguments(int *argc, char ***argv, const struct option *options)
{
	char *opt, *arg;
	size_t len = 0;
	const struct option *p;

	progname = **argv;

	/* Keep going until we run out of arguments */
	for (--*argc, ++*argv; (opt = **argv) != NULL; --*argc, ++*argv) {

		/* We're only interested in options */
		if (*opt != '-') break;

		/* Check for long options */
		if (*++opt == '-') {
			++opt;

			/* If just --, then terminate processing */
			if (!*opt) {
				--*argc;
				++*argv;
				break;
			}

			/* Otherwise, look for options */
			for (p = options; p->name; ++p) {
				len = strlen(p->name);
				if (!strncmp(p->name, opt, len))
					break;
			}

			/* If we didn't find it, alert the user */
			if (!p->name) {
				arg = strchr(opt, '=');
				if (arg) *arg = '\0';
				die("unknown option --%s", opt);
			}

			/* Parse an argument if required */
			arg = NULL;
			if (p->argname) {
				if (opt[len] == '=' && opt[len+1]) {
					arg = opt + len + 1;
				}
				else {
					die("option --%s takes an argument",
						p->name);
				}
			}
			else {
				/* Otherwise, no trailing characters */
				if (strlen(opt) != len) {
					die("unknown option --%s", opt);
				}
			}

			/* Invoke the callback and move on */
			p->callback(arg);
		}
		/* Short-style options, might be in a cluster */
		else while (*opt) {

			/* Find the corresponding option */
			for (p = options; p->name; ++p) {
				if (*opt == p->letter) break;
			}
			if (!p->name) {
				die("unknown option -%c", *opt);
			}

			/* Compute the argument if required */
			arg = NULL;
			if (p->argname) {
				/* Rest of this opt, or next */
				if (!*++opt) {
					opt = *++*argv;
					--*argc;
				}
				if (!opt) {
					die("option -%c takes an argument",
						p->letter);
				}
				arg = opt;
				opt += strlen(opt);
			}
			else {
				++opt;
			}

			/* Invoke the callback */
			p->callback(arg);
		}
	}
}
