/*
 * platform.h
 * Platform dependencies
 * by Paul T. Darga
 */

#ifndef SAUCY_PLATFORM_H
#define SAUCY_PLATFORM_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef __GNUC__
#define __attribute__(x)
#endif

#define _inline static __attribute__((always_inline,unused))
#define __internal static __attribute__((unused))  /* "_internal" clashes with def in cygwin's signal.h */

typedef void platform_callback_t(void);
__internal platform_callback_t *__platform_timer_callback;
__internal platform_callback_t *__platform_user_signal_callback;
__internal long __platform_timer_secs;

#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#define gzFile FILE *
#define gzopen fopen
#define gzclose fclose
#define gzgetc getc
#define gzrewind rewind

#define PLATFORM_CLOCKS_PER_SEC 1000

_inline long platform_clock(void) {
	FILETIME create, exit, kernel, user;
	GetProcessTimes(GetCurrentProcess(), &create, &exit, &kernel, &user);

	/* Ignore the high-order bits and hope that's okay */
	return kernel.dwLowDateTime + user.dwLowDateTime;
}

static DWORD WINAPI __platform_timer_wrapper(LPVOID arg) {
	MSG msg;
	SetTimer(NULL, 0, __platform_timer_secs, NULL);
	GetMessage(&msg, NULL, WM_TIMER, WM_TIMER);

	/* The callback gets invoked by the secondary thread */
	__platform_timer_callback();
	return 0;
}

_inline void platform_set_timer(int seconds, platform_callback_t *f) {
	__platform_timer_secs = 1000 * seconds;
	__platform_timer_callback = f;

	/*
	 * Sadly, we need a message pump to get Windows to invoke our
	 * callback.  So we create another thread that lives in the
	 * background and does the message pumping.
	 */
	CreateThread(NULL, 0, __platform_timer_wrapper, NULL, 0, NULL);
}

_inline void platform_set_user_signal(platform_callback_t *f) { }

#else

#include <signal.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>

#define PLATFORM_CLOCKS_PER_SEC CLOCKS_PER_SEC

_inline long platform_clock(void) {
	return clock();
}

static void __platform_timer_wrapper(int sig) {
	__platform_timer_callback();
}

_inline void platform_set_timer(int seconds, platform_callback_t *f) {
	__platform_timer_callback = f;
	signal(SIGALRM, __platform_timer_wrapper);
	alarm(seconds);
}

static void __platform_user_signal_wrapper(int sig) {
	__platform_user_signal_callback();
	signal(SIGUSR1, __platform_user_signal_wrapper);
}

_inline void platform_set_user_signal(platform_callback_t *f) {
	__platform_user_signal_callback = f;
	signal(SIGUSR1, __platform_user_signal_wrapper);
}

#endif

#undef _inline
#undef __internal

#ifdef __cplusplus
}
#endif

#endif
