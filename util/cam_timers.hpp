#ifdef CAM_TIMERS
#include <gptl.h>

#define TIMER_START(X) GPTLstart(X)
#define TIMER_STOP(X) GPTLstop(X)
#else
#define TIMER_START(X)
#define TIMER_STOP(X)
#endif

