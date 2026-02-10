#ifndef IQMATH_COMPAT_H
#define IQMATH_COMPAT_H

/*
 * 统一的数值层：
 * - MCU 侧若可用 IQmathLib.h，则用 IQmath 完成核心乘除/开方；
 * - PC 侧自动退回 float，保持可编译调试。
 */
#if defined(__has_include)
#  if __has_include("IQmathLib.h")
#    include "IQmathLib.h"
#    define SIM_USE_VENDOR_IQMATH 1
#  else
#    define SIM_USE_VENDOR_IQMATH 0
#  endif
#else
#  define SIM_USE_VENDOR_IQMATH 0
#endif

#include <math.h>

typedef float sim_iq_t;

#if SIM_USE_VENDOR_IQMATH

static inline sim_iq_t sim_iq_mul(sim_iq_t a, sim_iq_t b) {
    return _IQtoF(_IQmpy(_IQ(a), _IQ(b)));
}

static inline sim_iq_t sim_iq_div(sim_iq_t a, sim_iq_t b) {
    return _IQtoF(_IQdiv(_IQ(a), _IQ(b)));
}

static inline sim_iq_t sim_iq_sqrt(sim_iq_t v) {
    return _IQtoF(_IQsqrt(_IQ(v)));
}

static inline int sim_iq_floor_to_int(sim_iq_t v) {
    return (int)_IQint(_IQ(v));
}

#else

static inline sim_iq_t sim_iq_mul(sim_iq_t a, sim_iq_t b) { return a * b; }
static inline sim_iq_t sim_iq_div(sim_iq_t a, sim_iq_t b) { return a / b; }
static inline sim_iq_t sim_iq_sqrt(sim_iq_t v) { return sqrtf(v); }
static inline int sim_iq_floor_to_int(sim_iq_t v) { return (int)floorf(v); }

#endif

#define SIM_IQ(v) ((sim_iq_t)(v))
#define SIM_TO_FLOAT(v) ((float)(v))
#define SIM_ADD(a,b) ((a) + (b))
#define SIM_SUB(a,b) ((a) - (b))
#define SIM_MUL(a,b) sim_iq_mul((a),(b))
#define SIM_DIV(a,b) sim_iq_div((a),(b))
#define SIM_FLOOR_TO_INT(v) sim_iq_floor_to_int((v))
#define SIM_SQRT(v) sim_iq_sqrt((v))

#endif
