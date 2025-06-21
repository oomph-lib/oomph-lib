//====================================================================
// Macs don't have fpu_control.h. Bugger. We're using this 
// empty replacement if it can't found... Suggestion comes from
//
// http://stackoverflow.com/questions/4766801/compiling-c-floating-point-arithmetic-on-osx-for-shewchuks-triangle-program
//
//====================================================================
#ifndef _FPU_CONTROL_H
#define _FPU_CONTROL_H	1

#warning("Using dummy fpu_control.h header. Triangle may struggle...")

#define _FPU_SETCW(cw) // nothing
#define _FPU_GETCW(cw) // nothing


#endif

