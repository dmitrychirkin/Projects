#ifndef _DEFS_H
#define _DEFS_H

/*********************************************************************/
/*          General Purpose Defines                                  */
/*********************************************************************/
#ifndef True
#define True		1
#define False		0
#endif
#ifndef TRUE
#define TRUE		True
#define FALSE		False
#endif
#define Yes		True
#define No		False
#define Empty		NULL
#ifndef None
#define None		-1
#endif
#ifndef FOUND
#define FOUND            1
#endif
#define NOT_FOUND_NEG   -1
#define EOL		EOF
#ifndef DEG2RAD
#define DEG2RAD	(double)(57.29578)
#endif
#define sround(x) ((int) (((x)<0) ? (x)-0.5 : (x)+0.5))
#define sround_uint(x) ((unsigned int) (((x)<0) ? (x)-0.5 : (x)+0.5))
#ifdef _WINDOWS
#define xor(a, b)  (!(a && b) && (a || b))
#endif
#define align_to_16(_v_)   ((((_v_)+15)>>4)<<4)
#define align_to_32(_v_) ((((_v_)+31)>>5)<<5)
#ifndef CHUNKS
#define CHUNKS          100
#endif

#endif /* !_DEFS_H */
