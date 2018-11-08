#ifndef _WSQ_H
#define _WSQ_H

//#ifndef _STDIO_H
//#include <stdio.h>
//#endif
///*
//#ifndef _IHEAD_H
//#include "ihead.h"
//#endif
//*/
//#ifndef _JPEGL_H
#include "jpegl.h"
//#endif
//
//#ifndef TRUE
//#define TRUE  1
//#define FALSE 0
//#endif

/* WSQ Marker Definitions */
#define SOI_WSQ 0xffa0
#define EOI_WSQ 0xffa1
#define SOF_WSQ 0xffa2
#define SOB_WSQ 0xffa3
#define DTT_WSQ 0xffa4
#define DQT_WSQ 0xffa5
#define DHT_WSQ 0xffa6
#define DRT_WSQ 0xffa7
#define COM_WSQ 0xffa8
/* Case for getting ANY marker. */
#define ANY_WSQ 0xffff
#define TBLS_N_SOB   (TBLS_N_SOF + 2)

/* Filter Bank Definitions */
#ifdef FILTBANK_EVEN_8X8_1
#define MAX_HIFILT   8
#define MAX_LOFILT   8
#else
#define MAX_HIFILT   7
#define MAX_LOFILT   9
#endif

/* Subband Definitions */
#define STRT_SUBBAND_2      19
#define STRT_SUBBAND_3      52
#define MAX_SUBBANDS        64
#define NUM_SUBBANDS        60
#define STRT_SUBBAND_DEL    (NUM_SUBBANDS)
#define STRT_SIZE_REGION_2  4
#define STRT_SIZE_REGION_3  51

#define MIN_IMG_DIM         256

#define WHITE               255
#define BLACK               0

#define COEFF_CODE          0
#define RUN_CODE            1

#define RAW_IMAGE           1
#define IHEAD_IMAGE         0

#define VARIANCE_THRESH     1.01

typedef struct quantization {
   float q;  /* quantization level */
   float cr; /* compression ratio */
   float r;  /* compression bitrate */
//   float qbss_t[MAX_SUBBANDS];
   double qbss[MAX_SUBBANDS];
   double qzbs[MAX_SUBBANDS];
   double var[MAX_SUBBANDS];
   //float qbss[MAX_SUBBANDS];
   //float qzbs[MAX_SUBBANDS];
   //float var[MAX_SUBBANDS];
} QUANT_VALS;

typedef struct wavlet_tree {
   int x;
   int y;
   int lenx;
   int leny;
   int inv_rw;
   int inv_cl;
} W_TREE;
#define W_TREELEN 20

typedef struct quant_tree {
   short x;	/* UL corner of block */
   short y;
   short lenx;  /* block size */
   short leny;  /* block size */
} Q_TREE;
#define Q_TREELEN 64

typedef struct table_dtt {
   double *lofilt;
   double *hifilt;
   //float *lofilt;
   //float *hifilt;
   unsigned char losz;
   unsigned char hisz;
   char lodef;
   char hidef;
} DTT_TABLE;

typedef struct table_dqt {
   double bin_center;
   double q_bin[MAX_SUBBANDS];
   double z_bin[MAX_SUBBANDS];
   //float bin_center;
   //float q_bin[MAX_SUBBANDS];
   //float z_bin[MAX_SUBBANDS];
   char dqt_def;
} DQT_TABLE;

#define MAX_DHT_TABLES  8

/* Defined in jpegl.h */
// #define MAX_HUFFBITS      16  //DO NOT CHANGE THIS CONSTANT!! 
#define MAX_HUFFCOUNTS_WSQ  256  /* Length of code table: change as needed */
                                 /* but DO NOT EXCEED 256 */
#define MAX_HUFFCOEFF        74  /* -73 .. +74 */
#define MAX_HUFFZRUN        100

#define MAX_HUFFCOEFF1        45
#define MAX_HUFFZRUN1         30
#define MAX_HUFFCOEFF23       10
#define MAX_HUFFZRUN23        75


typedef struct table_dht {
   unsigned char tabdef;
   unsigned char huffbits[MAX_HUFFBITS];
   unsigned char huffvalues[MAX_HUFFCOUNTS_WSQ+1];
} DHT_TABLE;

typedef struct header_frm {
   unsigned char black;
   unsigned char white;
   unsigned short width;
   unsigned short height;
   float m_shift; 
   float r_scale;
   unsigned char wsq_encoder;
   unsigned short software;
} FRM_HEADER_WSQ;

#endif /* !_WSQ_H */
