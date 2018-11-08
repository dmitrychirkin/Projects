void conv_img_2_flt(
   float *fip,         /* output float image data  */
   float *m_shift,     /* shifting parameter       */
   float *r_scale,     /* scaling parameter        */
   unsigned char *data,        /* input unsigned char data */
   const int num_pix);  /* num pixels in image      */
void conv_img_2_uchar(
   unsigned char *data,                   /* uchar image pointer    */
   float *img,                    /* image pointer          */
   const int width,               /* image width            */
   const int height,              /* image height           */
   const float m_shift,           /* shifting parameter     */
   const float r_scale);           /* scaling parameter      */
void variance(
   QUANT_VALS *quant_vals, /* quantization parameters */
   Q_TREE q_tree[],        /* quantization "tree"     */
   const int q_treelen,    /* length of q_tree        */
   float *fip,             /* image pointer           */
   const int width,        /* image width             */
   const int height);       /* image height            */
int quantize(
   short **osip,           /* quantized output             */
   int *ocmp_siz,          /* size of quantized output     */
   QUANT_VALS *quant_vals, /* quantization parameters      */
   Q_TREE q_tree[],        /* quantization "tree"          */
   const int q_treelen,    /* size of q_tree               */
   float *fip,             /* floating point image pointer */
   const int width,        /* image width                  */
   const int height);       /* image height                 */
void quant_block_sizes(int *oqsize1, int *oqsize2, int *oqsize3,
                 QUANT_VALS *quant_vals,
                 W_TREE w_tree[], const int w_treelen,
                 Q_TREE q_tree[], const int q_treelen);
int unquantize(
   float **ofip,         /* floating point image pointer         */
   const DQT_TABLE *dqt_table, /* quantization table structure   */
   Q_TREE q_tree[],      /* quantization table structure         */
   const int q_treelen,  /* size of q_tree                       */
   short *sip,           /* quantized image pointer              */
   const int width,      /* image width                          */
   const int height);     /* image height                         */
int wsq_decompose(float *fdata, const int width, const int height,
                  W_TREE w_tree[], const int w_treelen,
                  double *hifilt, const int hisz, double *lofilt, const int losz);
//                  float *hifilt, const int hisz, float *lofilt, const int losz);
int wsq_reconstruct(float *fdata, const int width, const int height,
                  W_TREE w_tree[], const int w_treelen,
                  const DTT_TABLE *dtt_table);
int int_sign(
   const int power);  /* "sign" power */
int image_size(
   const int blocklen,  /* length of the compressed blocks */
   short *huffbits1,    /* huffman table parameters */
   short *huffbits2);

void get_lets(
   float *ew,     /* image pointers for creating subband splits */
   float *old,
   const int len1,       /* temporary length parameters */
   const int len2,
   const int pitch,      /* pitch gives next row_col to filter */
   const int  stride,    /*           stride gives next pixel to filter */
   float *hi,
   const int hsz,   /* NEW */
   float *lo,      /* filter coefficients */
   const int lsz,   /* NEW */
   const int inv);        /* spectral inversion? */

