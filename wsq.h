#ifndef WSQ_31_H_
#define WSQ_31_H_

// Wsq.h

#pragma once

//namespace WSQ {

	class Wsq
	{
		FILE *stdf;
		int debug;
		double *hifilt;//[MAX_HIFILT];
		double *lofilt;//[MAX_LOFILT];
		unsigned char code;   /*next byte of data*/
		unsigned char code2;  /*stuffed byte of data*/
		unsigned char bit_mask[9]; // = {0x00,0x01,0x03,0x07,0x0f,0x1f,0x3f,0x7f,0xff};
		FRM_HEADER_WSQ frm_header_wsq;
		W_TREE w_tree[W_TREELEN];
		Q_TREE q_tree[Q_TREELEN];
	public:
		Wsq()
		{
#ifdef _DEBUG
			debug = 4; 
#else
			debug = 0; 
#endif
			code = 0;
			code2 = 0;
//			bit_mask[9] = {0x00,0x01,0x03,0x07,0x0f,0x1f,0x3f,0x7f,0xff};
			bit_mask[0] = 0x00;
			bit_mask[1] = 0x01;
			bit_mask[2] = 0x03;
			bit_mask[3] = 0x07;
			bit_mask[4] = 0x0f;
			bit_mask[5] = 0x1f;
			bit_mask[6] = 0x3f;
			bit_mask[7] = 0x7f;  
			bit_mask[8] = 0xff;
			hifilt = new double [MAX_HIFILT];
			lofilt = new double [MAX_LOFILT];
#ifdef FILTBANK_EVEN_8X8_1
			hifilt[0] =  0.03226944131446922,
			hifilt[1] = -0.05261415011924844;
			hifilt[2] = -0.18870142780632693;
			hifilt[3] =  0.60328894481393847;
			hifilt[4] = -0.60328894481393847;
			hifilt[5] =  0.18870142780632693;
			hifilt[6] =  0.05261415011924844;
			hifilt[7] = -0.03226944131446922;

			lofilt[0] =  0.07565691101399093;
			lofilt[1] = -0.12335584105275092;
			lofilt[2] = -0.09789296778409587;
			lofilt[3] =  0.85269867900940344;
			lofilt[4] =  0.85269867900940344;
			lofilt[5] = -0.09789296778409587;
			lofilt[6] = -0.12335584105275092;
			lofilt[7] =  0.07565691101399093;
#else
			hifilt[0] =   0.064538882628938;
			hifilt[1] =  -0.040689417609558;
			hifilt[2] =  -0.41809227322221;
			hifilt[3] =   0.78848561640566;
			hifilt[4] =  -0.41809227322221;
			hifilt[5] =  -0.040689417609558;
			hifilt[6] =	  0.064538882628938;

			lofilt[0] =   0.037828455506995;
			lofilt[1] =  -0.02384946501938;
			lofilt[2] =  -0.11062440441842;
			lofilt[3] =   0.37740285561265;
			lofilt[4] =   0.85269867900940;
			lofilt[5] =   0.37740285561265;
			lofilt[6] =  -0.11062440441842;
			lofilt[7] =  -0.02384946501938;
			lofilt[8] =	  0.037828455506995;
#endif
		}
		~Wsq()
		{
			delete [] lofilt;
			delete [] hifilt;
		}
		int encode(unsigned char *idata, unsigned char *odata, int width, int height, float ratio, int chksum);
		int decode(unsigned char *odata, unsigned char *idata, int ilen, int *cols, int *rows, int chksum);
	private:
		int wsq_decode_mem(unsigned char *odata, int *ow, int *oh, int *od, int *oppi, int chksum, unsigned char *idata, const int ilen);
		int huffman_decode_data_mem(
			short *ip,               /* image pointer */
			DTT_TABLE *dtt_table,    /*transform table pointer */
			DQT_TABLE *dqt_table,    /* quantization table */
			DHT_TABLE *dht_table,    /* huffman table */
			unsigned char **cbufptr, /* points to current byte in input buffer */
			unsigned char *ebufptr);  /* points to end of input buffer */
		int decode_data_mem(
			int *onodeptr,       /* returned huffman code category        */
			int *mincode,        /* points to minimum code value for      */
					             /*    a given code length                */
			int *maxcode,        /* points to maximum code value for      */
			                     /*    a given code length                */
			int *valptr,         /* points to first code in the huffman   */
			                     /*    code table for a given code length */
			unsigned char *huffvalues,   /* defines order of huffman code          */
			                             /*    lengths in relation to code sizes   */
			unsigned char **cbufptr,     /* points to current byte in input buffer */
			unsigned char *ebufptr,      /* points to end of input buffer          */
			int *bit_count,              /* marks the bit to receive from the input byte */
			unsigned short *marker);
		int getc_nextbits_wsq(
			unsigned short *obits,       /* returned bits */
			unsigned short *marker,      /* returned marker */
			unsigned char **cbufptr,     /* points to current byte in input buffer */
			unsigned char *ebufptr,      /* points to end of input buffer */
			int *bit_count,              /* marks the bit to receive from the input byte */
			const int bits_req);         /* number of bits requested */
		int wsq_encode_mem(unsigned char *wsq_data, int *olen, const float r_bitrate,
			unsigned char *idata, const int w, const int h,
			const int d, const int ppi, char *comment_text);
		int gen_hufftable_wsq(HUFFCODE **ohufftable, unsigned char **ohuffbits,
			unsigned char **ohuffvalues, short *sip, const int *block_sizes,
			const int num_sizes);
		int compress_block(
			unsigned char *outbuf,      /* compressed output buffer            */
			int   *obytes,				/* number of compressed bytes          */
			short *sip,					/* quantized image                     */
			const int sip_siz,			/* size of quantized image to compress */
			const int MaxCoeff,			/* Maximum values for coefficients     */
			const int MaxZRun,			/* Maximum zero runs                   */
			HUFFCODE *codes);			/* huffman code table                  */
		int count_block(
			int **ocounts,				/* output count for each huffman catetory */
			const int max_huffcounts,	/* maximum number of counts */
			short *sip,					/* quantized data */
			const int sip_siz,			/* size of block being compressed */
			const int MaxCoeff,			/* maximum values for coefficients */
			const int MaxZRun);			/* maximum zero runs */
		int getc_byte(
			unsigned char *ochar_dat,	/* pointer to returned byte       */
			unsigned char **cbufptr,	/* pointer to next byte in buffer */
			unsigned char *ebufptr);    /* pointer to end of buffer       */
		int getc_bytes(
			unsigned char **ochar_dat,	/* pointer to returned bytes      */
			const int ilen,				/* number of bytes to be returned */
			unsigned char **cbufptr,	/* pointer to next byte in buffer */
			unsigned char *ebufptr);    /* pointer to end of buffer       */
		int putc_byte(
			const unsigned char idata,	/* input byte */
			unsigned char *odata,		/* output buffer of bytes */
			const int oalloc,			/* allocated size of output buffer */
			int *olen);                 /* filled length of output buffer  */
		int putc_bytes(
			unsigned char *idata,	/* input buffer of bytes           */
			const int ilen,			/* bytes to be copied              */
			unsigned char *odata,	/* output buffer of bytes          */
			const int oalloc,		/* allocated size of output buffer */
			int *olen);             /* filled length of output buffer  */
		int getc_ushort(
			unsigned short *oshrt_dat,  /* pointer to returned unsigned short */
			unsigned char  **cbufptr,   /* pointer to next byte in buffer     */
			unsigned char  *ebufptr);   /* pointer to end of buffer           */
		int putc_ushort(
			unsigned short ishort,     /* input unsigned short     */
			unsigned char *odata,      /* output byte buffer       */
			const int oalloc,          /* allocated size of buffer */
			int *olen);                /* filled length of buffer  */
		int getc_uint(
			unsigned int  *oint_dat,  /* pointer to returned unsigned int */
			unsigned char **cbufptr,  /* pointer to next byte in buffer   */
			unsigned char *ebufptr);  /* pointer to end of buffer         */
		int putc_uint(
			unsigned int iint,        /* input unsigned int       */
			unsigned char *odata,     /* output byte buffer       */
			const int oalloc,         /* allocated size of buffer */
			int *olen);               /* filled length of buffer  */
		void write_bits(
			unsigned char **outbuf,    /* output data buffer                          */
			const unsigned short code, /* info to write into buffer                   */
			const short size,          /* numbers bits of code to write into buffer   */
			int *outbit,               /* current bit location in out buffer byte     */
			unsigned char *bits,       /* byte to write to output buffer              */
			int *bytes);               /* count of number bytes written to the buffer */
		void flush_bits(
			unsigned char **outbuf,		/* output data buffer */
			int *outbit,				/* current bit location in out buffer byte */
			unsigned char *bits,		/* byte to write to output buffer */
			int *bytes);				/* count of number bytes written to the buffer */
		int check_huffcodes_wsq(HUFFCODE *hufftable, int last_size);
		int getc_huffman_table(unsigned char *otable_id, unsigned char **ohuffbits,
			unsigned char **ohuffvalues, const int max_huffcounts,
			unsigned char **cbufptr, unsigned char *ebufptr,
			const int read_table_len, int *bytes_left);
		int putc_huffman_table(
			unsigned short *factor, 
			const unsigned char table_id,   /* huffman table indicator  */
			unsigned char *huffbits,		/* huffman table parameters */
			unsigned char *huffvalues,
			unsigned char *outbuf,			/* output byte buffer       */
			const int outalloc,				/* allocated size of buffer */
			int   *outlen);					/* filled length of buffer  */
		int find_huff_sizes(int **ocodesize, int *freq, const int max_huffcounts);
		void find_least_freq(int *value1, int *value2, int *freq, const int max_huffcounts);
		int find_num_huff_sizes(unsigned char **obits, int *adjust, int *codesize, const int max_huffcounts);
		int sort_huffbits(unsigned char *bits);
		int sort_code_sizes(unsigned char **ovalues, int *codesize,	const int max_huffcounts);
		int build_huffcode_table(HUFFCODE **ohuffcode_table, HUFFCODE *in_huffcode_table, const int last_size, 
			unsigned char *values, const int max_huffcounts);
		int build_huffsizes(HUFFCODE **ohuffcode_table, int *temp_size,	unsigned char *huffbits, const int max_huffcounts);
		void build_huffcodes(HUFFCODE *huffcode_table);
		void gen_decode_table(HUFFCODE *huffcode_table,	int *maxcode, int *mincode, int *valptr, unsigned char *huffbits);
		int wsqCheckSum(unsigned char *src,  int lensrc);
		float GetBitRate(float press);
		int putc_comment(
			const unsigned short marker,
			unsigned char *comment,			/* comment */
			const int cs,					/* comment size */
			unsigned char *odata,			/* output byte buffer       */
			const int oalloc,				/* allocated size of buffer */
			int   *olen);					/* filled length of buffer  */
		int getc_marker_wsq(
			unsigned short *omarker,		/* marker read */
			const int type,					/* type of markers that could be found */
			unsigned char **cbufptr,		/* current byte in input buffer */
			unsigned char *ebufptr);		/* end of input buffer */
		int getc_table_wsq(
			unsigned short marker,			/* WSQ marker */
			DTT_TABLE *dtt_table,			/* transform table structure */
			DQT_TABLE *dqt_table,			/* quantization table structure */
			DHT_TABLE *dht_table,			/* huffman table structure */
			unsigned char **cbufptr,		/* current byte in input buffer */
			unsigned char *ebufptr);		/* end of input buffer */
		int getc_transform_table(
			DTT_TABLE *dtt_table,			/* transform table structure */
			unsigned char **cbufptr,		/* current byte in input buffer */
			unsigned char *ebufptr);		/* end of input buffer */
		int putc_transform_table(
			double *lofilt,					/* filter coefficients      */
			const int losz,
			double *hifilt,
			const int hisz,
			unsigned char *odata,			/* output byte buffer       */
			const int oalloc,				/* allocated size of buffer */
			int   *olen);					/* filled length of buffer  */
		int getc_quantization_table(
			DQT_TABLE *dqt_table,			/* quatization table structure */
			unsigned char **cbufptr,		/* current byte in input buffer */
			unsigned char *ebufptr);		/* end of input buffer */
		int putc_quantization_table(
			QUANT_VALS *quant_vals,			/* quantization parameters  */
			unsigned char *odata,           /* output byte buffer       */
			const int oalloc,				/* allocated size of buffer */
			int   *olen);					/* filled length of buffer  */
		int getc_huffman_table_wsq(
			DHT_TABLE *dht_table,			/* huffman table structure */
			unsigned char **cbufptr,		/* current byte in input buffer */
			unsigned char *ebufptr);		/* end of input buffer */
		int getc_frame_header_wsq(
			FRM_HEADER_WSQ *frm_header,		/* frame header structure */
			unsigned char **cbufptr,		/* current byte in input buffer */
			unsigned char *ebufptr);		/* end of input buffer */
		int putc_frame_header_wsq(
			const int width,				/* image width              */
			const int height,				/* image height             */
			const float m_shift,			/* image shifting parameter */
			const float r_scale,			/* image scaling parameter  */
			unsigned char *odata,			/* output byte buffer       */
			const int oalloc,				/* allocated size of buffer */
			int   *olen);					/* filled length of buffer  */
		int getc_block_header(
			unsigned char *huff_table,		/* huffman table indicator */
			unsigned char **cbufptr,		/* current byte in input buffer */
			unsigned char *ebufptr);		/* end of input buffer */
		int putc_block_header(
			const int table,				/* huffman table indicator  */
			unsigned char *odata,			/* output byte buffer       */
			const int oalloc,				/* allocated size of buffer */
			int   *olen);					/* filled length of buffer  */
		int getc_comment(
			unsigned char **ocomment,
			unsigned char **cbufptr,		/* current byte in input buffer */
			unsigned char *ebufptr);		/* end of input buffer */
		void build_wsq_trees(W_TREE w_tree[], const int w_treelen,
			Q_TREE q_tree[], const int q_treelen,
			const int width, const int height);
		void build_w_tree(
			W_TREE w_tree[],				/* wavelet tree structure */
			const int width,				/* image width            */
			const int height);				/* image height           */
		void w_tree4(
			W_TREE w_tree[],    /* wavelet tree structure                      */
			const int start1,   /* w_tree locations to start calculating       */
			const int start2,   /*    subband split locations and sizes        */
			const int lenx,     /* (temp) subband split location and sizes     */
			const int leny,
			const int x,
			const int y,
			const int stop1);		/* 0 normal operation, 1 used to avoid marking */
									/*    size and location of subbands 60-63      */
		void build_q_tree(
			W_TREE *w_tree,			/* wavelet tree structure */
			Q_TREE *q_tree);		/* quantization tree structure */
		void q_tree16(
			Q_TREE *q_tree,			/* quantization tree structure */
			const int start,		/* q_tree location of first subband        */
									/*   in the subband group being calculated */
			const int lenx,			/* (temp) subband location and sizes */
			const int leny,
			const int x,
			const int y,
			const int rw,	/* NEW */   /* spectral invert 1st row/col splits */
			const int cl);  /* NEW */
		void q_tree4(
			Q_TREE *q_tree,		/* quantization tree structure */
			const int start,	/* q_tree location of first subband         */
								/*    in the subband group being calculated */
			const int lenx,		/* (temp) subband location and sizes */
			const int leny,
			const int x,
			const int  y);        
		void conv_img_2_flt(
			float *fip,					/* output float image data  */
			float *m_shift,				/* shifting parameter       */
			float *r_scale,				/* scaling parameter        */
			unsigned char *data,        /* input unsigned char data */
			const int num_pix);			/* num pixels in image      */
		void conv_img_2_uchar(
			unsigned char *data,           /* uchar image pointer    */
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
			const int height);      /* image height            */
		int quantize(
			short **osip,           /* quantized output             */
			int *ocmp_siz,          /* size of quantized output     */
			QUANT_VALS *quant_vals, /* quantization parameters      */
			Q_TREE q_tree[],        /* quantization "tree"          */
			const int q_treelen,    /* size of q_tree               */
			float *fip,             /* floating point image pointer */
			const int width,        /* image width                  */
			const int height);      /* image height                 */
		void quant_block_sizes(int *oqsize1, int *oqsize2, int *oqsize3,
			QUANT_VALS *quant_vals,
			W_TREE w_tree[], const int w_treelen,
			Q_TREE q_tree[], const int q_treelen);
		int unquantize(
			float **ofip,				/* floating point image pointer         */
			const DQT_TABLE *dqt_table, /* quantization table structure   */
			Q_TREE q_tree[],			/* quantization table structure         */
			const int q_treelen,		/* size of q_tree                       */
			short *sip,					/* quantized image pointer              */
			const int width,			/* image width                          */
			const int height);			/* image height                         */
		int wsq_decompose(float *fdata, const int width, const int height,
			W_TREE w_tree[], const int w_treelen,
			double *hifilt, const int hisz, double *lofilt, const int losz);
		void get_lets(
			float *new1,			/* image pointers for creating subband splits */
			float *old,
			const int len1,			/* temporary length parameters */
			const int len2,
			const int pitch,		/* pitch gives next row_col to filter */
			const int  stride,		/* stride gives next pixel to filter */
			double *hi,
			const int hsz,			/* NEW */
			double *lo,				/* filter coefficients */
			const int lsz,			/* NEW */
			const int inv);			/* spectral inversion? */
		int wsq_reconstruct(float *fdata, const int width, const int height,
			W_TREE w_tree[], const int w_treelen,
			const DTT_TABLE *dtt_table);
		void  join_lets(
			float *new1,			/* image pointers for creating subband splits */
			float *old,
			const int len1,			/* temporary length parameters */
			const int len2,
			const int pitch,		/* pitch gives next row_col to filter */
			const int  stride,		/* stride gives next pixel to filter */
			double *hi,
//			float *hi,
			const int hsz,				/* NEW */
			double *lo,					/* filter coefficients */
//			float *lo,					/* filter coefficients */
			const int lsz,				/* NEW */
			const int inv);				/* spectral inversion? */
		int int_sign(const int power);  /* "sign" power */
		void init_wsq_decoder_resources(DTT_TABLE *dtt_table);
		void free_wsq_decoder_resources(DTT_TABLE *dtt_table);
	};
//}

#endif // WSQ_31_H_
