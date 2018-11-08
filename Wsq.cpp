// This is the main DLL file.
#include "stdafx.h"
#include "wsq.h"

//namespace WSQ {
/************************************************************************/
/*              This is an implementation based on the Crinimal         */
/*              Justice Information Services (CJIS) document            */
/*              "WSQ Gray-scale Fingerprint Compression                 */
/*              Specification", Dec. 1997.                              */
/***************************************************************************/
/* WSQ Decoder routine.  Takes an WSQ compressed memory buffer and decodes */
/* it, returning the reconstructed pixmap.                                 */
/***************************************************************************/
int Wsq::wsq_decode_mem(unsigned char *odata, int *ow, int *oh, int *od, int *oppi,
                   int chksum, unsigned char *idata, const int ilen)
{
   int ret, i; 
   unsigned short marker;         /* WSQ marker */
   int num_pix;             /* image size and counter */
   int width, height; //, ppi;        /* image parameters */
   float *fdata; //, *fdata1, *fdata_bse;   /* image pointers */
   short *qdata;                  /* image pointers */
   unsigned char *cbufptr;        /* points to current byte in buffer */
   unsigned char *ebufptr;        /* points to end of buffer */

   /* Set memory buffer pointers. */
   if (!idata) return -1;
	if (wsqCheckSum(idata,ilen) != chksum )  return ERROR_CHECK_SUM;

	code = 0;
	code2 = 0;

   DTT_TABLE dtt_table;
   DQT_TABLE dqt_table;
   DHT_TABLE dht_table[MAX_DHT_TABLES];

   cbufptr = idata;
   ebufptr = idata + ilen;

   /* Added by MDG on 02-24-05 */
   init_wsq_decoder_resources(&dtt_table);

   /* Init DHT Tables to 0. */
   for(i = 0; i < MAX_DHT_TABLES; i++)
      (dht_table + i)->tabdef = 0;

   /* Read the SOI marker. */
   if(ret = getc_marker_wsq(&marker, SOI_WSQ, &cbufptr, ebufptr)){
      free_wsq_decoder_resources(&dtt_table);
      return(ret);
   }

   /* Read in supporting tables up to the SOF marker. */
   if(ret = getc_marker_wsq(&marker, TBLS_N_SOF, &cbufptr, ebufptr)){
      free_wsq_decoder_resources(&dtt_table);
      return(ret);
   }
   while(marker != SOF_WSQ) {
      if(ret = getc_table_wsq(marker, &dtt_table, &dqt_table, dht_table,
                          &cbufptr, ebufptr)){
      free_wsq_decoder_resources(&dtt_table);
         return(ret);
      }
      if(ret = getc_marker_wsq(&marker, TBLS_N_SOF, &cbufptr, ebufptr)){
      free_wsq_decoder_resources(&dtt_table);
         return(ret);
      }
   }

   /* Read in the Frame Header. */
   if (ret = getc_frame_header_wsq(&frm_header_wsq, &cbufptr, ebufptr)){
      free_wsq_decoder_resources(&dtt_table);
      return(ret);
   }
   width = frm_header_wsq.width;
   height = frm_header_wsq.height;

   if (!odata) 
   {
	   *ow = width;
	   *oh = height;
	   *od = 8;
	   *oppi = 500; //ppi;
      free_wsq_decoder_resources(&dtt_table);
	   return EMPTY_SOURCE;
   }

   num_pix = width * height;

   if(debug > 0)
      fprintf(stdf, "SOI, tables, and frame header read\n\n");

	/* Build WSQ decomposition trees. */
//	W_TREE w_tree[W_TREELEN];
	//Q_TREE q_tree[Q_TREELEN];
   build_wsq_trees(w_tree, W_TREELEN, q_tree, Q_TREELEN, width, height);

   if(debug > 0)
      fprintf(stdf, "Tables for wavelet decomposition finished\n\n");

   /* Allocate working memory. */
   qdata = (short *) calloc(num_pix,  sizeof(short));
   if(qdata == (short *)NULL) {
//      fprintf(stdf,"ERROR: wsq_decode_mem : malloc : qdata1\n");
      free_wsq_decoder_resources(&dtt_table);
      return(-20);
   }
   /* Decode the Huffman encoded data blocks. */
   if(ret = huffman_decode_data_mem(qdata, &dtt_table, &dqt_table, dht_table,&cbufptr, ebufptr))
	{
      free_wsq_decoder_resources(&dtt_table);
      free(qdata);
      return(ret);
   }
#ifdef _WINDOWS
   if(debug > 0)
   {
	   int i=1;
	   char sfile[100];

	   do {
		   sprintf(sfile,"binind%02d.out",i++);
	   } while (!access(sfile,00));
	   FILE *f = fopen(sfile,"wb");
	   if (!f) return -1;
	   fwrite(qdata,num_pix * sizeof(short),1,f);
	   fclose(f);

      fprintf(stdf,
         "Quantized WSQ subband data blocks read and Huffman decoded\n\n");
   }
#endif
//t3 = clock();
 
   /* Decode the quantize wavelet subband data. */
   if(ret = unquantize(&fdata, &dqt_table, q_tree, Q_TREELEN,qdata, width, height))
	{
      free_wsq_decoder_resources(&dtt_table);
      free(qdata);
      return(ret);
   }

   if(debug > 0)
      fprintf(stdf, "WSQ subband data blocks unquantized\n\n");

   /* Done with quantized wavelet subband data. */
   free(qdata);
//t4 = clock();

   if (ret = wsq_reconstruct(fdata, width, height, w_tree, W_TREELEN,&dtt_table))
	{
      free_wsq_decoder_resources(&dtt_table);
      free(fdata);
      return(ret);
   }

   if(debug > 0)
      fprintf(stdf, "WSQ reconstruction of image finished\n\n");
//t5 = clock();
/*
   cdata = (unsigned char *)malloc(num_pix * sizeof(unsigned char));
   if(cdata == (unsigned char *)NULL) {
      free(fdata);
      fprintf(stdf,"ERROR: wsq_decode_mem : malloc : cdata\n");
      return(-21);
   }
  */
   /* Convert floating point pixels to unsigned char pixels. */
   conv_img_2_uchar(odata, fdata, width, height,frm_header_wsq.m_shift, frm_header_wsq.r_scale);

   /* Done with floating point pixels. */
   free(fdata);

   free_wsq_decoder_resources(&dtt_table);

   if(debug > 0)
      fprintf(stdf, "Doubleing point pixels converted to unsigned char\n\n");

/*
t6 = clock();
t6-=t5;
t5-=t4;
t4-=t3;
t3-=t2; 
t2-=t1;
t1-=t0;
*/
   /* Assign reconstructed pixmap and attributes to output pointers. */
//   *odata = cdata;
   *ow = width;
   *oh = height;
   *od = 8;
   *oppi = 500; //ppi;

   /* Return normally. */
	return  frm_header_wsq.software;
}

/***************************************************************************/
/* Routine to decode an entire "block" of encoded data from memory buffer. */
/***************************************************************************/
int Wsq::huffman_decode_data_mem(
   short *ip,               /* image pointer */
   DTT_TABLE *dtt_table,    /*transform table pointer */
   DQT_TABLE *dqt_table,    /* quantization table */
   DHT_TABLE *dht_table,    /* huffman table */
   unsigned char **cbufptr, /* points to current byte in input buffer */
   unsigned char *ebufptr)  /* points to end of input buffer */
{
   int ret;
   int blk = 0;           /* block number */
   unsigned short marker; /* WSQ markers */
   int bit_count;         /* bit count for getc_nextbits_wsq routine */
   int n;                 /* zero run count */
   int nodeptr;           /* pointers for decoding */
   int last_size;         /* last huffvalue */
   unsigned char hufftable_id;    /* huffman table number */
   HUFFCODE *hufftable;   /* huffman code structure */
   int maxcode[MAX_HUFFBITS+1]; /* used in decoding data */
   int mincode[MAX_HUFFBITS+1]; /* used in decoding data */
   int valptr[MAX_HUFFBITS+1];     /* used in decoding data */
   unsigned short tbits;
   int ipc, ipc_mx, ipc_q;   /* image byte count adjustment parameters */


   if(ret = getc_marker_wsq(&marker, TBLS_N_SOB, cbufptr, ebufptr))
      return(ret);

   bit_count = 0;
   ipc = 0;
   ipc_q = 0;
   ipc_mx = frm_header_wsq.width * frm_header_wsq.height;

   while(marker != EOI_WSQ) {

      if(marker != 0) {
         blk++;
         while(marker != SOB_WSQ) {
            if(ret = getc_table_wsq(marker, dtt_table, dqt_table,
                                dht_table, cbufptr, ebufptr))
               return(ret);
            if(ret = getc_marker_wsq(&marker, TBLS_N_SOB, cbufptr, ebufptr))
               return(ret);
         }
         if(dqt_table->dqt_def && !ipc_q) {
            for(n = 0; n < 64; n++)
               if(dqt_table->q_bin[n] == 0.0)
                  ipc_mx -= q_tree[n].lenx*q_tree[n].leny;

            ipc_q = 1;
         }
         if(ret = getc_block_header(&hufftable_id, cbufptr, ebufptr))
            return(ret);

         if((dht_table+hufftable_id)->tabdef != 1) {
//            fprintf(stdf, "ERROR : huffman_decode_data_mem : ");
  //          fprintf(stdf, "huffman table {%d} undefined.\n", hufftable_id);
            return(-51);
         }

         /* the next two routines reconstruct the huffman tables */
         if(ret = build_huffsizes(&hufftable, &last_size,
                                  (dht_table+hufftable_id)->huffbits,
                                  MAX_HUFFCOUNTS_WSQ))
            return(ret);

         build_huffcodes(hufftable);
         if (ret = check_huffcodes_wsq(hufftable, last_size))
			{
				if (debug)
            fprintf(stdf, "         hufftable_id = %d\n", hufftable_id);
			}

         /* this routine builds a set of three tables used in decoding */
         /* the compressed data*/
         gen_decode_table(hufftable, maxcode, mincode, valptr, (dht_table+hufftable_id)->huffbits);
         free(hufftable);
         bit_count = 0;
         marker = 0;
      }

      /* get next huffman category code from compressed input data stream */
      if(ret = decode_data_mem(&nodeptr, mincode, maxcode, valptr, (dht_table+hufftable_id)->huffvalues, cbufptr, ebufptr, &bit_count, &marker))
         return(ret);

      if(nodeptr == -1) {
         while(marker == COM_WSQ && blk == 3) {
            if((ret = getc_table_wsq(marker, dtt_table, dqt_table,
                                dht_table, cbufptr, ebufptr)))
               return(ret);
            if((ret = getc_marker_wsq(&marker, ANY_WSQ, cbufptr, ebufptr)))
               return(ret);
         }
         continue;
      }

      if(ipc > ipc_mx) {
		  if (debug)
		  {
         fprintf(stdf, "ERROR : huffman_decode_data_mem [1]: ");
         fprintf(stdf, "Decoded data extends past image buffer. ");
         fprintf(stdf, "Encoded data appears corrupt or non-standard.\n");
		  }
         return(-510);
      }

      if(nodeptr > 0 && nodeptr <= 100) {
         ipc += nodeptr;
         if(ipc > ipc_mx) {
			 if (debug)
			 {
            fprintf(stdf, "ERROR : huffman_decode_data_mem [2]: ");
            fprintf(stdf, "Decoded data extends past image buffer. ");
            fprintf(stdf, "Encoded data appears corrupt or non-standard.\n");
			 }
            return(-511);
         }
         for(n = 0; n < nodeptr; n++)
            *ip++ = 0; /* z run */
      }
      else if(nodeptr > 106 && nodeptr < 0xff){
         *ip++ = nodeptr - 180;
         ipc++;
      }
      else if(nodeptr == 101){
         if(ret = getc_nextbits_wsq(&tbits, &marker, cbufptr, ebufptr, &bit_count, 8))
            return(ret);
         *ip++ = tbits;
         ipc++;
      }
      else if(nodeptr == 102){
         if(ret = getc_nextbits_wsq(&tbits, &marker, cbufptr, ebufptr, &bit_count, 8))
            return(ret);
         *ip++ = -tbits;
         ipc++;
      }
      else if(nodeptr == 103){
         if(ret = getc_nextbits_wsq(&tbits, &marker, cbufptr, ebufptr, &bit_count, 16))
            return(ret);
         *ip++ = tbits;
         ipc++;
      }
      else if(nodeptr == 104){
         if(ret = getc_nextbits_wsq(&tbits, &marker, cbufptr, ebufptr, &bit_count, 16))
            return(ret);
         *ip++ = -tbits;
         ipc++;
      }
      else if(nodeptr == 105) {
         if(ret = getc_nextbits_wsq(&tbits, &marker, cbufptr, ebufptr, &bit_count, 8))
            return(ret);
         ipc+=tbits;
         if(ipc > ipc_mx) {
			 if (debug)
			 {
            fprintf(stdf, "ERROR : huffman_decode_data_mem [3]: ");
            fprintf(stdf, "Decoded data extends past image buffer. ");
            fprintf(stdf, "Encoded data appears corrupt or non-standard.\n");
			 }
            return(-512);
         }
         n = tbits;
         while(n--)
            *ip++ = 0;
      }
      else if(nodeptr == 106) {
         if(ret = getc_nextbits_wsq(&tbits, &marker, cbufptr, ebufptr, &bit_count, 16))
            return(ret);
         ipc+=tbits;
         if(ipc > ipc_mx) {
			 if (debug)
			 {
            fprintf(stdf, "ERROR : huffman_decode_data_mem [4]: ");
            fprintf(stdf, "Decoded data extends past image buffer. ");
            fprintf(stdf, "Encoded data appears corrupt or non-standard.\n");
			 }
            return(0); //-513);
         }
         n = tbits;
         while(n--)
            *ip++ = 0;
      }
      else {
//         fprintf(stdf, 
  //              "ERROR: huffman_decode_data_mem : Invalid code %d (%x).\n",
    //            nodeptr, nodeptr);
         return(-52);
      }

   }

   return(0);
}


/**********************************************************/
/* Routine to decode the encoded data from memory buffer. */
/**********************************************************/
int Wsq::decode_data_mem(
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
   int *bit_count,      /* marks the bit to receive from the input byte */
   unsigned short *marker)
{
   int ret;
   int inx, inx2;       /*increment variables*/
   unsigned short code, tbits;  /* becomes a huffman code word
                                   (one bit at a time)*/

   if(ret = getc_nextbits_wsq(&code, marker, cbufptr, ebufptr, bit_count, 1))
      return(ret);

   if(*marker != 0){
      *onodeptr = -1;
      return(0);
   }

   for(inx = 1; (int)code > maxcode[inx] && inx<=MAX_HUFFBITS; inx++) {
      if(ret = getc_nextbits_wsq(&tbits, marker, cbufptr, ebufptr, bit_count, 1))
         return(ret);

      code = (code << 1) + tbits;
      if(*marker != 0){
         *onodeptr = -1;
         return(0);
      }
   }
//   if (code==65535)
   if (inx==17)
   {
	   //FILE *f = fopen("minmaxvalptr.txt","w");
	   //for (int i=0; i<=MAX_HUFFBITS; i++)
		  // fprintf(f,"min=%d max=%d valptr=%d\n",mincode[i],maxcode[i],valptr[i]);
	   //fclose(f);
	   return(0);
   }
   inx2 = valptr[inx]; 
   inx2 = inx2 + code - mincode[inx];

   *onodeptr = huffvalues[inx2];
   return(0);
}

/****************************************************************/
/* Routine to get nextbit(s) of data stream from memory buffer. */
/****************************************************************/
int Wsq::getc_nextbits_wsq(
   unsigned short *obits,       /* returned bits */
   unsigned short *marker,      /* returned marker */
   unsigned char **cbufptr,     /* points to current byte in input buffer */
   unsigned char *ebufptr,      /* points to end of input buffer */
   int *bit_count,      /* marks the bit to receive from the input byte */
   const int bits_req)  /* number of bits requested */
{
   int ret;
   //static unsigned char code;   /*next byte of data*/
   //static unsigned char code2;  /*stuffed byte of data*/
   unsigned short bits, tbits;  /*bits of current data byte requested*/
   int bits_needed;     /*additional bits required to finish request*/

                              /*used to "mask out" n number of
                                bits from data stream*/
//   static unsigned char bit_mask[9] = {0x00,0x01,0x03,0x07,0x0f,0x1f,0x3f,0x7f,0xff};
   if(*bit_count == 0) {
      if(ret = getc_byte(&code, cbufptr, ebufptr)){
         return(ret);
      }
      *bit_count = 8;
      if(code == 0xFF) {
         if(ret = getc_byte(&code2, cbufptr, ebufptr)){
            return(ret);
         }
         if(code2 != 0x00 && bits_req == 1) {
            *marker = (code << 8) | code2;
            *obits = 1;
            return(0);
         }
         if(code2 != 0x00) {
//            fprintf(stdf, "ERROR: getc_nextbits_wsq : No stuffed zeros\n");
            return(-41);
         }
      }
   }
   if(bits_req <= *bit_count) {
      bits = (code >>(*bit_count - bits_req)) & (bit_mask[bits_req]);
      *bit_count -= bits_req;
      code &= bit_mask[*bit_count];
   }
   else {
      bits_needed = bits_req - *bit_count;
      bits = code << bits_needed;
      *bit_count = 0;
      if(ret = getc_nextbits_wsq(&tbits, (unsigned short *)NULL, cbufptr, ebufptr, bit_count, bits_needed))
         return(ret);
      bits |= tbits;
   }

   *obits = bits;
   return(0);
}
/************************************************************************/
/*              This is an implementation based on the Crinimal         */
/*              Justice Information Services (CJIS) document            */
/*              "WSQ Gray-scale Fingerprint Compression                 */
/*              Specification", Dec. 1997.                              */
/************************************************************************/
/* WSQ encodes/compresses an image pixmap.                              */
/************************************************************************/
int Wsq::wsq_encode_mem(unsigned char *wsq_data, int *olen, const float r_bitrate,
                   unsigned char *idata, const int w, const int h,
                   const int d, const int ppi, char *comment_text)
{
   int ret, num_pix;
   float *fdata;                 /* floating point pixel image  */
   float m_shift, r_scale;       /* shift/scale parameters      */
   short *qdata;                 /* quantized image pointer     */
   int qsize, qsize1, qsize2, qsize3;  /* quantized block sizes */
   unsigned char *huffbits, *huffvalues; /* huffman code parameters     */
   HUFFCODE *hufftable;          /* huffcode table              */
   HUFFCODE *hufftable23;          /* huffcode table              */
   unsigned char *huff_buf;      /* huffman encoded buffer      */ 
   int hsize, hsize1, hsize2, hsize3; /* Huffman coded blocks sizes */
//   unsigned char *wsq_data;      /* compressed data buffer      */
   int wsq_alloc, wsq_len;       /* number of bytes in buffer   */
   int block_sizes[3];

   //DTT_TABLE dtt_table;
   //DQT_TABLE dqt_table;
   //DHT_TABLE dht_table[MAX_DHT_TABLES];

   /* Compute the total number of pixels in image. */
   num_pix = w * h;

   /* Allocate floating point pixmap. */
   if((fdata = (float *) malloc(num_pix*sizeof(float))) == NULL) {
//      fprintf(stdf,"ERROR : wsq_encode_1 : malloc : fdata\n");
      return(-10);
   }

   /* Convert image pixels to floating point. */
   conv_img_2_flt(fdata, &m_shift, &r_scale, idata, num_pix);

   if(debug > 0)
      fprintf(stdf, "Input image pixels converted to floating point\n\n");

   /* Build WSQ decomposition trees */
   build_wsq_trees(w_tree, W_TREELEN, q_tree, Q_TREELEN, w, h);

   if(debug > 0)
      fprintf(stdf, "Tables for wavelet decomposition finished\n\n");

   /* WSQ decompose the image */
   if(ret = wsq_decompose(fdata, w, h, w_tree, W_TREELEN,
                            hifilt, MAX_HIFILT, lofilt, MAX_LOFILT)){
      free(fdata);
      return(ret);
   }

   if(debug > 0)
      fprintf(stdf, "WSQ decomposition of image finished\n\n");

   /* Set compression ratio and 'q' to zero. */
   QUANT_VALS quant_vals;
   quant_vals.cr = 0;
   quant_vals.q = 0.0;
   /* Assign specified r-bitrate into quantization structure. */
   quant_vals.r = r_bitrate;
   /* Compute subband variances. */
   variance(&quant_vals, q_tree, Q_TREELEN, fdata, w, h);

   if(debug > 0)
      fprintf(stdf, "Subband variances computed\n\n");

   /* Quantize the floating point pixmap. */
   if(ret = quantize(&qdata, &qsize, &quant_vals, q_tree, Q_TREELEN,
                      fdata, w, h)){
      free(fdata);
      return(ret);
   }

   /* Done with floating point wsq subband data. */
   free(fdata);
#ifdef _WINDOWS
   if(debug > 0)
   {
	   int i=1;
	   char sfile[100];
	   do {
		   sprintf(sfile,"binind%02d.in",i++);
	   } while (!access(sfile,00));
	   FILE *f = fopen(sfile,"wb");
	   if (!f) return -1;
	   fwrite(qdata, w * h * sizeof(short),1,f);
	   fclose(f);

      fprintf(stdf, "WSQ subband decomposition data quantized\n\n");
   }
#endif
   /* Compute quantized WSQ subband block sizes */
   quant_block_sizes(&qsize1, &qsize2, &qsize3, &quant_vals,
                           w_tree, W_TREELEN, q_tree, Q_TREELEN);

   if(qsize != qsize1+qsize2+qsize3){
//      fprintf(stdf,
  //            "ERROR : wsq_encode_1 : problem w/quantization block sizes\n");
      return(-11);
   }

   /* Allocate a WSQ-encoded output buffer.  Allocate this buffer */
   /* to be the size of the original pixmap.  If the encoded data */
   /* exceeds this buffer size, then throw an error because we do */
   /* not want our compressed data to be larger than the original */
   /* image data.                                                 */
   /*
   wsq_data = (unsigned char *)malloc(num_pix);
   if(wsq_data == (unsigned char *)NULL){
      free(qdata);
      fprintf(stdf, "ERROR : wsq_encode_1 : malloc : wsq_data\n");
      return(-12);
   }
   */
   wsq_alloc = num_pix;
   wsq_len = 0;

   /* Add a Start Of Image (SOI) marker to the WSQ buffer. */
   if(ret = putc_ushort(SOI_WSQ, wsq_data, wsq_alloc, &wsq_len)){
      free(qdata);
//      free(wsq_data);
      return(ret);
   }
/*
   if(ret = putc_nistcom_wsq(comment_text, w, h, d, ppi, 1,
                             r_bitrate, wsq_data, wsq_alloc, &wsq_len)){
      free(qdata);
//      free(wsq_data);
      return(ret);
   }
*/
   if(ret = putc_comment(COM_WSQ, (unsigned char *)comment_text, (const int)strlen(comment_text),
                            wsq_data, wsq_alloc, &wsq_len))
	{
      free(qdata);
      return(ret);
	}
	int check_sum_offset = wsq_len-1;

   /* Store the Wavelet filter taps to the WSQ buffer. */
   if(ret = putc_transform_table(lofilt, MAX_LOFILT, hifilt, MAX_HIFILT, wsq_data, wsq_alloc, &wsq_len))
   {
      free(qdata);
//      free(wsq_data);
      return(ret);
   }
 
   /* Store the quantization parameters to the WSQ buffer. */
   if(ret = putc_quantization_table(&quant_vals, wsq_data, wsq_alloc, &wsq_len))
   {
      free(qdata);
//      free(wsq_data);
      return(ret);
   }

   /* Allocate a temporary buffer for holding compressed block data.    */
   /* This buffer is allocated to the size of the original input image, */
   /* and it is "assumed" that the compressed blocks will not exceed    */
   /* this buffer size.                                                 */
   huff_buf = (unsigned char *)malloc(num_pix);
   if(huff_buf == (unsigned char *)NULL) {
      free(qdata);
//      free(wsq_data);
//      fprintf(stdf, "ERROR : wsq_encode_1 : malloc : huff_buf\n");
      return(-13);
   }

   /******************/
   /* ENCODE Block 1 */
   /******************/
   /* Compute Huffman table for Block 1. */
   if(ret = gen_hufftable_wsq(&hufftable, &huffbits, &huffvalues, qdata, &qsize1, 1))
   {
      free(qdata);
//      free(wsq_data);
      free(huff_buf);
      return(ret);
   }

   /* Store a frame header to the WSQ buffer. */
   if(ret = putc_frame_header_wsq(w, h, m_shift, r_scale, wsq_data, wsq_alloc, &wsq_len))
   {
      free(qdata);
//      free(wsq_data);
      return(ret);
   }

   if(debug > 0)
      fprintf(stdf, "SOI, tables, and frame header written\n\n");

   /* DHT */
   if(ret = putc_ushort(DHT_WSQ, wsq_data, wsq_alloc, &wsq_len))
      return(ret);
	int DHT_offset = wsq_len;

	/* Store Huffman table for Block 1 to WSQ buffer. */
   unsigned short factor = 1;
	if(ret = putc_huffman_table(&factor, 0, huffbits, huffvalues,
                               wsq_data, wsq_alloc, &wsq_len)){
      free(qdata);
//      free(wsq_data);
      free(huff_buf);
      free(huffbits);
      free(huffvalues);
      free(hufftable);
      return(ret);
   }
	int DHT_length = factor;
   free(huffbits);
   free(huffvalues);
   /* Compute  Huffman table for Blocks 2 & 3. */
   block_sizes[0] = qsize2;
   block_sizes[1] = qsize3;
   if(ret = gen_hufftable_wsq(&hufftable23, &huffbits, &huffvalues,
                          qdata+qsize1, block_sizes, 2)){
      free(qdata);
//      free(wsq_data);
      free(huff_buf);
      return(ret);
   }
	factor = 0;
   // Store Huffman table for Blocks 2 & 3 to WSQ buffer. 
   if(ret = putc_huffman_table(&factor, 1, huffbits, huffvalues,
                               wsq_data, wsq_alloc, &wsq_len)){
      free(qdata);
//      free(wsq_data);
      free(huff_buf);
      free(huffbits);
      free(huffvalues);
      free(hufftable23);
      free(hufftable);
      return(ret);
   }
	DHT_length += factor-2; // 3+MAX_HUFFBITS twice
	putc_ushort(DHT_length, wsq_data, wsq_alloc, &DHT_offset);

   free(huffbits);
   free(huffvalues);

   if(debug > 0)
      fprintf(stdf, "Huffman code Table 1 generated and written\n\n");

   /* Compress Block 1 data. */

   if(ret = compress_block(huff_buf, &hsize1, qdata, qsize1,
                           MAX_HUFFCOEFF, MAX_HUFFZRUN, hufftable)){
      free(qdata);
//      free(wsq_data);
      free(huff_buf);
      free(hufftable);
      return(ret);
   }
   /* Done with current Huffman table. */
   free(hufftable);

   /* Accumulate number of bytes compressed. */
   hsize = hsize1;

   /* Store Block 1's header to WSQ buffer. */
   if(ret = putc_block_header(0, wsq_data, wsq_alloc, &wsq_len)){
      free(qdata);
//      free(wsq_data);
      free(huff_buf);
      return(ret);
   }

   /* Store Block 1's compressed data to WSQ buffer. */
   if(ret = putc_bytes(huff_buf, hsize1, wsq_data, wsq_alloc, &wsq_len)){
      free(qdata);
  //    free(wsq_data);
      free(huff_buf);
      return(ret);
   }

   if(debug > 0)
      fprintf(stdf, "Block 1 compressed and written\n\n");

   // Compress Block 2 data. 
   if(ret = compress_block(huff_buf, &hsize2, qdata+qsize1, qsize2,
                           MAX_HUFFCOEFF, MAX_HUFFZRUN, hufftable23)){
      free(qdata);
//      free(wsq_data);
      free(huff_buf);
      free(hufftable23);
      return(ret);
   }

   /* Accumulate number of bytes compressed. */
   hsize += hsize2;

   /* Store Block 2's header to WSQ buffer. */
   if(ret = putc_block_header(1, wsq_data, wsq_alloc, &wsq_len)){
      free(qdata);
//      free(wsq_data);
      free(huff_buf);
      free(hufftable23);
      return(ret);
   }

   /* Store Block 2's compressed data to WSQ buffer. */
   if(ret = putc_bytes(huff_buf, hsize2, wsq_data, wsq_alloc, &wsq_len)){
      free(qdata);
//      free(wsq_data);
      free(huff_buf);
      free(hufftable23);
      return(ret);
   }

   if(debug > 0)
      fprintf(stdf, "Block 2 compressed and written\n\n");

   /******************/
   /* ENCODE Block 3 */
   /******************/
   /* Compress Block 3 data. */
   if(ret = compress_block(huff_buf, &hsize3, qdata+qsize1+qsize2, qsize3,
                           MAX_HUFFCOEFF, MAX_HUFFZRUN, hufftable23)){
      free(qdata);
//      free(wsq_data);
      free(huff_buf);
      free(hufftable23);
      return(ret);
   }
   /* Done with current Huffman table. */
   free(hufftable23);

   /* Done with quantized image buffer. */
   free(qdata);

   /* Accumulate number of bytes compressed. */
   hsize += hsize3;

   /* Store Block 3's header to WSQ buffer. */
   if(ret = putc_block_header(1, wsq_data, wsq_alloc, &wsq_len)){
//      free(wsq_data);
      free(huff_buf);
      return(ret);
   }

   /* Store Block 3's compressed data to WSQ buffer. */
   if(ret = putc_bytes(huff_buf, hsize3, wsq_data, wsq_alloc, &wsq_len)){
//      free(wsq_data);
      free(huff_buf);
      return(ret);
   }

   if(debug > 0)
      fprintf(stdf, "Block 3 compressed and written\n\n");

   /* Done with huffman compressing blocks, so done with buffer. */
   free(huff_buf);

   /* Add a End Of Image (EOI) marker to the WSQ buffer. */
   if(ret = putc_ushort(EOI_WSQ, wsq_data, wsq_alloc, &wsq_len)){
//      free(wsq_data);
      return(ret);
   }

   if(debug >= 1) {
      fprintf(stdf,
              "hsize1 = %d :: hsize2 = %d :: hsize3 = %d\n", hsize1, hsize2, hsize3);
      fprintf(stdf,"@ r = %.3f :: complen = %d :: ratio = %.1f\n",
              r_bitrate, hsize, (float)(num_pix)/(float)hsize);
   }
	
	unsigned char cs = 0;
	for (int i=0; i<wsq_len; i++) cs += wsq_data[i];
	wsq_data[check_sum_offset] = -cs;

//	cs = 0;
//	for (int i=0; i<wsq_len; i++) cs += wsq_data[i];

//   *odata = wsq_data;
   *olen = wsq_len;

   /* Return normally. */
   return(0);
}

/*************************************************************/
/* Generate a Huffman code table for a quantized data block. */
/*************************************************************/
int Wsq::gen_hufftable_wsq(HUFFCODE **ohufftable, unsigned char **ohuffbits,
               unsigned char **ohuffvalues, short *sip, const int *block_sizes,
               const int num_sizes)
{
   int i, j;
   int ret;
   int adjust;          /* tells if codesize is greater than MAX_HUFFBITS */
   int *codesize;       /* code sizes to use */
   int last_size;       /* last huffvalue */
   unsigned char *huffbits;     /* huffbits values */
   unsigned char *huffvalues;   /* huffvalues */
   int *huffcounts;     /* counts for each huffman category */
   int *huffcounts2;    /* counts for each huffman category */
   HUFFCODE *hufftable1, *hufftable2;  /* hufftables */

   if(ret = count_block(&huffcounts, MAX_HUFFCOUNTS_WSQ,
                        sip, block_sizes[0], MAX_HUFFCOEFF, MAX_HUFFZRUN))
      return(ret);
///
   for(i = 1; i < num_sizes; i++) {
      if(ret = count_block(&huffcounts2, MAX_HUFFCOUNTS_WSQ,
                           sip+block_sizes[i-1], block_sizes[i],
                           MAX_HUFFCOEFF, MAX_HUFFZRUN))
         return(ret);
/*
	int offset = block_sizes[0];
   for(i = 1; i < num_sizes; i++) {
      if(ret = count_block(&huffcounts2, MAX_HUFFCOUNTS_WSQ,
                           sip+offset, block_sizes[i],
                           MAX_HUFFCOEFF, MAX_HUFFZRUN))
         return(ret);
		offset += block_sizes[i];
*/
      for(j = 0; j < MAX_HUFFCOUNTS_WSQ; j++)
         huffcounts[j] += huffcounts2[j];

      free(huffcounts2);
   }

   if(ret = find_huff_sizes(&codesize, huffcounts, MAX_HUFFCOUNTS_WSQ))
   {
      free(huffcounts);
      return(ret);
   }
   free(huffcounts);

   if(ret = find_num_huff_sizes(&huffbits, &adjust, codesize, MAX_HUFFCOUNTS_WSQ))
   {
      free(codesize);
      return(ret);
   }

   if(adjust){
      if(ret = sort_huffbits(huffbits)){
         free(codesize);
         free(huffbits);
         return(ret);
      }
   }

   if(ret = sort_code_sizes(&huffvalues, codesize, MAX_HUFFCOUNTS_WSQ))
   {
      free(codesize);
      free(huffbits);
      return(ret);
   }
   free(codesize);

   if(ret = build_huffsizes(&hufftable1, &last_size,
                              huffbits, MAX_HUFFCOUNTS_WSQ)){
      free(huffbits);
      free(huffvalues);
      return(ret);
   }

   build_huffcodes(hufftable1);
   if(ret = check_huffcodes_wsq(hufftable1, last_size)){
//      fprintf(stdf, "ERROR: This huffcode warning is an error ");
  //    fprintf(stdf, "for the encoder.\n");
      free(huffbits);
      free(huffvalues);
      free(hufftable1);
      return(ret);
   }

   if(ret = build_huffcode_table(&hufftable2, hufftable1, last_size,
                                 huffvalues, MAX_HUFFCOUNTS_WSQ)){
      free(huffbits);
      free(huffvalues);
      free(hufftable1);
      return(ret);
   }

   free(hufftable1);

   *ohuffbits = huffbits;
   *ohuffvalues = huffvalues;
   *ohufftable = hufftable2;

   return(0);
}

/*****************************************************************/
/* Routine "codes" the quantized image using the huffman tables. */
/*****************************************************************/
int Wsq::compress_block(
   unsigned char *outbuf,       /* compressed output buffer            */
   int   *obytes,       /* number of compressed bytes          */
   short *sip,          /* quantized image                     */
   const int sip_siz,   /* size of quantized image to compress */
   const int MaxCoeff,  /* Maximum values for coefficients     */
   const int MaxZRun,   /* Maximum zero runs                   */
   HUFFCODE *codes)     /* huffman code table                  */
{
   unsigned char *optr;
   int LoMaxCoeff;        /* lower (negative) MaxCoeff limit */
   short pix;             /* temp pixel pointer */
   unsigned int rcnt = 0, state;  /* zero run count and if current pixel
                             is in a zero run or just a coefficient */
   int cnt;               /* pixel counter */
   int outbit, bytes;     /* parameters used by write_bits to */
   unsigned char bits;            /* output the "coded" image to the  */
                          /* output buffer                    */

   LoMaxCoeff = 1 - MaxCoeff;
   optr = outbuf;
   outbit = 7;
   bytes = 0;
   bits = 0;
   state = COEFF_CODE;
   for (cnt = 0; cnt < sip_siz; cnt++) {
      pix = *(sip + cnt);

      switch (state) {

         case COEFF_CODE:
            if (pix == 0) {
               state = RUN_CODE;
               rcnt = 1;
               break;
            }
            if (pix > MaxCoeff) { 
               if (pix > 255) {
                  /* 16bit pos esc */
                  write_bits( &optr, (unsigned short) codes[103].code,
                              codes[103].size, &outbit, &bits, &bytes );
                  write_bits( &optr, (unsigned short) pix, 16,
                              &outbit, &bits, &bytes);
               }
               else {
                  /* 8bit pos esc */
                  write_bits( &optr, (unsigned short) codes[101].code,
                              codes[101].size, &outbit, &bits, &bytes );
                  write_bits( &optr, (unsigned short) pix, 8,
                              &outbit, &bits, &bytes);
               }
            }
            else if (pix < LoMaxCoeff) {
               if (pix < -255) {
                  /* 16bit neg esc */
                  write_bits( &optr, (unsigned short) codes[104].code,
                              codes[104].size, &outbit, &bits, &bytes );
                  write_bits( &optr, (unsigned short) -pix, 16,
                              &outbit, &bits, &bytes);
               }
               else {
                  /* 8bit neg esc */
                  write_bits( &optr, (unsigned short) codes[102].code,
                              codes[102].size, &outbit, &bits, &bytes );
                  write_bits( &optr, (unsigned short) -pix, 8,
                              &outbit, &bits, &bytes);
               }
            }
            else {
               /* within table */
               write_bits( &optr, (unsigned short) codes[pix+180].code,
                           codes[pix+180].size, &outbit, &bits, &bytes);
            }
            break;

         case RUN_CODE:
            if (pix == 0  &&  rcnt < 0xFFFF) {
               ++rcnt;
               break;
            }
            if (rcnt <= MaxZRun) {
               /* log zero run length */
               write_bits( &optr, (unsigned short) codes[rcnt].code,
                           codes[rcnt].size, &outbit, &bits, &bytes );
            }
            else if (rcnt <= 0xFF) {
               /* 8bit zrun esc */
               write_bits( &optr, (unsigned short) codes[105].code,
                           codes[105].size, &outbit, &bits, &bytes );
               write_bits( &optr, (unsigned short) rcnt, 8,
                           &outbit, &bits, &bytes);
            }
            else if (rcnt <= 0xFFFF) {
               /* 16bit zrun esc */
               write_bits( &optr, (unsigned short) codes[106].code,
                           codes[106].size, &outbit, &bits, &bytes );
               write_bits( &optr, (unsigned short) rcnt, 16,
                           &outbit, &bits, &bytes);
            }
            else {
//               fprintf(stdf,
  //                    "ERROR : compress_block : zrun too large.\n");
               return(-47);
            }

            if(pix != 0) {
               if (pix > MaxCoeff) {
                  /** log current pix **/
                  if (pix > 255) {
                     /* 16bit pos esc */
                     write_bits( &optr, (unsigned short) codes[103].code,
                                 codes[103].size, &outbit, &bits, &bytes );
                     write_bits( &optr, (unsigned short) pix, 16,
                                 &outbit, &bits, &bytes);
                  }
                  else {
                     /* 8bit pos esc */
                     write_bits( &optr, (unsigned short) codes[101].code,
                                 codes[101].size, &outbit, &bits, &bytes );
                     write_bits( &optr, (unsigned short) pix, 8,
                                 &outbit, &bits, &bytes);
                  }
               }
               else if (pix < LoMaxCoeff) {
                  if (pix < -255) {
                     /* 16bit neg esc */
                     write_bits( &optr, (unsigned short) codes[104].code,
                                 codes[104].size, &outbit, &bits, &bytes );
                     write_bits( &optr, (unsigned short) -pix, 16,
                                 &outbit, &bits, &bytes);
                  }
                  else {
                     /* 8bit neg esc */
                     write_bits( &optr, (unsigned short) codes[102].code,
                                 codes[102].size, &outbit, &bits, &bytes );
                     write_bits( &optr, (unsigned short) -pix, 8,
                                 &outbit, &bits, &bytes);
                  }
               }
               else {
                  /* within table */
                  write_bits( &optr, (unsigned short) codes[pix+180].code,
                              codes[pix+180].size, &outbit, &bits, &bytes);
               }
               state = COEFF_CODE;
            }
            else {
               rcnt = 1;
               state = RUN_CODE;
            }
            break;
      }
   }
   if (state == RUN_CODE) {
      if (rcnt <= MaxZRun) {
         write_bits( &optr, (unsigned short) codes[rcnt].code,
                     codes[rcnt].size, &outbit, &bits, &bytes );
      }
      else if (rcnt <= 0xFF) {
         write_bits( &optr, (unsigned short) codes[105].code,
                     codes[105].size, &outbit, &bits, &bytes );
         write_bits( &optr, (unsigned short) rcnt, 8,
                     &outbit, &bits, &bytes);
      }
      else if (rcnt <= 0xFFFF) {
         write_bits( &optr, (unsigned short) codes[106].code,
                     codes[106].size, &outbit, &bits, &bytes );
         write_bits( &optr, (unsigned short) rcnt, 16,
                     &outbit, &bits, &bytes);
      }
      else {
//         fprintf(stdf, "ERROR : compress_block : zrun2 too large.\n");
         return(-48);
      }
   }

   flush_bits( &optr, &outbit, &bits, &bytes);

   *obytes = bytes;
   return(0);
}

/*****************************************************************/
/* This routine counts the number of occurences of each category */
/* in the huffman coding tables.                                 */
/*****************************************************************/
int Wsq::count_block(
   int **ocounts,     /* output count for each huffman catetory */
   const int max_huffcounts, /* maximum number of counts */
   short *sip,          /* quantized data */
   const int sip_siz,   /* size of block being compressed */
   const int MaxCoeff,  /* maximum values for coefficients */
   const int MaxZRun)   /* maximum zero runs */
{
   int *counts;         /* count for each huffman category */
   int LoMaxCoeff;        /* lower (negative) MaxCoeff limit */
   short pix;             /* temp pixel pointer */
   unsigned int rcnt = 0, state;  /* zero run count and if current pixel
                             is in a zero run or just a coefficient */
   int cnt;               /* pixel counter */

   /* Ininitalize vector of counts to 0. */
   counts = (int *)calloc(max_huffcounts+1, sizeof(int));
   if(counts == (int *)NULL){
//      fprintf(stdf,
  //    "ERROR : count_block : calloc : counts\n");
      return(-48);
   }
   /* Set last count to 1. */
   counts[max_huffcounts] = 1;

   LoMaxCoeff = 1 - MaxCoeff;
   state = COEFF_CODE;
   for(cnt = 0; cnt < sip_siz; cnt++) {
      pix = *(sip + cnt);
      switch(state) {

         case COEFF_CODE:   /* for runs of zeros */
            if(pix == 0) {
               state = RUN_CODE;
               rcnt = 1;
               break;
            }
            if(pix > MaxCoeff) { 
               if(pix > 255)
                  counts[103]++; /* 16bit pos esc */
               else
                  counts[101]++; /* 8bit pos esc */
            }
            else if (pix < LoMaxCoeff) {
               if(pix < -255)
                  counts[104]++; /* 16bit neg esc */
               else
                  counts[102]++; /* 8bit neg esc */
            }
            else
               counts[pix+180]++; /* within table */
            break;

         case RUN_CODE:  /* get length of zero run */
            if(pix == 0  &&  rcnt < 0xFFFF) {
               ++rcnt;
               break;
            }
               /* limit rcnt to avoid EOF problem in bitio.c */
            if(rcnt <= MaxZRun)
               counts[rcnt]++;  /** log zero run length **/
            else if(rcnt <= 0xFF)
               counts[105]++;
            else if(rcnt <= 0xFFFF)
               counts[106]++; /* 16bit zrun esc */
            else {
//               fprintf(stdf,
  //             "ERROR: count_block : Zrun to long in count block.\n");
               return(-49);
            }

            if(pix != 0) {
               if(pix > MaxCoeff) { /** log current pix **/
                  if(pix > 255)
                     counts[103]++; /* 16bit pos esc */
                  else
                     counts[101]++; /* 8bit pos esc */
               }
               else if(pix < LoMaxCoeff) {
                  if(pix < -255)
                     counts[104]++; /* 16bit neg esc */
                  else
                     counts[102]++; /* 8bit neg esc */
               }
               else
                  counts[pix+180]++; /* within table */
               state = COEFF_CODE;
            }
            else {
               rcnt = 1;
               state = RUN_CODE;
            }
            break;
      }
   }
   if(state == RUN_CODE){ /** log zero run length **/
      if(rcnt <= MaxZRun)
         counts[rcnt]++;
      else if(rcnt <= 0xFF)
         counts[105]++;
      else if(rcnt <= 0xFFFF)
         counts[106]++; /* 16bit zrun esc */
      else {
//         fprintf(stdf,
  //       "ERROR: count_block : Zrun to long in count block.\n");
         return(-50);
      }
   }

   *ocounts = counts;
   return(0);
}

/*******************************************************/
/* Read a byte of data from the input buffer.          */
/*******************************************************/
int Wsq::getc_byte(
   unsigned char *ochar_dat,  /* pointer to returned byte       */
   unsigned char **cbufptr,   /* pointer to next byte in buffer */
   unsigned char *ebufptr)    /* pointer to end of buffer       */
{
   /* If at end of buffer ... */
   if((*cbufptr) >= ebufptr){
//      fprintf(stdf, "ERROR : getc_byte : premature End Of Buffer\n");
      return(-39);
   }

   /* Assign next byte in buffer. */
   *ochar_dat = **cbufptr;
   /* Bump buffer pointer. */
   (*cbufptr)++;

   return(0);
}

/*******************************************************/
/* Read a sequence of data bytes from the input buffer.*/
/*******************************************************/
int Wsq::getc_bytes(
   unsigned char **ochar_dat, /* pointer to returned bytes      */
   const int ilen,            /* number of bytes to be returned */
   unsigned char **cbufptr,   /* pointer to next byte in buffer */
   unsigned char *ebufptr)    /* pointer to end of buffer       */
{
   /* If at end of buffer ... */
   if((*cbufptr) >= (ebufptr-(ilen-1))){
//      fprintf(stdf, "ERROR : getc_bytes : premature End Of Buffer\n");
      return(-40);
   }

   /* Copy sequence of bytes from buffer. */
   memcpy(*ochar_dat, *cbufptr, ilen);

   /* Bump buffer pointer. */
   (*cbufptr)+=ilen;

   return(0);
}

/***********************************************/
/* Stores a byte of data to the output buffer. */
/***********************************************/
int Wsq::putc_byte(
   const unsigned char idata, /* input byte */
   unsigned char *odata,      /* output buffer of bytes */
   const int oalloc,          /* allocated size of output buffer */
   int *olen)                 /* filled length of output buffer  */
{
   /* olen points to next position in odata */
   /* If output buffer is out of space ...  */
   if((*olen) >= oalloc){
//      fprintf(stdf,
  //    "ERROR : putc_byte : buffer overlow : alloc = %d, request = %d\n",
    //  oalloc, *olen);
      return(-32);
   }

   *(odata+(*olen)) = idata;
   (*olen)++;

   return(0);
}

/**********************************************************************/
/* Stores a vector of bytes of specified length to the output buffer. */
/**********************************************************************/
int Wsq::putc_bytes(
   unsigned char *idata,  /* input buffer of bytes           */
   const int ilen,        /* bytes to be copied              */
   unsigned char *odata,  /* output buffer of bytes          */
   const int oalloc,      /* allocated size of output buffer */
   int *olen)             /* filled length of output buffer  */
{
   /* olen points to next position in odata */
   /* If output buffer is out of space ...  */
   if(((*olen)+ilen) > oalloc){
//      fprintf(stdf,
  //    "ERROR : putc_bytes : buffer overlow : alloc = %d, request = %d\n",
    //  oalloc, (*olen)+ilen);
      return(-33);
   }

   memcpy(odata+(*olen), idata, ilen);
   (*olen) += ilen;

   return(0);
}

/************************************************************/
/* Routine to read an unsigned short from the input buffer. */
/************************************************************/
int Wsq::getc_ushort(
   unsigned short *oshrt_dat,  /* pointer to returned unsigned short */
   unsigned char  **cbufptr,   /* pointer to next byte in buffer     */
   unsigned char  *ebufptr)    /* pointer to end of buffer           */
{
   int ret;
   unsigned short shrt_dat;
   unsigned char  *cptr;

   cptr = (unsigned char *)(&shrt_dat);

   if(ret = getc_bytes(&cptr, sizeof(unsigned short), cbufptr, ebufptr))
      return(ret);

#ifdef __i386__
   swap_short_bytes(shrt_dat);
#endif

   *oshrt_dat = shrt_dat;
   return(0);
}

/***************************************************************/
/* This routine stores an unsigned short to the output buffer. */
/***************************************************************/
int Wsq::putc_ushort(
   unsigned short ishort,     /* input unsigned short     */
   unsigned char *odata,      /* output byte buffer       */
   const int oalloc,          /* allocated size of buffer */
   int *olen)                 /* filled length of buffer  */
{
   int ret;
   unsigned char *cptr;

#ifdef __i386__
   swap_short_bytes(ishort);
#endif

   cptr = (unsigned char *)(&ishort);

   if(ret = putc_bytes(cptr, sizeof(unsigned short), odata, oalloc, olen))
      return(ret);

   return(0);
}

/**********************************************************/
/* Routine to read an unsigned int from the input buffer. */
/**********************************************************/
int Wsq::getc_uint(
   unsigned int  *oint_dat,  /* pointer to returned unsigned int */
   unsigned char **cbufptr,  /* pointer to next byte in buffer   */
   unsigned char *ebufptr)   /* pointer to end of buffer         */
{
   int ret;
   unsigned int int_dat;
   unsigned char  *cptr;

   cptr = (unsigned char *)(&int_dat);
   if(ret = getc_bytes(&cptr, sizeof(unsigned int), cbufptr, ebufptr))
      return(ret);

#ifdef __i386__
   swap_int_bytes(int_dat);
#endif

   *oint_dat = int_dat;
   return(0);
}

/****************************************************/
/* Stores an unsigned integer to the output buffer. */
/****************************************************/
int Wsq::putc_uint(
   unsigned int iint,        /* input unsigned int       */
   unsigned char *odata,     /* output byte buffer       */
   const int oalloc,         /* allocated size of buffer */
   int *olen)                /* filled length of buffer  */
{
   int ret;
   unsigned char *cptr;

#ifdef __i386__
   swap_int_bytes(iint);
#endif

   cptr = (unsigned char *)(&iint);

   if(ret = putc_bytes(cptr, sizeof(unsigned int), odata, oalloc, olen))
      return(ret);

   return(0);
}

/*******************************************************/
/*Routine to write "compressed" bits to output buffer. */
/*******************************************************/
void Wsq::write_bits(
   unsigned char **outbuf,    /* output data buffer                          */
   const unsigned short code, /* info to write into buffer                   */
   const short size,          /* numbers bits of code to write into buffer   */
   int *outbit,               /* current bit location in out buffer byte     */
   unsigned char *bits,       /* byte to write to output buffer              */
   int *bytes)                /* count of number bytes written to the buffer */
{
   short num;

   num = size;

   for(--num; num >= 0; num--) {
      *bits <<= 1;
      *bits |= ((unsigned char)((code >> num) & 0x0001));

      if(--(*outbit) < 0) {
         **outbuf = *bits;
         (*outbuf)++;
         if(*bits == 0xFF) {
            **outbuf = 0;
            (*outbuf)++;
            (*bytes)++;
         }
         (*bytes)++;
         *outbit = 7;
         *bits = 0;
      }
   }
   return;
}


/*********************************************/
/* Routine to "flush" left over bits in last */
/* byte after compressing a block.           */
/*********************************************/
void Wsq::flush_bits(
   unsigned char **outbuf, /* output data buffer */
   int *outbit,            /* current bit location in out buffer byte */
   unsigned char *bits,    /* byte to write to output buffer */
   int *bytes)             /* count of number bytes written to the buffer */
{
   int cnt;              /* temp counter */

   if(*outbit != 7) {
      for(cnt = *outbit; cnt >= 0; cnt--) {
         *bits <<= 1;
         *bits |= 0x01;
      }

      **outbuf = *bits;
      (*outbuf)++;
      if(*bits == 0xFF) {
         *bits = 0;
         **outbuf = 0;
         (*outbuf)++;
         (*bytes)++;
      }
      (*bytes)++;
      *outbit = 7;
      *bits = 0;
   }
   return;
}

/***********************************************************************
      LIBRARY: WSQ - Grayscale Image Compression

      FILE:    HUFF.C
      AUTHORS: Craig Watson
               Michael Garris
      DATE:    12/02/1999

      Checks that huffman codes are WSQ compliant. The specification
      does not allow for an all 1's code in the code table.

      ROUTINES:
#cat: check_huffcodes_wsq - Checks for an all 1's code in the code table.

***********************************************************************/

int Wsq::check_huffcodes_wsq(HUFFCODE *hufftable, int last_size)
{
   int i, k;
   int all_ones;

   for(i = 0; i < last_size; i++){
	   all_ones = 1;
	   for(k = 0; (k < (hufftable+i)->size) && all_ones; k++)
		   all_ones = (all_ones && (((hufftable+i)->code >> k) & 0x0001));
	   if(all_ones) {
		   //         fprintf(stdf, "WARNING: A code in the hufftable contains an ");
		   //       fprintf(stdf, "all 1's code.\n         This image may still be ");
		   //     fprintf(stdf, "decodable.\n         It is not compliant with ");
		   //   fprintf(stdf, "the WSQ specification.\n");
		   return(-1);
	   }
   }
   return(0);
}

/*****************************************************/
/* Reads huffman table from compressed memory buffer */
/*****************************************************/
int Wsq::getc_huffman_table(unsigned char *otable_id, unsigned char **ohuffbits,
                       unsigned char **ohuffvalues, const int max_huffcounts,
                       unsigned char **cbufptr, unsigned char *ebufptr,
                       const int read_table_len, int *bytes_left)
{
   int ret, i;
   unsigned short table_len = 0;
   unsigned char table_id;
   unsigned char *huffbits, *huffvalues;
   unsigned short num_hufvals;

   if(debug > 0)
      fprintf(stdf, "Start reading huffman table.\n");

   /* table_len */
   if(read_table_len){
      if(ret = getc_ushort(&table_len, cbufptr, ebufptr))
         return(ret);
      *bytes_left = table_len - 2;
   }

   /* If no bytes left ... */ 
   if((*bytes_left) <= 0){
      fprintf(stdf, "ERROR : getc_huffman_table : ");
      fprintf(stdf, "no huffman table bytes remaining\n");
      return(-2);
   }

   /* Table ID */
   if(ret = getc_byte(&table_id, cbufptr, ebufptr))
      return(ret);
   (*bytes_left)--;

   huffbits = (unsigned char *)calloc(MAX_HUFFBITS, sizeof(unsigned char));
   if(huffbits == (unsigned char *)NULL){
      fprintf(stdf,
              "ERROR : getc_huffman_table : calloc : huffbits\n");
      return(-3);
   }

   num_hufvals = 0;
   /* L1 ... L16 */
   for(i = 0; i < MAX_HUFFBITS; i++){
      if(ret = getc_byte(&(huffbits[i]), cbufptr, ebufptr)){
         free(huffbits);
         return(ret);
      }
      num_hufvals += huffbits[i];
   }
   (*bytes_left) -= MAX_HUFFBITS;

   if(num_hufvals > max_huffcounts+1){
      fprintf(stdf, "ERROR : getc_huffman_table : ");
      fprintf(stdf, "num_hufvals (%d) is larger", num_hufvals);
      fprintf(stdf, "than MAX_HUFFCOUNTS (%d)\n", max_huffcounts+1);
      free(huffbits);
      return(-4);
   }

   /* Could allocate only the amount needed ... then we wouldn't */
   /* need to pass MAX_HUFFCOUNTS. */
   huffvalues = (unsigned char *)calloc(max_huffcounts+1,
                                        sizeof(unsigned char));
   if(huffvalues == (unsigned char *)NULL){
      fprintf(stdf,
              "ERROR : getc_huffman_table : calloc : huffvalues\n");
      free(huffbits);
      return(-5);
   }

   /* V1,1 ... V16,16 */
   for(i = 0; i < num_hufvals; i++){
      if(ret = getc_byte(&(huffvalues[i]), cbufptr, ebufptr)){
         free(huffbits);
         free(huffvalues);
         return(ret);
      }
   }
   (*bytes_left) -= num_hufvals;

   if(debug > 1){
      fprintf(stdf, "Table Len = %d\n", table_len);
      fprintf(stdf, "Table ID = %d\n", table_id);
      for(i = 0; i < MAX_HUFFBITS; i++)
         fprintf(stdf, "bits[%d] = %d\n", i, huffbits[i]);
      for(i = 0; i < num_hufvals; i++)
         fprintf(stdf, "values[%d] = %d\n", i, huffvalues[i]);
   }
   
   if(debug > 0)
      fprintf(stdf, "Finished reading huffman table.\n");

   *otable_id = table_id;
   *ohuffbits = huffbits;
   *ohuffvalues = huffvalues;

   return(0);
}

/********************************************************/
/* Writes huffman table to the compressed memory buffer */
/********************************************************/
int Wsq::putc_huffman_table(
   unsigned short *factor, 
   const unsigned char table_id,   /* huffman table indicator  */
   unsigned char *huffbits,      /* huffman table parameters */
   unsigned char *huffvalues,
   unsigned char *outbuf,        /* output byte buffer       */
   const int outalloc,   /* allocated size of buffer */
   int   *outlen)        /* filled length of buffer  */
{
   int i, ret;
   unsigned short table_len, values_offset;

   if(debug > 0)
      fprintf(stdf, "Start writing huffman table.\n");

   /* DHT */
//   if(ret = putc_ushort(marker, outbuf, outalloc, outlen))
  //    return(ret);

   /* "value(2) + table id(1) + bits(16)" */
   table_len = values_offset = 3 + MAX_HUFFBITS;
   for(i = 0; i < MAX_HUFFBITS; i++)
      table_len += huffbits[i];   /* values size */

   if(debug > 1){
      fprintf(stdf, "Table Len = %d\n", table_len);
      fprintf(stdf, "Table ID = %d\n", table_id);
      for(i = 0; i < MAX_HUFFBITS; i++)
         fprintf(stdf, "bits[%d] = %d\n", i, huffbits[i]);
      for(i = 0; i < table_len-values_offset; i++)
         fprintf(stdf, "values[%d] = %d\n", i, huffvalues[i]);
   }

   /* Table Len */
	if (*factor)
	{
		if(ret = putc_ushort(table_len, outbuf, outalloc, outlen))
      return(ret);
	}
	*factor = table_len;

   /* Table ID */
   if(ret = putc_byte(table_id, outbuf, outalloc, outlen))
      return(ret);

   /* Huffbits (MAX_HUFFBITS) */
   for(i = 0; i < MAX_HUFFBITS; i++){
      if(ret = putc_byte(huffbits[i], outbuf, outalloc, outlen))
         return(ret);
   }

   /* Huffvalues (MAX_HUFFCOUNTS) */
   for(i = 0; i < table_len-values_offset; i++){
      if(ret = putc_byte(huffvalues[i], outbuf, outalloc, outlen))
         return(ret);
   }

   if(debug > 0)
      fprintf(stdf, "Finished writing huffman table.\n\n");

   return(0);
}

/******************************************************************/
/*routine to optimize code sizes by frequency of difference values*/
/******************************************************************/
int Wsq::find_huff_sizes(int **ocodesize, int *freq, const int max_huffcounts)
{
	int *codesize;       /*codesizes for each category*/
	int *others;         /*pointer used to generate codesizes*/
	int value1;          /*smallest and next smallest frequency*/
	int value2;          /*of difference occurrence in the largest
								difference category*/
	int i;               /*increment variable*/


	codesize = (int *)calloc(max_huffcounts+1, sizeof(int));
	if(codesize == (int *)NULL){
		fprintf(stdf, "ERROR : find_huff_sizes : calloc : codesize\n");
		return(-2);
	}
	others = (int *)malloc((max_huffcounts+1) * sizeof(int));
	if(others == (int *)NULL){
		fprintf(stdf, "ERROR : find_huff_sizes : malloc : others\n");
		return(-3);
	}

	for (i = 0; i <= max_huffcounts; i++) 
		others[i] = -1;

	while(1) {

		find_least_freq(&value1, &value2, freq, max_huffcounts);

		if(value2 == -1) {
			free(others);
			if(debug > 2){
				for (i = 0; i <= max_huffcounts; i++) 
					fprintf(stdf, "codesize[%d] = %d\n", i, codesize[i]);
			}
			break;
		}

		freq[value1] += freq[value2];
		freq[value2] = 0;

		codesize[value1]++;
		while(others[value1] != -1) {
			value1 = others[value1];
			codesize[value1]++;
		}
		others[value1] = value2;
		codesize[value2]++;

		while(others[value2] != -1) {
			value2 = others[value2];
			codesize[value2]++;
		}
	}

	*ocodesize = codesize;
	return(0);
}

/***********************************************************************/
/*routine to find the largest difference with the least frequency value*/
/***********************************************************************/
void Wsq::find_least_freq(int *value1, int *value2, int *freq,
                     const int max_huffcounts)
{
	int i;               /*increment variable*/
	int code_temp;       /*store code*/
	int value_temp;      /*store size*/
	int code2 = 0;           /*next smallest frequency in largest diff category*/
	int code1 = 0;           /*smallest frequency in largest difference category*/
	int set = 1;         /*flag first two non-zero frequency values*/

	*value1 = -1;
	*value2 = -1;

	for(i = 0; i <= max_huffcounts; i++) {
		if(freq[i] == 0)
			continue;
		if(set == 1) {
			code1 = freq[i];
			*value1 = i;
			set++;
			continue;
		}
		if(set == 2) {
			code2 = freq[i];
			*value2 = i;
			set++;
		}
		code_temp = freq[i];
		value_temp = i;
		if(code1 < code_temp && code2 < code_temp)
			continue;
		if((code_temp < code1) || (code_temp == code1 && value_temp > *value1)) {
			code2 = code1;
			*value2 = *value1;
			code1 = code_temp;
			*value1 = value_temp;
			continue;
		}
		if((code_temp < code2) || (code_temp == code2 && value_temp > *value2)) {
			code2 = code_temp;
			*value2 = value_temp;
		}
	}
}

/**********************************************/
/*routine to find number of codes of each size*/
/**********************************************/
int Wsq::find_num_huff_sizes(unsigned char **obits, int *adjust, int *codesize,
                        const int max_huffcounts)
{
   unsigned char *bits;    /*defines number of codes for each size*/
   int i;          /*increment variable*/

   *adjust = 0;

   /* Allocate 2X desired number of bits due to possible codesize. */
   bits = (unsigned char *)calloc((MAX_HUFFBITS<<1), sizeof(unsigned char));
   if(bits == (unsigned char *)NULL){
      fprintf(stdf, "ERROR : find_num_huff_sizes : calloc : bits\n");
      return(-2);
   }

   for(i = 0; i < max_huffcounts; i++) {
      if(codesize[i] != 0)
	 bits[(codesize[i] - 1)]++;
         if(codesize[i] > MAX_HUFFBITS)
            *adjust = 1;
   }

   if(debug > 2){
      for(i = 0; i < MAX_HUFFBITS<<1; i++)
	 fprintf(stdf, "bits[%d] = %d\n", i, bits[i]);
      fprintf(stdf, "ADJUST = %d\n", *adjust);
   }

   *obits = bits;
   return(0);
}

/****************************************************************/
/*routine to insure that no huffman code size is greater than 16*/
/****************************************************************/
int Wsq::sort_huffbits(unsigned char *bits)
{
   int i, j;
   int l1, l2, l3;
   short *tbits;

   l3 = MAX_HUFFBITS<<1;       /* 32 */
   l1 = l3 - 1;                /* 31 */
   l2 = MAX_HUFFBITS - 1;      /* 15 */

   tbits = (short *)malloc(l3*sizeof(short));
   if(tbits == (short *)NULL){
      fprintf(stdf, "ERROR : sort_huffbits : malloc : tbits\n");
      return(-2);
   }

   for(i = 0; i < MAX_HUFFBITS<<1; i++) tbits[i] = bits[i];

   for(i = l1; i > l2; i--) 
   {
	   while(tbits[i] > 0) 
	   {
		   j = i - 2;
		   while(tbits[j] == 0) j--;
		   tbits[i] -= 2;
		   tbits[i - 1] += 1;
		   tbits[j + 1] += 2;
		   tbits[j] -= 1;
	   }
	   tbits[i] = 0;
   }

   while(tbits[i] == 0)
      i--;

   tbits[i] -= 1;

   for(i = 0; i < MAX_HUFFBITS<<1; i++)
      bits[i] = tbits[i];
   free(tbits);

   for(i = MAX_HUFFBITS; i < l3; i++){
      if(bits[i] > 0){
         fprintf(stdf,
            "ERROR : sort_huffbits : Code length of %d is greater than 16.\n",
            i);
         return(-3);
      }
   }

   if(debug > 1){
      fprintf(stdf, "Huffbits after sorting.\n");
      for(i = 0; i < MAX_HUFFBITS<<1; i++)
         fprintf(stdf,"sort_bits[%d] = %d\n", i, bits[i]);
   }

   return(0);
}

/****************************************/
/*routine to sort the huffman code sizes*/
/****************************************/
int Wsq::sort_code_sizes(unsigned char **ovalues, int *codesize,
         const int max_huffcounts)
{
   unsigned char *values;      /*defines order of huffman codelengths in
                         relation to the code sizes*/
   int i, i2 = 0, i3;  /*increment variables*/


   values = (unsigned char *)calloc(max_huffcounts+1, sizeof(unsigned char));
   if(values == (unsigned char *)NULL){
      fprintf(stdf, "ERROR : sort_code_sizes : calloc : value\n");
      return(-2);
   }

   for(i = 1; i <= (MAX_HUFFBITS<<1); i++) {
      for(i3 = 0; i3 < max_huffcounts; i3++) {
	 if(codesize[i3] == i) {
	    values[i2] = i3;
	    i2++;
	 }
      }
   }

   if(debug > 2){
      for(i = 0; i <= max_huffcounts; i++)
	 fprintf(stdf, "values[%d] = %d\n", i, values[i]);
   }

   *ovalues = values;
   return(0);
}

/*****************************************/
/*routine to sort huffman codes and sizes*/
/*****************************************/
int Wsq::build_huffcode_table(HUFFCODE **ohuffcode_table,
          HUFFCODE *in_huffcode_table, const int last_size,
          unsigned char *values, const int max_huffcounts)
{
   int size;      /*huffman code size variable*/
   HUFFCODE *new_huffcode_table; /*pointer to a huffman code structure*/

   new_huffcode_table = (HUFFCODE *)calloc(max_huffcounts+1, sizeof(HUFFCODE));
   if(new_huffcode_table == (HUFFCODE *)NULL){
      fprintf(stdf,
      "ERROR : build_huffcode_table : calloc : new_huffcode_table\n");
      return(-2);
   }

   for(size = 0; size < last_size; size++) {
      (new_huffcode_table+values[size])->code = (in_huffcode_table+size)->code;
      (new_huffcode_table+values[size])->size = (in_huffcode_table+size)->size;
   }

   if(debug > 3){
      for(size = 0; size <= max_huffcounts; size++) {
         fprintf(stdf, "huff_size[%d] = %d\n", size,
                 new_huffcode_table[size].size);
         fprintf(stdf, "huff_code[%d] = %d\n", size,
                 new_huffcode_table[size].code);
      }
   }

   *ohuffcode_table = new_huffcode_table;
   return(0);
}

/**************************************************************************/
/*This routine defines the huffman code sizes for each difference category*/
/**************************************************************************/
int Wsq::build_huffsizes(HUFFCODE **ohuffcode_table, int *temp_size,
                    unsigned char *huffbits, const int max_huffcounts)
{
   HUFFCODE *huffcode_table;    /*table of huffman codes and sizes*/
   int code_size;               /*code sizes*/
   int number_of_codes = 1;     /*the number codes for a given code size*/

   huffcode_table = (HUFFCODE *)calloc(max_huffcounts+1, sizeof(HUFFCODE));
   if(huffcode_table == (HUFFCODE *)NULL){
      fprintf(stdf, "ERROR : build_huffsizes : calloc : huffcode_table\n");
      return(-2);
   }

   *temp_size = 0;

   for(code_size = 1; code_size <= MAX_HUFFBITS; code_size++) {
      while(number_of_codes <= huffbits[code_size - 1]) {
	 (huffcode_table + *temp_size)->size = code_size;
	 (*temp_size)++;
	 number_of_codes++;
      }
      number_of_codes = 1;
   }
   (huffcode_table+(*temp_size))->size = 0;

   if(debug > 2){
      int ii;
      fprintf(stdf, "In build_huffsizes:\n");
      for(ii = 0; ii < max_huffcounts+1; ii++)
         fprintf(stdf, "hf_sz[%d] = %d\n", ii, huffcode_table[ii].size);
      fflush(stdf);
   }

   *ohuffcode_table = huffcode_table;
   return(0);
}

/****************************************************************************/
/*This routine defines the huffman codes needed for each difference category*/
/****************************************************************************/
void Wsq::build_huffcodes(HUFFCODE *huffcode_table)
{
   int pointer = 0;                     /*pointer to code word information*/
   unsigned short temp_code = 0;        /*used to construct code word*/
   short  temp_size;                    /*used to construct code size*/

   temp_size = huffcode_table->size;
   if((huffcode_table+pointer)->size == 0)
      return;

   do {
	   do {
		   (huffcode_table+pointer)->code = temp_code;
		   temp_code++;
		   pointer++;
	   } while((huffcode_table+pointer)->size == temp_size);

	   if((huffcode_table+pointer)->size == 0)
		   return;

	   do {
		   temp_code <<= 1;
		   temp_size++;
	   } while((huffcode_table+pointer)->size != temp_size);
   } while((huffcode_table+pointer)->size == temp_size);
}

/*********************************************/
/*routine to generate tables needed to decode*/
/*********************************************/
void Wsq::gen_decode_table(HUFFCODE *huffcode_table,
            int *maxcode, int *mincode, int *valptr, unsigned char *huffbits)
{
   int i, i2 = 0;                   /*increment variables*/

   for(i = 0; i <= MAX_HUFFBITS; i++) {
      maxcode[i] = 0;
      mincode[i] = 0;
      valptr[i] = 0;
   }

   for(i = 1; i <= MAX_HUFFBITS; i++) {
	   if(huffbits[i-1] == 0) {
		   maxcode[i] = -1;
		   continue;
	   }
	   valptr[i] = i2;
	   mincode[i] = (huffcode_table + i2)->code;
	   i2 = i2 + huffbits[i - 1] - 1;
	   maxcode[i] = (huffcode_table + i2)->code;
	   i2++;
   }
   //if (debug)
   //{
	  // for (int i=0; i<=MAX_HUFFBITS; i++)
		 //  fprintf(stdf,"min=%d max=%d valptr=%d\n",mincode[i],maxcode[i],valptr[i]);
   //}
}

int Wsq::wsqCheckSum(unsigned char *src,  int lensrc)
{
   int sum=0;
   for(int i=0; i<lensrc; i+=10) sum+=(*(src+i)+(i<<1));
   return sum;
}

float Wsq::GetBitRate(float press)
{
  float btr[49],cmpr[49],bit_rate;

  btr[48] = (float)0.01;        cmpr[48] = (float)295.95;
  btr[47] = (float)0.02;        cmpr[47] = (float)233.88;
  btr[46] = (float)0.03;        cmpr[46] = (float)194.00;
  btr[45] = (float)0.04;        cmpr[45] = (float)173.16;
  btr[44] = (float)0.05;        cmpr[44] = (float)144.15;
  btr[43] = (float)0.10;        cmpr[43] = (float)80.20;
  btr[42] = (float)0.15;        cmpr[42] = (float)55.73;
  btr[41] = (float)0.20;        cmpr[41] = (float)45.10;
  btr[40] = (float)0.25;        cmpr[40] = (float)40.73;
  btr[39] = (float)0.30;        cmpr[39] = (float)24.15;
  btr[38] = (float)0.35;        cmpr[38] = (float)21.57;
  btr[37] = (float)0.40;        cmpr[37] = (float)21.57;
  btr[36] = (float)0.45;        cmpr[36] = (float)18.81;
  btr[35] = (float)0.50;        cmpr[35] = (float)18.04;
  btr[34] = (float)0.55;        cmpr[34] = (float)16.76;
  btr[33] = (float)0.60;        cmpr[33] = (float)17.35;
  btr[32] = (float)0.65;        cmpr[32] = (float)17.91;
  btr[31] = (float)0.70;        cmpr[31] = (float)16.24;
  btr[30] = (float)0.75;        cmpr[30] = (float)16.15;
  btr[29] = (float)0.80;        cmpr[29] = (float)16.45;
  btr[28] = (float)0.85;        cmpr[28] = (float)15.39;
  btr[27] = (float)0.90;        cmpr[27] = (float)15.09;
  btr[26] = (float)0.95;        cmpr[26] = (float)14.62;
  btr[25] = (float)1.00;        cmpr[25] = (float)13.93;
  btr[24] = (float)1.05;        cmpr[24] = (float)13.51;
  btr[23] = (float)1.10;        cmpr[23] = (float)12.83;
  btr[22] = (float)1.15;        cmpr[22] = (float)12.31;
  btr[21] = (float)1.20;        cmpr[21] = (float)11.68;
  btr[20] = (float)1.25;        cmpr[20] = (float)11.11;
  btr[19] = (float)1.30;        cmpr[19] = (float)10.57;
  btr[18] = (float)1.35;        cmpr[18] = (float)10.05;
  btr[17] = (float)1.40;        cmpr[17] = (float)9.57;
  btr[16] = (float)1.45;        cmpr[16] = (float)9.08;
  btr[15] = (float)1.50;        cmpr[15] = (float)8.64;
  btr[14] = (float)1.55;        cmpr[14] = (float)8.21;
  btr[13] = (float)1.60;        cmpr[13] = (float)7.79;
  btr[12] = (float)1.65;        cmpr[12] = (float)7.47;
  btr[11] = (float)1.70;        cmpr[11] = (float)7.18;
  btr[10] = (float)1.75;        cmpr[10] = (float)6.90;
  btr[ 9] = (float)1.80;        cmpr[ 9] = (float)6.64;
  btr[ 8] = (float)1.85;        cmpr[ 8] = (float)6.40;
  btr[ 7] = (float)1.90;        cmpr[ 7] = (float)6.17;
  btr[ 6] = (float)1.95;        cmpr[ 6] = (float)5.96;
  btr[ 5] = (float)2.00;        cmpr[ 5] = (float)5.77;
  btr[ 4] = (float)2.25;        cmpr[ 4] = (float)4.92;
  btr[ 3] = (float)2.50;        cmpr[ 3] = (float)4.24;
  btr[ 2] = (float)3.00;        cmpr[ 2] = (float)4.39;
  btr[ 1] = (float)4.00;        cmpr[ 1] = (float)2.25;
  btr[ 0] = (float)5.00;        cmpr[ 0] = (float)1.66;

	if( press<cmpr[ 0] )  bit_rate = btr[ 0];
	if( press>cmpr[48] )  bit_rate = btr[48];
	if( press>=cmpr[ 0] && press<=cmpr[48] )
	{
		for(int i=1; i<49; i++)
		{
			if( press>=cmpr[i-1] && press<=cmpr[i] )
			{
				float alf = (press-cmpr[i-1])/(cmpr[i]-cmpr[i-1]);
				bit_rate = btr[i-1]+alf*(btr[i]-btr[i-1]);
				break;
			}
		}
	}
	return bit_rate;
}
 
int Wsq::encode(unsigned char *idata, unsigned char *odata, int width, int height, float ratio, int chksum)
//int nist_encode(unsigned char *idata, unsigned char *odata, int width, int height, float ratio, int chksum)
{
   int ret; 
   char *outext="wsq";
   int rawflag=1;             /* input image flag: 0 == Raw, 1 == IHead */
   int depth=8, ppi=500;
   int olen;                 /* Number of bytes in output data. */
//   char *comment_text="XXXXX XXX.";
   char *comment_text="Sonda Lab.";
	float r_bitrate;

	if (ratio==0) return ERROR_BITRATE;
//	if (ratio==0 || ratio<-1) return ERROR_BITRATE;

	if (ratio<0)  
		r_bitrate = -ratio;
	else
		r_bitrate = GetBitRate(ratio);
#ifdef _WINDOWS
   if (debug>0) 
   {
	   int i=1;
	   char sfile[100];
	   do {
		   sprintf(sfile,"stdfin%02d.dbg",i++);
	   } while (!access(sfile,00));
	   stdf = fopen(sfile,"w");
	   if (!stdf) return -1;
   }
#endif
   if (wsqCheckSum(idata,(width*height)) != chksum )  return ERROR_CHECK_SUM;

//	for (int i=0;i<3; i++)
	for (int i=0;i<7; i++)
	{
		ret = wsq_encode_mem(odata,&olen,r_bitrate,idata,width,height,depth,ppi,comment_text);
		if (!ret) break;
		if (ret==-86 || ret==-87) 
			r_bitrate*=(float)1.08;
//			r_bitrate*=(float)1.05;
   }
	if (ret) return ret;

   if(debug > 0)
      fprintf(stdf, "Image data encoded, compressed byte length = %d\n",olen);

   if (debug>0) fclose(stdf);
   return olen;
}

//int nist_decode(unsigned char *odata, unsigned char *idata, int ilen, int *cols, int *rows, int chksum)
int Wsq::decode(unsigned char *odata, unsigned char *idata, int ilen, int *cols, int *rows, int chksum)
{
   int ret;
   int width, height;
   int depth, ppi;

#ifdef _WINDOWS
   if (debug>0)
   {
	   int i=1;
	   char sfile[100];
	   do {
		   sprintf(sfile,"stdfout%02d.dbg",i++);
	   } while (!access(sfile,00));
	   stdf = fopen(sfile,"w");
   if (!stdf) return -1;
   }
#endif
	code = 0;
	code2 = 0;
	ret = wsq_decode_mem(odata, &width, &height, &depth, &ppi,chksum, idata, ilen);
	*cols = width;
	*rows = height; 
   if(debug > 1) 
      fprintf(stdf, "Image pixmap constructed\n");

   if (debug>0) fclose(stdf);
   return ret;
}

int Wsq::putc_comment(
   const unsigned short marker,
   unsigned char *comment,        /* comment */
   const int cs,      /* comment size */
   unsigned char *odata,      /* output byte buffer       */
   const int oalloc,  /* allocated size of buffer */
   int   *olen)       /* filled length of buffer  */
{
   int ret, i;
   unsigned short hdr_size;              /* header size */

   if(debug > 0)
      fprintf(stdf, "Writing Comment Field to Buffer.\n");

   if(ret = putc_ushort(marker, odata, oalloc, olen))
      return(ret);
   /* comment size */
   hdr_size = 3 + cs;
//   hdr_size = 2 + cs;
   if(ret = putc_ushort(hdr_size, odata, oalloc, olen))
      return(ret);
   for(i = 0; i < cs; i++)
      if(ret = putc_byte(comment[i], odata, oalloc, olen))
         return(ret);
	// insert check sum to comment
	unsigned char check_sum = 0; 
	putc_byte(check_sum, odata, oalloc, olen);

   if(debug > 0)
      fprintf(stdf, "Finished Writing Comment Field to Buffer.\n");

   return(0);
}

/******************************************************/
/* Routine to read in WSQ markers from memory buffer. */
/******************************************************/
int Wsq::getc_marker_wsq(
   unsigned short *omarker,  /* marker read */
   const int type,   /* type of markers that could be found */
   unsigned char **cbufptr,  /* current byte in input buffer */
   unsigned char *ebufptr)   /* end of input buffer */
{
   int ret;
   unsigned short marker;    /* WSQ marker */

   if(ret = getc_ushort(&marker, cbufptr, ebufptr))
      return(ret);

   switch(type){
   case SOI_WSQ:
      if(marker != SOI_WSQ) {
			if (debug)
         fprintf(stdf,
         "ERROR : getc_marker_wsq : No SOI marker. {%04X}\n", marker);
         return(-88);
      }
      break;
   case TBLS_N_SOF:
      if(marker != DTT_WSQ && marker != DQT_WSQ && marker != DHT_WSQ
         && marker != SOF_WSQ && marker != COM_WSQ) {
			if (debug)
         fprintf(stdf,
         "ERROR : getc_marker_wsq : No SOF, Table, or comment markers.\n");
         return(-89);
      }
      break;
   case TBLS_N_SOB:
      if(marker != DTT_WSQ && marker != DQT_WSQ && marker != DHT_WSQ
         && marker != SOB_WSQ && marker != COM_WSQ) {
			if (debug)
         fprintf(stdf,
         "ERROR : getc_marker_wsq : No SOB, Table, or comment markers.{%04X}\n",
                 marker);
         return(-90);
      }
      break;
   case ANY_WSQ:
      if((marker & 0xff00) != 0xff00){
			if (debug)
	fprintf(stdf,"ERROR : getc_marker_wsq : no marker found {%04X}\n",
                marker);
         return(-91);
      }
      break;
   default:
			if (debug)
      fprintf(stdf,
      "ERROR : getc_marker_wsq : Invalid marker -> {%4X}\n", marker);
      return(-92);
   }

   *omarker = marker;
   return(0);
}

/*******************************************************/
/* Routine to read specified table from memory buffer. */
/*******************************************************/
int Wsq::getc_table_wsq(
   unsigned short marker,         /* WSQ marker */
   DTT_TABLE *dtt_table,  /* transform table structure */
   DQT_TABLE *dqt_table,  /* quantization table structure */
   DHT_TABLE *dht_table,  /* huffman table structure */
   unsigned char **cbufptr,  /* current byte in input buffer */
   unsigned char *ebufptr)   /* end of input buffer */
{
   int ret;
   unsigned char *comment;

   switch(marker){
   case DTT_WSQ:
      if(ret = getc_transform_table(dtt_table, cbufptr, ebufptr))
         return(ret);
      break;
   case DQT_WSQ:
      if(ret = getc_quantization_table(dqt_table, cbufptr, ebufptr))
         return(ret);
      break;
   case DHT_WSQ:
      if(ret = getc_huffman_table_wsq(dht_table, cbufptr, ebufptr))
         return(ret);
      break;
   case COM_WSQ:
      if(ret = getc_comment(&comment, cbufptr, ebufptr))
         return(ret);
#ifdef PRINT_COMMENT
      fprintf(stdf, "COMMENT:\n%s\n\n", comment);
#endif
      free(comment);
      break;
   default:
			if (debug)
      fprintf(stdf,"ERROR: getc_table_wsq : Invalid table defined -> {%u}\n",
              marker);
      return(-93);
   }

   return(0);
}

/*********************************************************************/
/* Routine to read in transform table parameters from memory buffer. */
/*********************************************************************/
int Wsq::getc_transform_table(
   DTT_TABLE *dtt_table,  /* transform table structure */
   unsigned char **cbufptr,  /* current byte in input buffer */
   unsigned char *ebufptr)   /* end of input buffer */
{
   int ret;
   unsigned short hdr_size;              /* header size */
   double *a_lofilt, *a_hifilt;  /* unexpanded filter coefficients */
//   float *a_lofilt, *a_hifilt;  /* unexpanded filter coefficients */
   unsigned char a_size;                 /* size of unexpanded coefficients */
   unsigned int cnt, shrt_dat;           /* counter and temp short data */
   unsigned char scale, sign;            /* scaling and sign parameters */

   if(debug > 0)
      fprintf(stdf, "Reading transform table.\n");

   if(ret = getc_ushort(&hdr_size, cbufptr, ebufptr))
      return(ret);
   if(ret = getc_byte(&(dtt_table->hisz), cbufptr, ebufptr))
      return(ret);
   if(ret = getc_byte(&(dtt_table->losz), cbufptr, ebufptr))
      return(ret);


   if(debug > 2) {
      fprintf(stdf, "losize = %d\n", dtt_table->losz);
      fprintf(stdf, "hisize = %d\n", dtt_table->hisz);
   }

   /* Added 02-24-05 by MDG */
   /* If lofilt member previously allocated ... */
   if(dtt_table->lofilt != (double *)NULL){
      /* Deallocate the member prior to new allocation */
      free(dtt_table->lofilt);
      dtt_table->lofilt = (double *)NULL;
   }

   dtt_table->lofilt = (double *)calloc(dtt_table->losz,sizeof(double));
   if(dtt_table->lofilt == (double *)NULL) {
			if (debug)
      fprintf(stdf,
      "ERROR : getc_transform_table : calloc : lofilt\n");
      return(-94);
   }

   /* Added 02-24-05 by MDG */
   /* If hifilt member previously allocated ... */
   if(dtt_table->hifilt != (double *)NULL){
      /* Deallocate the member prior to new allocation */
      free(dtt_table->hifilt);
      dtt_table->hifilt = (double *)NULL;
   }

   dtt_table->hifilt = (double *)calloc(dtt_table->hisz,sizeof(double));
   if(dtt_table->hifilt == (double *)NULL) {
      free(dtt_table->lofilt);
			if (debug)
      fprintf(stdf,
      "ERROR : getc_transform_table : calloc : hifilt\n");
      return(-95);
   }

   if(dtt_table->hisz % 2)
      a_size = (dtt_table->hisz + 1) / 2;
   else
      a_size = dtt_table->hisz / 2;

   a_lofilt = (double *) calloc(a_size, sizeof(double));
   if(a_lofilt == (double *)NULL) {
      free(dtt_table->lofilt);
      free(dtt_table->hifilt);
			if (debug)
      fprintf(stdf,
      "ERROR : getc_transform_table : calloc : a_lofilt\n");
      return(-96);
   }

   a_size--;
   for(cnt = 0; cnt <= a_size; cnt++) {
      if(ret = getc_byte(&sign, cbufptr, ebufptr)){
         free(dtt_table->lofilt);
         free(dtt_table->hifilt);
         free(a_lofilt);
         return(ret);
      }
      if(ret = getc_byte(&scale, cbufptr, ebufptr)){
         free(dtt_table->lofilt);
         free(dtt_table->hifilt);
         free(a_lofilt);
         return(ret);
      }
      if(ret = getc_uint(&shrt_dat, cbufptr, ebufptr)){
         free(dtt_table->lofilt);
         free(dtt_table->hifilt);
         free(a_lofilt);
         return(ret);
      }
      a_lofilt[cnt] = (float)shrt_dat;
      while(scale > 0) {
         a_lofilt[cnt] /= 10.0;
         scale--; 
      }
      if(sign != 0)
         a_lofilt[cnt] *= -1.0;

      if(debug > 3)
         fprintf(stdf, "lofilt[%d] = %.15f\n", cnt, a_lofilt[cnt]);

      if(dtt_table->hisz % 2) {
         dtt_table->hifilt[cnt + a_size] = (double)((double)int_sign(cnt)
                                         * a_lofilt[cnt]);
         if(cnt > 0)
            dtt_table->hifilt[a_size - cnt] = dtt_table->hifilt[cnt + a_size];
      }
      else {
         dtt_table->hifilt[cnt + a_size + 1] = (double)((double)int_sign(cnt)
                                             * a_lofilt[cnt]);
         dtt_table->hifilt[a_size - cnt] = (double)-1.0 *
                                    dtt_table->hifilt[cnt + a_size + 1];
      }

   }
   free(a_lofilt);

   if(dtt_table->losz % 2)
      a_size = (dtt_table->losz + 1) / 2;
   else
      a_size = dtt_table->losz / 2;

   a_hifilt = (double *) calloc(a_size, sizeof(double));
   if(a_hifilt == (double *)NULL) {
      free(dtt_table->lofilt);
      free(dtt_table->hifilt);
			if (debug)
      fprintf(stdf,
      "ERROR : getc_transform_table : calloc : a_hifilt\n");
      return(-97);
   }

   a_size--;
   for(cnt = 0; cnt <= a_size; cnt++) {
      if(ret = getc_byte(&sign, cbufptr, ebufptr)){
         free(dtt_table->lofilt);
         free(dtt_table->hifilt);
         free(a_hifilt);
         return(ret);
      }
      if(ret = getc_byte(&scale, cbufptr, ebufptr)){
         free(dtt_table->lofilt);
         free(dtt_table->hifilt);
         free(a_hifilt);
         return(ret);
      }
      if(ret = getc_uint(&shrt_dat, cbufptr, ebufptr)){
         free(dtt_table->lofilt);
         free(dtt_table->hifilt);
         free(a_hifilt);
         return(ret);
      }
      a_hifilt[cnt] = (float)shrt_dat;
      while(scale > 0) {
         a_hifilt[cnt] /= 10.0;
         scale--;
      }
      if(sign != 0)
         a_hifilt[cnt] *= -1.0;

      if(debug > 2)
         fprintf(stdf, "hifilt[%d] = %.15f\n", cnt, a_hifilt[cnt]);

      if(dtt_table->losz % 2) {
         dtt_table->lofilt[cnt + a_size] = (float)((float)int_sign(cnt)
                                         * a_hifilt[cnt]);
         if(cnt > 0)
            dtt_table->lofilt[a_size - cnt] = dtt_table->lofilt[cnt + a_size];
      }
      else {
         dtt_table->lofilt[cnt + a_size + 1] = (float)((float)int_sign(cnt+1)
                                             * a_hifilt[cnt]);
         dtt_table->lofilt[a_size - cnt] = dtt_table->lofilt[cnt + a_size + 1];
      }


   }
   free(a_hifilt);

   dtt_table->lodef = 1;
   dtt_table->hidef = 1;

   if(debug > 0)
      fprintf(stdf, "Finished reading transform table.\n\n");

   return(0);
}

/************************************************/
/* Stores transform table to the output buffer. */
/************************************************/
int Wsq::putc_transform_table(
   double *lofilt,     /* filter coefficients      */
//   float *lofilt,     /* filter coefficients      */
   const int losz,
   double *hifilt,
//   float *hifilt,
   const int hisz,
   unsigned char *odata,      /* output byte buffer       */
   const int oalloc,  /* allocated size of buffer */
   int   *olen)       /* filled length of buffer  */
{
   int ret;
   int coef;           /* filter coefficient indicator */
   unsigned int int_dat;        /* temp variable */
   double dbl_tmp;       /* temp variable */
   char scale_ex, sign; /* exponent scaling and sign parameters */

   if(debug > 0)
      fprintf(stdf, "Writing transform table.\n");

   if(ret = putc_ushort(DTT_WSQ, odata, oalloc, olen))
      return(ret);
   /* table size */
   if(ret = putc_ushort(58, odata, oalloc, olen))
      return(ret);
   /* number analysis lowpass coefficients */
   if(ret = putc_byte(losz, odata, oalloc, olen))
      return(ret);
   /* number analysis highpass coefficients */
   if(ret = putc_byte(hisz, odata, oalloc, olen))
      return(ret);

   for(coef = (losz>>1); coef < losz; coef++) {
      dbl_tmp = lofilt[coef];
      if(dbl_tmp >= 0.0) {
         sign = 0;
      }
      else {
         sign = 1;
         dbl_tmp *= -1.0;
      }
      scale_ex = 0;
      if(dbl_tmp == 0.0)
         int_dat = 0;
      else if(dbl_tmp < 4294967295.0) {
         while(dbl_tmp < 4294967295.0) {
            scale_ex += 1;
            dbl_tmp *= 10.0;
         }
         scale_ex -= 1;
         int_dat = (unsigned int)sround_uint(dbl_tmp / 10.0);
      }
      else {
         dbl_tmp = lofilt[coef];
			if (debug)
         fprintf(stdf,
         "ERROR: putc_transform_table : lofilt[%d] to high at %f\n",
         coef, dbl_tmp);
         return(-82);
      }

      if(debug > 2) {
         fprintf(stdf, "lo[%d] = %u\n", coef, int_dat);
         fprintf(stdf, "lof[%d] = %0.15f\n", coef, lofilt[coef]);
      }

      if(ret = putc_byte(sign, odata, oalloc, olen))
         return(ret);
      if(ret = putc_byte(scale_ex, odata, oalloc, olen))
         return(ret);
      if(ret = putc_uint(int_dat, odata, oalloc, olen))
         return(ret);
   }

   for(coef = (hisz>>1); coef < hisz; coef++) {
      dbl_tmp = hifilt[coef];
      if(dbl_tmp >= 0.0) {
         sign = 0;
      }
      else {
         sign = 1;
         dbl_tmp *= -1.0;
      }
      scale_ex = 0;
      if(dbl_tmp == 0.0)
         int_dat = 0;
      else if(dbl_tmp < 4294967295.0) {
         while(dbl_tmp < 4294967295.0) {
            scale_ex += 1;
            dbl_tmp *= 10.0;
         }
         scale_ex -= 1;
         int_dat = (unsigned int)sround_uint(dbl_tmp / 10.0);
      }
      else {
         dbl_tmp = hifilt[coef];
			if (debug)
         fprintf(stdf,
         "ERROR: putc_transform_table : hifilt[%d] to high at %f\n",
         coef, dbl_tmp);
         return(-83);
      }

      if(debug > 2) {
         fprintf(stdf, "hi[%d] = %u\n", coef, int_dat);
         fprintf(stdf, "hif[%d] = %0.15f\n", coef, hifilt[coef]);
      }

      if(ret = putc_byte(sign, odata, oalloc, olen))
         return(ret);
      if(ret = putc_byte(scale_ex, odata, oalloc, olen))
         return(ret);
      if(ret = putc_uint(int_dat, odata, oalloc, olen))
         return(ret);
   }

   if(debug > 0)
      fprintf(stdf, "Finished writing transform table.\n\n");

   return(0);
}

/************************************************************************/
/* Routine to read in quantization table parameters from memory buffer. */
/************************************************************************/
int Wsq::getc_quantization_table(
   DQT_TABLE *dqt_table,  /* quatization table structure */
   unsigned char **cbufptr,  /* current byte in input buffer */
   unsigned char *ebufptr)   /* end of input buffer */
{
   int ret;
   unsigned short hdr_size;       /* header size */
   unsigned short cnt, shrt_dat;  /* counter and temp short data */
   unsigned char scale;           /* scaling parameter */

   if(debug > 0)
      fprintf(stdf, "Reading quantization table.\n");

   if(ret = getc_ushort(&hdr_size, cbufptr, ebufptr))
      return(ret);
   if(ret = getc_byte(&scale, cbufptr, ebufptr))
      return(ret);
   if(ret = getc_ushort(&shrt_dat, cbufptr, ebufptr))
      return(ret);
   dqt_table->bin_center = (double)shrt_dat;
//   dqt_table->bin_center = (float)shrt_dat;
   while(scale > 0) {
      dqt_table->bin_center /= 10.0;
      scale--;
   }

   for(cnt = 0; cnt < 64; cnt++) {
      if(ret = getc_byte(&scale, cbufptr, ebufptr))
         return(ret);
      if(ret = getc_ushort(&shrt_dat, cbufptr, ebufptr))
         return(ret);
      dqt_table->q_bin[cnt] = (double)shrt_dat;
//      dqt_table->q_bin[cnt] = (float)shrt_dat;
      while(scale > 0) {
         dqt_table->q_bin[cnt] /= 10.0;
         scale--;
      }
      if(ret = getc_byte(&scale, cbufptr, ebufptr))
         return(ret);
      if(ret = getc_ushort(&shrt_dat, cbufptr, ebufptr))
         return(ret);
      dqt_table->z_bin[cnt] = (double)shrt_dat;
//      dqt_table->z_bin[cnt] = (float)shrt_dat;
      while(scale > 0) {
         dqt_table->z_bin[cnt] /= 10.0;
         scale--;
      }

      if(debug > 2)
         fprintf(stdf, "q[%d] = %f :: z[%d] = %f\n",
         cnt, dqt_table->q_bin[cnt], cnt, dqt_table->z_bin[cnt]);

   }
   dqt_table->dqt_def = 1;

   if(debug > 0)
      fprintf(stdf, "Finished reading quantization table.\n\n");

   return(0);
}


/***************************************************/
/* Stores quantization table in the output buffer. */
/***************************************************/
int Wsq::putc_quantization_table(
   QUANT_VALS *quant_vals,   /* quantization parameters  */
   unsigned char *odata,             /* output byte buffer       */
   const int oalloc,         /* allocated size of buffer */
   int   *olen)              /* filled length of buffer  */
{
   int ret, sub;               /* subband indicators */
   char scale_ex, scale_ex2;   /* exponent scaling parameters */
   unsigned short shrt_dat, shrt_dat2; /* temp variables */
   double flt_tmp;              /* temp variable */
//   float flt_tmp;              /* temp variable */


   if(debug > 0)
      fprintf(stdf, "Writing quantization table.\n");

   if(ret = putc_ushort(DQT_WSQ, odata, oalloc, olen))
      return(ret);
   /* table size */
   if(ret = putc_ushort(389, odata, oalloc, olen))
      return(ret);
   /* exponent scaling value */
   if(ret = putc_byte(2, odata, oalloc, olen))
      return(ret);
   /* quantizer bin center parameter */
   if(ret = putc_ushort(44, odata, oalloc, olen))
      return(ret);

   for(sub = 0; sub < 64; sub ++) 
   {
	   if(sub >= 0 && sub < 60) 
	   {
		   if(quant_vals->qbss[sub] != 0.0) 
		   {
			   flt_tmp = quant_vals->qbss[sub];
			   scale_ex = 0;
			   if(flt_tmp < 65535) 
			   {
				   while(flt_tmp < 65535) 
				   {
					   scale_ex += 1;
					   flt_tmp *= 10.;
				   }
				   scale_ex -= 1;
				   shrt_dat = (unsigned short)sround(flt_tmp / 10.0);
			   }
			   else 
			   {
				   flt_tmp = quant_vals->qbss[sub];
				   if (debug)
					   fprintf(stdf,
					   "ERROR : putc_quantization_table : Q[%d] to high at %f\n",
					   sub, flt_tmp);
				   return(-86);
			   }

			   flt_tmp = quant_vals->qzbs[sub];
			   scale_ex2 = 0;
			   if(flt_tmp < 65535) 
			   {
				   while(flt_tmp < 65535) 
				   {
					   scale_ex2 += 1;
					   flt_tmp *= 10.;
				   }
				   scale_ex2 -= 1;
				   shrt_dat2 = (unsigned short)sround(flt_tmp / 10.0);
			   }
			   else 
			   {
				   flt_tmp = quant_vals->qzbs[sub];
				   if (debug)
					   fprintf(stdf,
					   "ERROR : putc_quantization_table : Z[%d] to high at %f\n",
					   sub, flt_tmp);
				   return(-87);
			   }
		   }
		   else {
			   scale_ex = 0;
			   scale_ex2 = 0;
			   shrt_dat = 0;
			   shrt_dat2 = 0;
		   }
	   }
	   else {
		   scale_ex = 0;
		   scale_ex2 = 0;
		   shrt_dat = 0;
		   shrt_dat2 = 0;
	   }

	   if(debug > 2) {
		   fprintf(stdf,
			   "qi[%d] = %d    ::  zi[%d] = %d\n", sub, shrt_dat, sub, shrt_dat2);
		   fprintf(stdf,
			   "q[%d] = %5.7f  ::  z[%d] = %5.7f\n", sub, quant_vals->qbss[sub],
			   sub, quant_vals->qzbs[sub]);
	   }

	   if(ret = putc_byte(scale_ex, odata, oalloc, olen))
		   return(ret);
	   if(ret = putc_ushort(shrt_dat, odata, oalloc, olen))
		   return(ret);
	   if(ret = putc_byte(scale_ex2, odata, oalloc, olen))
		   return(ret);
	   if(ret = putc_ushort(shrt_dat2, odata, oalloc, olen))
		   return(ret);
   }

   if(debug > 0)
      fprintf(stdf, "Finished writing quantization table.\n\n");

   return(0);
}

/*******************************************************************/
/* Routine to read in huffman table parameters from memory buffer. */
/*******************************************************************/
int Wsq::getc_huffman_table_wsq(
   DHT_TABLE *dht_table,  /* huffman table structure */
   unsigned char **cbufptr,  /* current byte in input buffer */
   unsigned char *ebufptr)   /* end of input buffer */
{
   int ret;
   unsigned char table_id;        /* huffman table indicator */
   unsigned char *huffbits;
   unsigned char *huffvalues;
   int bytes_left;

   /* First time, read table len. */
   if(ret = getc_huffman_table(&table_id, &huffbits, &huffvalues,
                               MAX_HUFFCOUNTS_WSQ, cbufptr, ebufptr,
                               READ_TABLE_LEN, &bytes_left))
      return(ret);

   /* Store table into global structure list. */
   memcpy((dht_table+table_id)->huffbits, huffbits, MAX_HUFFBITS);
   memcpy((dht_table+table_id)->huffvalues, huffvalues,
          MAX_HUFFCOUNTS_WSQ+1);
   (dht_table+table_id)->tabdef = 1;
   free(huffbits);
   free(huffvalues);

   while(bytes_left){
      /* Read next table without rading table len. */
      if(ret = getc_huffman_table(&table_id, &huffbits, &huffvalues,
                                  MAX_HUFFCOUNTS_WSQ, cbufptr, ebufptr,
                                  NO_READ_TABLE_LEN, &bytes_left))
         return(ret);

      /* If table is already defined ... */
      if((dht_table+table_id)->tabdef){
         free(huffbits);
         free(huffvalues);
			if (debug)
			{
         fprintf(stdf, "ERROR : getc_huffman_table_wsq : ");
         fprintf(stdf, "huffman table ID = %d already defined\n", table_id);
			}
         return(-2);
      }

      /* Store table into global structure list. */
      memcpy((dht_table+table_id)->huffbits, huffbits, MAX_HUFFBITS);
      memcpy((dht_table+table_id)->huffvalues, huffvalues,
             MAX_HUFFCOUNTS_WSQ+1);
      (dht_table+table_id)->tabdef = 1;
      free(huffbits);
      free(huffvalues);
   }

   return(0);
}

/******************************************************************/
/* Routine to read in frame header parameters from memory buffer. */
/******************************************************************/
int Wsq::getc_frame_header_wsq(
   FRM_HEADER_WSQ *frm_header,  /* frame header structure */
   unsigned char **cbufptr,  /* current byte in input buffer */
   unsigned char *ebufptr)   /* end of input buffer */
{
   int ret;
   unsigned short hdr_size, shrt_dat;  /* header size and data pointer */
   unsigned char scale;                /* exponent scaling parameter */

   if(debug > 0)
      fprintf(stdf, "Reading frame header.\n");

   if(ret = getc_ushort(&hdr_size, cbufptr, ebufptr))
      return(ret);
   if(ret = getc_byte(&(frm_header->black), cbufptr, ebufptr))
      return(ret);
   if(ret = getc_byte(&(frm_header->white), cbufptr, ebufptr))
      return(ret);
   if(ret = getc_ushort(&(frm_header->height), cbufptr, ebufptr))
      return(ret);
   if(ret = getc_ushort(&(frm_header->width), cbufptr, ebufptr))
      return(ret);
   if(ret = getc_byte(&scale, cbufptr, ebufptr))
      return(ret);
   if(ret = getc_ushort(&shrt_dat, cbufptr, ebufptr))
      return(ret);
   frm_header->m_shift = (float) shrt_dat;
   while(scale > 0) {
      frm_header->m_shift /= 10.0;
      scale--;
   }
   if(ret = getc_byte(&scale, cbufptr, ebufptr))
      return(ret);
   if(ret = getc_ushort(&shrt_dat, cbufptr, ebufptr))
      return(ret);
   frm_header->r_scale = (float) shrt_dat;
   while(scale > 0) {
      frm_header->r_scale /= 10.0;
      scale--;
   }

   if(ret = getc_byte(&(frm_header->wsq_encoder), cbufptr, ebufptr))
      return(ret);
   if(ret = getc_ushort(&(frm_header->software), cbufptr, ebufptr))
      return(ret);

   if(debug > 2) {
       fprintf(stdf, "black = %d :: white = %u\n",
               frm_header->black, frm_header->white);
       fprintf(stdf, "w = %d :: h = %d\n",
               frm_header->width, frm_header->height);
       fprintf(stdf, "m_shift = %f :: r_scale = %f\n",
               frm_header->m_shift,frm_header->r_scale);
       fprintf(stdf, "WSQ_encoder = %d\n",
                        frm_header->wsq_encoder);
       fprintf(stdf, "Software = %d\n", frm_header->software);
   }
   if(debug > 0)
      fprintf(stdf, "Finished reading frame header.\n\n");
   
   return(0);
}

/*********************************************/
/* Stores frame header to the output buffer. */
/*********************************************/
int Wsq::putc_frame_header_wsq(
   const int width,       /* image width              */
   const int height,      /* image height             */
   const float m_shift,   /* image shifting parameter */
   const float r_scale,   /* image scaling parameter  */
   unsigned char *odata,          /* output byte buffer       */
   const int oalloc,      /* allocated size of buffer */
   int   *olen)           /* filled length of buffer  */
{
   int ret;
   float flt_tmp;         /* temp variable */
   char scale_ex;         /* exponent scaling parameter */
   unsigned short shrt_dat;       /* temp variable */

   if(debug > 0)
      fprintf(stdf, "Writing frame header.\n");

   if(ret = putc_ushort(SOF_WSQ, odata, oalloc, olen))
      return(ret);
   /* size of frame header */
   if(ret = putc_ushort(17, odata, oalloc, olen))
      return(ret);
   /* black pixel */
   if(ret = putc_byte(0, odata, oalloc, olen))
      return(ret);
   /* white pixel */
   if(ret = putc_byte(255, odata, oalloc, olen))
      return(ret);
   if(ret = putc_ushort(height, odata, oalloc, olen))
      return(ret);
   if(ret = putc_ushort(width, odata, oalloc, olen))
      return(ret);

   if(debug > 2)
      fprintf(stdf,
              "m_shift = %f  :: r_scale = %f\n", m_shift, r_scale);

   flt_tmp = m_shift;
   scale_ex = 0;
   if(flt_tmp != 0.0) {
      while(flt_tmp < 65535) {
         scale_ex += 1;
         flt_tmp *= 10;
      }
      scale_ex -= 1;
      shrt_dat = (unsigned short)sround(flt_tmp / 10.0);
   }
   else
      shrt_dat = 0;
   if(ret = putc_byte(scale_ex, odata, oalloc, olen))
      return(ret);
   if(ret = putc_ushort(shrt_dat, odata, oalloc, olen))
      return(ret);

   flt_tmp = r_scale;
   scale_ex = 0;
   if(flt_tmp != 0.0) {
      while(flt_tmp < 65535) {
         scale_ex += 1;
         flt_tmp *= 10;
      }
      scale_ex -= 1;
      shrt_dat = (unsigned short)sround(flt_tmp / 10.0);
   }
   else
      shrt_dat = 0;
   if(ret = putc_byte(scale_ex, odata, oalloc, olen))
      return(ret);
   if(ret = putc_ushort(shrt_dat, odata, oalloc, olen))
      return(ret);
   if(ret = putc_byte(0x02, odata, oalloc, olen))  // Ev
//   if(ret = putc_byte(0, odata, oalloc, olen))
      return(ret);
//   if(ret = putc_ushort(0xb091, odata, oalloc, olen)) // SSoftware for 2.0
   if(ret = putc_ushort(0xb095, odata, oalloc, olen)) // SSoftware for 3.1
      return(ret);

   if(debug > 0)
      fprintf(stdf, "Finished writing frame header.\n\n");

   return(0);
}

/******************************************************************/
/* Routine to read in block header parameters from memory buffer. */
/******************************************************************/
int Wsq::getc_block_header(
   unsigned char *huff_table,   /* huffman table indicator */
   unsigned char **cbufptr,  /* current byte in input buffer */
   unsigned char *ebufptr)   /* end of input buffer */
{
   int ret;
   unsigned short hdr_size;     /* block header size */

   if(debug > 0)
      fprintf(stdf, "Reading block header.\n");

   if(ret = getc_ushort(&hdr_size, cbufptr, ebufptr))
      return(ret);
   if(ret = getc_byte(huff_table, cbufptr, ebufptr))
      return(ret);

   if(debug > 2)
      fprintf(stdf, "huff_table = %d\n", *huff_table);
   if(debug > 0)
      fprintf(stdf, "Finished reading block header.\n\n");

   return(0);
}

/*********************************************/
/* Stores block header to the output buffer. */
/*********************************************/
int Wsq::putc_block_header(
   const int table,   /* huffman table indicator  */
   unsigned char *odata,      /* output byte buffer       */
   const int oalloc,  /* allocated size of buffer */
   int   *olen)       /* filled length of buffer  */
{
   int ret;

   if(debug > 0)
      fprintf(stdf, "Writing block header.\n");

   if(ret = putc_ushort(SOB_WSQ, odata, oalloc, olen))
      return(ret);
   /* block header size */
   if(ret = putc_ushort(3, odata, oalloc, olen))
      return(ret);
   if(ret = putc_byte((unsigned char)table, odata, oalloc, olen))
      return(ret);

   if(debug > 0)
      fprintf(stdf, "Finished writing block header.\n\n");

   return(0);
}

int Wsq::getc_comment(
   unsigned char **ocomment,
   unsigned char **cbufptr,  /* current byte in input buffer */
   unsigned char *ebufptr)   /* end of input buffer */
{
   int ret, cs;
   unsigned short hdr_size;              /* header size */
   unsigned char *comment;

   if(debug > 0)
      fprintf(stdf, "Reading Comment Field.\n");

   if(ret = getc_ushort(&hdr_size, cbufptr, ebufptr))
      return(ret);

   /* cs = hdr_size - sizeof(length value) */
   cs = hdr_size - 2;

   /* Allocate including a possible NULL terminator. */
   comment = (unsigned char *)calloc(cs+1, sizeof(unsigned char));
   if(comment == (unsigned char *)NULL){
			if (debug)
      fprintf(stdf, "ERROR : getc_comment : malloc : comment\n");
      return(-2);
   }

   /* Read only the number of bytes as specified in the header length. */
   if(ret = getc_bytes(&comment, cs, cbufptr, ebufptr)){
      free(comment);
      return(ret);
   }

   /* If comment did not explicitly contain a NULL terminator, it will */
   /* have one here by default due to the calloc of one extra byte at  */
   /* the end. */

   if(debug > 0)
      fprintf(stdf, "Comment =  %s", comment);

   *ocomment = comment;
   return(0);
}

/************************************************************************/
/* Build WSQ decomposition trees.                                       */
/************************************************************************/
void Wsq::build_wsq_trees(W_TREE w_tree[], const int w_treelen,
                     Q_TREE q_tree[], const int q_treelen,
                     const int width, const int height)
{
   /* Build a W-TREE structure for the image. */
   build_w_tree(w_tree, width, height);
   /* Build a Q-TREE structure for the image. */
   build_q_tree(w_tree, q_tree);
}

/********************************************************************/
/* Routine to obtain subband "x-y locations" for creating wavelets. */
/********************************************************************/
void Wsq::build_w_tree(
   W_TREE w_tree[],   /* wavelet tree structure */
   const int width,   /* image width            */
   const int height)  /* image height           */
{
   int lenx, lenx2, leny, leny2;  /* starting lengths of sections of
                                     the image being split into subbands */
   int node;

   for(node = 0; node < 20; node++) {
      w_tree[node].inv_rw = 0;
      w_tree[node].inv_cl = 0;
   }
   w_tree[2].inv_rw = 1;
   w_tree[4].inv_rw = 1;
   w_tree[7].inv_rw = 1;
   w_tree[9].inv_rw = 1;
   w_tree[11].inv_rw = 1;
   w_tree[13].inv_rw = 1;
   w_tree[16].inv_rw = 1;
   w_tree[18].inv_rw = 1;
   w_tree[3].inv_cl = 1;
   w_tree[5].inv_cl = 1;
   w_tree[8].inv_cl = 1;
   w_tree[9].inv_cl = 1;
   w_tree[12].inv_cl = 1;
   w_tree[13].inv_cl = 1;
   w_tree[17].inv_cl = 1;
   w_tree[18].inv_cl = 1;

   w_tree4(w_tree, 0, 1, width, height, 0, 0, 1);

   if((w_tree[1].lenx % 2) == 0) {
      lenx = w_tree[1].lenx / 2;
      lenx2 = lenx;
   }
   else {
      lenx = (w_tree[1].lenx + 1) / 2;
      lenx2 = lenx - 1;
   }

   if((w_tree[1].leny % 2) == 0) {
      leny = w_tree[1].leny / 2;
      leny2 = leny;
   }
   else {
      leny = (w_tree[1].leny + 1) / 2;
      leny2 = leny - 1;
   }

   w_tree4(w_tree, 4, 6, lenx2, leny, lenx, 0, 0);
   w_tree4(w_tree, 5, 10, lenx, leny2, 0, leny, 0);
   w_tree4(w_tree, 14, 15, lenx, leny, 0, 0, 0);

   w_tree[19].x = 0;
   w_tree[19].y = 0;
   if((w_tree[15].lenx % 2) == 0)
      w_tree[19].lenx = w_tree[15].lenx / 2;
   else
      w_tree[19].lenx = (w_tree[15].lenx + 1) / 2;

   if((w_tree[15].leny % 2) == 0)
      w_tree[19].leny = w_tree[15].leny / 2;
   else
      w_tree[19].leny = (w_tree[15].leny + 1) / 2;

   if(debug > 1) {
      for(node = 0; node < 20; node++)
         fprintf(stdf,
         "t%d -> x = %d  y = %d : dx = %d  dy = %d : ir = %d  ic = %d\n",
         node, w_tree[node].x, w_tree[node].y,
         w_tree[node].lenx, w_tree[node].leny,
         w_tree[node].inv_rw, w_tree[node].inv_cl);
      fprintf(stdf, "\n\n");
   }

   return;
}

/***************************************************************/
/* Gives location and size of subband splits for build_w_tree. */
/***************************************************************/
void Wsq::w_tree4(
   W_TREE w_tree[],    /* wavelet tree structure                      */
   const int start1,   /* w_tree locations to start calculating       */
   const int start2,   /*    subband split locations and sizes        */
   const int lenx,     /* (temp) subband split location and sizes     */
   const int leny,
   const int x,
   const int y,
   const int stop1)    /* 0 normal operation, 1 used to avoid marking */
                       /*    size and location of subbands 60-63      */
{
   int evenx, eveny;   /* Check length of subband for even or odd */
   int p1, p2;         /* w_tree locations for storing subband sizes and
                          locations */
   
   p1 = start1;
   p2 = start2;

   evenx = lenx % 2;
   eveny = leny % 2;

   w_tree[p1].x = x;
   w_tree[p1].y = y;
   w_tree[p1].lenx = lenx;
   w_tree[p1].leny = leny;
   
   w_tree[p2].x = x;
   w_tree[p2+2].x = x;
   w_tree[p2].y = y;
   w_tree[p2+1].y = y;

   if(evenx == 0) {
      w_tree[p2].lenx = lenx / 2;
      w_tree[p2+1].lenx = w_tree[p2].lenx;
   }
   else {
      if(p1 == 4) {
         w_tree[p2].lenx = (lenx - 1) / 2;
         w_tree[p2+1].lenx = w_tree[p2].lenx + 1;
      }
      else {
         w_tree[p2].lenx = (lenx + 1) / 2;
         w_tree[p2+1].lenx = w_tree[p2].lenx - 1;
      }
   }
   w_tree[p2+1].x = w_tree[p2].lenx + x;
   if(stop1 == 0) {
      w_tree[p2+3].lenx = w_tree[p2+1].lenx;
      w_tree[p2+3].x = w_tree[p2+1].x;
   }
   w_tree[p2+2].lenx = w_tree[p2].lenx;


   if(eveny == 0) {
      w_tree[p2].leny = leny / 2;
      w_tree[p2+2].leny = w_tree[p2].leny;
   }
   else {
      if(p1 == 5) {
         w_tree[p2].leny = (leny - 1) / 2;
         w_tree[p2+2].leny = w_tree[p2].leny + 1;
      }
      else {
         w_tree[p2].leny = (leny + 1) / 2;
         w_tree[p2+2].leny = w_tree[p2].leny - 1;
      }
   }
   w_tree[p2+2].y = w_tree[p2].leny + y;
   if(stop1 == 0) {
      w_tree[p2+3].leny = w_tree[p2+2].leny;
      w_tree[p2+3].y = w_tree[p2+2].y;
   }
   w_tree[p2+1].leny = w_tree[p2].leny;
}

/****************************************************************/
void Wsq::build_q_tree(
   W_TREE *w_tree,  /* wavelet tree structure */
   Q_TREE *q_tree)   /* quantization tree structure */
{
   int node;

   q_tree16(q_tree,3,w_tree[14].lenx,w_tree[14].leny,
              w_tree[14].x,w_tree[14].y, 0, 0);
   q_tree16(q_tree,19,w_tree[4].lenx,w_tree[4].leny,
              w_tree[4].x,w_tree[4].y, 0, 1);
   q_tree16(q_tree,48,w_tree[0].lenx,w_tree[0].leny,
              w_tree[0].x,w_tree[0].y, 0, 0);
   q_tree16(q_tree,35,w_tree[5].lenx,w_tree[5].leny,
              w_tree[5].x,w_tree[5].y, 1, 0);
   q_tree4(q_tree,0,w_tree[19].lenx,w_tree[19].leny,
             w_tree[19].x,w_tree[19].y);

   if(debug > 1) {
      for(node = 0; node < 60; node++)
         fprintf(stdf, "t%d -> x = %d  y = %d : lx = %d  ly = %d\n",
         node, q_tree[node].x, q_tree[node].y,
         q_tree[node].lenx, q_tree[node].leny);
      fprintf(stdf, "\n\n");
   }
   return;
}

/*****************************************************************/
void Wsq::q_tree16(
   Q_TREE *q_tree,   /* quantization tree structure */
   const int start,  /* q_tree location of first subband        */
                     /*   in the subband group being calculated */
   const int lenx,   /* (temp) subband location and sizes */
   const int leny,
   const int x,
   const int y,
   const int rw,  /* NEW */   /* spectral invert 1st row/col splits */
   const int cl)  /* NEW */
{
   int tempx, temp2x;   /* temporary x values */
   int tempy, temp2y;   /* temporary y values */
   int evenx, eveny;    /* Check length of subband for even or odd */
   int p;               /* indicates subband information being stored */

   p = start;
   evenx = lenx % 2;
   eveny = leny % 2;

   if(evenx == 0) {
      tempx = lenx / 2;
      temp2x = tempx;
   }
   else {
      if(cl) {
         temp2x = (lenx + 1) / 2;
         tempx = temp2x - 1;
      }
      else  {
        tempx = (lenx + 1) / 2;
        temp2x = tempx - 1;
      }
   }

   if(eveny == 0) {
      tempy = leny / 2;
      temp2y = tempy;
   }
   else {
      if(rw) {
         temp2y = (leny + 1) / 2;
         tempy = temp2y - 1;
      }
      else {
        tempy = (leny + 1) / 2;
        temp2y = tempy - 1;
      }
   }

   evenx = tempx % 2;
   eveny = tempy % 2;

   q_tree[p].x = x;
   q_tree[p+2].x = x;
   q_tree[p].y = y;
   q_tree[p+1].y = y;
   if(evenx == 0) {
      q_tree[p].lenx = tempx / 2;
      q_tree[p+1].lenx = q_tree[p].lenx;
      q_tree[p+2].lenx = q_tree[p].lenx;
      q_tree[p+3].lenx = q_tree[p].lenx;
   }
   else {
      q_tree[p].lenx = (tempx + 1) / 2;
      q_tree[p+1].lenx = q_tree[p].lenx - 1;
      q_tree[p+2].lenx = q_tree[p].lenx;
      q_tree[p+3].lenx = q_tree[p+1].lenx;
   }
   q_tree[p+1].x = x + q_tree[p].lenx;
   q_tree[p+3].x = q_tree[p+1].x;
   if(eveny == 0) {
      q_tree[p].leny = tempy / 2;
      q_tree[p+1].leny = q_tree[p].leny;
      q_tree[p+2].leny = q_tree[p].leny;
      q_tree[p+3].leny = q_tree[p].leny;
   }
   else {
      q_tree[p].leny = (tempy + 1) / 2;
      q_tree[p+1].leny = q_tree[p].leny;
      q_tree[p+2].leny = q_tree[p].leny - 1;
      q_tree[p+3].leny = q_tree[p+2].leny;
   }
   q_tree[p+2].y = y + q_tree[p].leny;
   q_tree[p+3].y = q_tree[p+2].y;


   evenx = temp2x % 2;

   q_tree[p+4].x = x + tempx;
   q_tree[p+6].x = q_tree[p+4].x;
   q_tree[p+4].y = y;
   q_tree[p+5].y = y;
   q_tree[p+6].y = q_tree[p+2].y;
   q_tree[p+7].y = q_tree[p+2].y;
   q_tree[p+4].leny = q_tree[p].leny;
   q_tree[p+5].leny = q_tree[p].leny;
   q_tree[p+6].leny = q_tree[p+2].leny;
   q_tree[p+7].leny = q_tree[p+2].leny;
   if(evenx == 0) {
      q_tree[p+4].lenx = temp2x / 2;
      q_tree[p+5].lenx = q_tree[p+4].lenx;
      q_tree[p+6].lenx = q_tree[p+4].lenx;
      q_tree[p+7].lenx = q_tree[p+4].lenx;
   }
   else {
      q_tree[p+5].lenx = (temp2x + 1) / 2;
      q_tree[p+4].lenx = q_tree[p+5].lenx - 1;
      q_tree[p+6].lenx = q_tree[p+4].lenx;
      q_tree[p+7].lenx = q_tree[p+5].lenx;
   }
   q_tree[p+5].x = q_tree[p+4].x + q_tree[p+4].lenx;
   q_tree[p+7].x = q_tree[p+5].x;


   eveny = temp2y % 2;

   q_tree[p+8].x = x;
   q_tree[p+9].x = q_tree[p+1].x;
   q_tree[p+10].x = x;
   q_tree[p+11].x = q_tree[p+1].x;
   q_tree[p+8].y = y + tempy;
   q_tree[p+9].y = q_tree[p+8].y;
   q_tree[p+8].lenx = q_tree[p].lenx;
   q_tree[p+9].lenx = q_tree[p+1].lenx;
   q_tree[p+10].lenx = q_tree[p].lenx;
   q_tree[p+11].lenx = q_tree[p+1].lenx;
   if(eveny == 0) {
      q_tree[p+8].leny = temp2y / 2;
      q_tree[p+9].leny = q_tree[p+8].leny;
      q_tree[p+10].leny = q_tree[p+8].leny;
      q_tree[p+11].leny = q_tree[p+8].leny;
   }
   else {
      q_tree[p+10].leny = (temp2y + 1) / 2;
      q_tree[p+11].leny = q_tree[p+10].leny;
      q_tree[p+8].leny = q_tree[p+10].leny - 1;
      q_tree[p+9].leny = q_tree[p+8].leny;
   }
   q_tree[p+10].y = q_tree[p+8].y + q_tree[p+8].leny;
   q_tree[p+11].y = q_tree[p+10].y;


   q_tree[p+12].x = q_tree[p+4].x;
   q_tree[p+13].x = q_tree[p+5].x;
   q_tree[p+14].x = q_tree[p+4].x;
   q_tree[p+15].x = q_tree[p+5].x;
   q_tree[p+12].y = q_tree[p+8].y;
   q_tree[p+13].y = q_tree[p+8].y;
   q_tree[p+14].y = q_tree[p+10].y;
   q_tree[p+15].y = q_tree[p+10].y;
   q_tree[p+12].lenx = q_tree[p+4].lenx;
   q_tree[p+13].lenx = q_tree[p+5].lenx;
   q_tree[p+14].lenx = q_tree[p+4].lenx;
   q_tree[p+15].lenx = q_tree[p+5].lenx;
   q_tree[p+12].leny = q_tree[p+8].leny;
   q_tree[p+13].leny = q_tree[p+8].leny;
   q_tree[p+14].leny = q_tree[p+10].leny;
   q_tree[p+15].leny = q_tree[p+10].leny;
}

/********************************************************************/
void Wsq::q_tree4(
   Q_TREE *q_tree,   /* quantization tree structure */
   const int start,  /* q_tree location of first subband         */
                     /*    in the subband group being calculated */
   const int lenx,   /* (temp) subband location and sizes */
   const int leny,
   const int x,
   const int  y)        
{
   int evenx, eveny;    /* Check length of subband for even or odd */
   int p;               /* indicates subband information being stored */


   p = start;
   evenx = lenx % 2;
   eveny = leny % 2;


   q_tree[p].x = x;
   q_tree[p+2].x = x;
   q_tree[p].y = y;
   q_tree[p+1].y = y;
   if(evenx == 0) {
      q_tree[p].lenx = lenx / 2;
      q_tree[p+1].lenx = q_tree[p].lenx;
      q_tree[p+2].lenx = q_tree[p].lenx;
      q_tree[p+3].lenx = q_tree[p].lenx;
   }
   else {
      q_tree[p].lenx = (lenx + 1) / 2;
      q_tree[p+1].lenx = q_tree[p].lenx - 1;
      q_tree[p+2].lenx = q_tree[p].lenx;
      q_tree[p+3].lenx = q_tree[p+1].lenx;
   }
   q_tree[p+1].x = x + q_tree[p].lenx;
   q_tree[p+3].x = q_tree[p+1].x;
   if(eveny == 0) {
      q_tree[p].leny = leny / 2;
      q_tree[p+1].leny = q_tree[p].leny;
      q_tree[p+2].leny = q_tree[p].leny;
      q_tree[p+3].leny = q_tree[p].leny;
   }
   else {
      q_tree[p].leny = (leny + 1) / 2;
      q_tree[p+1].leny = q_tree[p].leny;
      q_tree[p+2].leny = q_tree[p].leny - 1;
      q_tree[p+3].leny = q_tree[p+2].leny;
   }
   q_tree[p+2].y = y + q_tree[p].leny;
   q_tree[p+3].y = q_tree[p+2].y;
}

/******************************************************************/
/* This routine converts the unsigned char data to float.  In the */
/* process it shifts and scales the data so the values range from */
/* +/- 128.0                                                      */
/******************************************************************/
void Wsq::conv_img_2_flt(
   float *fip,         /* output float image data  */
   float *m_shift,     /* shifting parameter       */
   float *r_scale,     /* scaling parameter        */
   unsigned char *data,        /* input unsigned char data */
   const int num_pix)  /* num pixels in image      */

{
   int cnt;                     /* pixel cnt */
   unsigned int sum;                     /* sum of pixel values */
   float mean;                  /* mean pixel value */
   int low, high;               /* low/high pixel values */
   float low_diff, high_diff;   /* new low/high pixels values shifting */

   sum = 0;
   low = 255;
   high = 0;
   for(cnt = 0; cnt < num_pix; cnt++) 
	{
      if(data[cnt] > high)    high = data[cnt];
      if(data[cnt] < low)     low = data[cnt];
      sum += data[cnt];
   }

   mean = (float) sum / (float)num_pix;
   *m_shift = mean;

   low_diff = *m_shift - low;
   high_diff = high - *m_shift;

   if(low_diff >= high_diff)  *r_scale = low_diff; else  *r_scale = high_diff;

   *r_scale /= (float)128.0;
	float scale = (float)1./(*r_scale);
#ifdef SSE
	__m128 shift, scl;
	shift = _mm_load1_ps(m_shift);
	scl   = _mm_load1_ps(&scale);
	float *p = (float *)_aligned_malloc(4*sizeof(float),16);
   for(cnt = 0; cnt < num_pix; cnt+=4) 
	{
		p[0] = (float)data[cnt];
		p[1] = (float)data[cnt+1];
		p[2] = (float)data[cnt+2];
		p[3] = (float)data[cnt+3];
		__m128 d = _mm_load_ps(p);
		__m128 dif = _mm_sub_ps(d,shift);
		d = _mm_mul_ps(dif,scl);
//		_mm_store_ps(fip+cnt,d);
		_mm_stream_ps(fip+cnt,d);
   }
	_aligned_free(p);
#else
   for(cnt = 0; cnt < num_pix; cnt++) {
      fip[cnt] = ((float)data[cnt] - *m_shift) * scale;
   }
#endif
   return;
}

/*********************************************************/
/* Routine to convert image from float to unsigned char. */
/*********************************************************/
void Wsq::conv_img_2_uchar(
   unsigned char *data,                   /* uchar image pointer    */
   float *img,                    /* image pointer          */
   const int width,               /* image width            */
   const int height,              /* image height           */
   const float m_shift,           /* shifting parameter     */
   const float r_scale)           /* scaling parameter      */
{
   int r, c;       /* row/column counters */
   float img_tmp;  /* temp image data store */

   for (r = 0; r < height; r++) {
      for (c = 0; c < width; c++) {
         img_tmp = (*img * r_scale) + m_shift;
         img_tmp += 0.5;
         if (img_tmp < 0.0)
            *data = 0; /* neg pix poss after quantization */
         else if (img_tmp > 255.0)
            *data = 255;
         else
            *data = (unsigned char)img_tmp;

         ++img;
         ++data;
      }
   }
}

/**********************************************************/
/* This routine calculates the variances of the subbands. */
/**********************************************************/
void Wsq::variance(
   QUANT_VALS *quant_vals, /* quantization parameters */
   Q_TREE q_tree[],        /* quantization "tree"     */
   const int q_treelen,    /* length of q_tree        */
   float *fip,             /* image pointer           */
   const int width,        /* image width             */
   const int height)       /* image height            */
{
   float *fp;              /* temp image pointer */
   int cvr;                /* subband counter */
   int lenx = 0, leny = 0; /* dimensions of area to calculate variance */
   int skipx, skipy;       /* pixels to skip to get to area for
                              variance calculation */
   int row, col;           /* dimension counters */
   float ssq;             /* sum of squares */
   float sum2;            /* variance calculation parameter */
   float sum_pix;         /* sum of pixels */
   float vsum;            /* variance sum for subbands 0-3 */


   vsum = 0.0;
   for(cvr = 0; cvr < 4; cvr++) {
      fp = fip + (q_tree[cvr].y * width) + q_tree[cvr].x;
      ssq = 0.0;
      sum_pix = 0.0;

      skipx = q_tree[cvr].lenx / 8;
      skipy = (9 * q_tree[cvr].leny)/32;

      lenx = (3 * q_tree[cvr].lenx)/4;
      leny = (7 * q_tree[cvr].leny)/16;

      fp += (skipy * width) + skipx;
      for(row = 0; row < leny; row++, fp += (width - lenx)) {
         for(col = 0; col < lenx; col++) {
            sum_pix += *fp;
            ssq += *fp * *fp;
            fp++;
         }
      }
      sum2 = (sum_pix * sum_pix)/(lenx * leny);
      quant_vals->var[cvr] = (float)((ssq - sum2)/((lenx * leny)-1.0));
      vsum += quant_vals->var[cvr];
   }

   if(vsum < 20000.0) {
      for(cvr = 0; cvr < NUM_SUBBANDS; cvr++) {
         fp = fip + (q_tree[cvr].y * width) + q_tree[cvr].x;
         ssq = 0.0;
         sum_pix = 0.0;

         lenx = q_tree[cvr].lenx;
         leny = q_tree[cvr].leny;

         for(row = 0; row < leny; row++, fp += (width - lenx)) {
            for(col = 0; col < lenx; col++) {
               sum_pix += *fp;
               ssq += *fp * *fp;
               fp++;
            }
         }
         sum2 = (sum_pix * sum_pix)/(lenx * leny);
         quant_vals->var[cvr] = (float)((ssq - sum2)/((lenx * leny)-1.0));
      }
   }
   else {
      for(cvr = 4; cvr < NUM_SUBBANDS; cvr++) {
         fp = fip + (q_tree[cvr].y * width) + q_tree[cvr].x;
         ssq = 0.0;
         sum_pix = 0.0;

         skipx = q_tree[cvr].lenx / 8;
         skipy = (9 * q_tree[cvr].leny)/32;

         lenx = (3 * q_tree[cvr].lenx)/4;
         leny = (7 * q_tree[cvr].leny)/16;

         fp += (skipy * width) + skipx;
         for(row = 0; row < leny; row++, fp += (width - lenx)) {
            for(col = 0; col < lenx; col++) {
               sum_pix += *fp;
               ssq += *fp * *fp;
               fp++;
            }
         }
         sum2 = (sum_pix * sum_pix)/(lenx * leny);
         quant_vals->var[cvr] = (float)((ssq - sum2)/((lenx * leny)-1.0));
      }
   }
}

/************************************************/
/* This routine quantizes the wavelet subbands. */
/************************************************/
int Wsq::quantize(
   short **osip,           /* quantized output             */
   int *ocmp_siz,          /* size of quantized output     */
   QUANT_VALS *quant_vals, /* quantization parameters      */
   Q_TREE q_tree[],        /* quantization "tree"          */
   const int q_treelen,    /* size of q_tree               */
   float *fip,             /* floating point image pointer */
   const int width,        /* image width                  */
   const int height)       /* image height                 */
{
   int i;                 /* temp counter */
   int j;                 /* interation index */
   float *fptr;           /* temp image pointer */
   short *sip, *sptr;     /* pointers to quantized image */
   int row, col;          /* temp image characteristic parameters */
   int cnt;               /* subband counter */
   double zbin;            /* zero bin size */
   double A[NUM_SUBBANDS]; /* subband "weights" for quantization */
   double m[NUM_SUBBANDS]; /* subband size to image size ratios */
                          /* (reciprocal of FBI spec for 'm')  */
   double m1, m2, m3;      /* reciprocal constants for 'm' */
   double sigma[NUM_SUBBANDS]; /* square root of subband variances */
   //float zbin;            /* zero bin size */
   //float A[NUM_SUBBANDS]; /* subband "weights" for quantization */
   //float m[NUM_SUBBANDS]; /* subband size to image size ratios */
   //                       /* (reciprocal of FBI spec for 'm')  */
   //float m1, m2, m3;      /* reciprocal constants for 'm' */
   //float sigma[NUM_SUBBANDS]; /* square root of subband variances */
   int K0[NUM_SUBBANDS];  /* initial list of subbands w/variance >= thresh */
   int K1[NUM_SUBBANDS];  /* working list of subbands */
   int *K, *nK;           /* pointers to sets of subbands */
   int NP[NUM_SUBBANDS];  /* current subbounds with nonpositive bit rates. */
   int K0len;             /* number of subbands in K0 */
   int Klen, nKlen;       /* number of subbands in other subband lists */
   int NPlen;             /* number of subbands flagged in NP */
   double S;               /* current frac of subbands w/positive bit rate */
   double q;               /* current proportionality constant */
   double P;               /* product of 'q/Q' ratios */
   //float S;               /* current frac of subbands w/positive bit rate */
   //float q;               /* current proportionality constant */
   //float P;               /* product of 'q/Q' ratios */

   /* Set up 'A' table. */
   for(cnt = 0; cnt < STRT_SUBBAND_3; cnt++)
      A[cnt] = 1.0;
   A[cnt++ /*52*/] = 1.32;
   A[cnt++ /*53*/] = 1.08;
   A[cnt++ /*54*/] = 1.42;
   A[cnt++ /*55*/] = 1.08;
   A[cnt++ /*56*/] = 1.32;
   A[cnt++ /*57*/] = 1.42;
   A[cnt++ /*58*/] = 1.08;
   A[cnt++ /*59*/] = 1.08;

   for(cnt = 0; cnt < MAX_SUBBANDS; cnt++) {
      quant_vals->qbss[cnt] = 0.0;
      quant_vals->qzbs[cnt] = 0.0;
   }

   /* Set up 'Q1' (prime) table. */
   for(cnt = 0; cnt < NUM_SUBBANDS; cnt++) {
      if(quant_vals->var[cnt] < VARIANCE_THRESH)
         quant_vals->qbss[cnt] = 0.0;
      else
         /* NOTE: q has been taken out of the denominator in the next */
         /*       2 formulas from the original code. */
         if(cnt < STRT_SIZE_REGION_2 /*4*/)
            quant_vals->qbss[cnt] = 1.0;
         else
            quant_vals->qbss[cnt] = 10.0 / (A[cnt] *
                                    (double)log(quant_vals->var[cnt]));
   }


   /* Set up output buffer. */
   if((sip = (short *) calloc(width*height, sizeof(short))) == NULL) {
      fprintf(stderr,"ERROR : quantize : calloc : sip\n");
      return(-90);
   }
   sptr = sip;

   /* Set up 'm' table (these values are the reciprocal of 'm' in */
   /* the FBI spec).                                              */
   m1 = 1.0/1024.0;
   m2 = 1.0/256.0;
   m3 = 1.0/16.0;
   for(cnt = 0; cnt < STRT_SIZE_REGION_2; cnt++)
      m[cnt] = m1;
   for(cnt = STRT_SIZE_REGION_2; cnt < STRT_SIZE_REGION_3; cnt++)
      m[cnt] = m2;
   for(cnt = STRT_SIZE_REGION_3; cnt < NUM_SUBBANDS; cnt++)
      m[cnt] = m3;

   j = 0;
   /* Initialize 'K0' and 'K1' lists. */
   K0len = 0;
   for(cnt = 0; cnt < NUM_SUBBANDS; cnt++){
      if(quant_vals->var[cnt] >= VARIANCE_THRESH){
         K0[K0len] = cnt;
         K1[K0len++] = cnt;
         /* Compute square root of subband variance. */
         sigma[cnt] = sqrt(quant_vals->var[cnt]);
      }
   }
   K = K1;
   Klen = K0len;

   while(1){
      /* Compute new 'S' */
      S = 0.0;
      for(i = 0; i < Klen; i++){
         /* Remeber 'm' is the reciprocal of spec. */
         S += m[K[i]];
      }

      /* Compute product 'P' */
      P = 1.0;
      for(i = 0; i < Klen; i++){
         /* Remeber 'm' is the reciprocal of spec. */
         P *= pow((sigma[K[i]] / quant_vals->qbss[K[i]]), m[K[i]]);
      }

      /* Compute new 'q' */
//      q = (pow(2,((quant_vals->r/S)-1.0))/2.5) / pow(P, (1.0/S));
//      q = (pow((float)2,(float)((quant_vals->r/S)-1.0))/2.5) / pow(P, (float)(1.0/S)); // 3 err
      q = (pow((double)2,((quant_vals->r/S)-1.0))/2.5) / pow((double)P, (1.0/S)); 

      /* Flag subbands with non-positive bitrate. */
      memset(NP, 0, NUM_SUBBANDS * sizeof(int));
      NPlen = 0;
      for(i = 0; i < Klen; i++){
         if((quant_vals->qbss[K[i]] / q) >= (5.0*sigma[K[i]])){
            NP[K[i]] = TRUE;
            NPlen++;
         }
      }

      /* If list of subbands with non-positive bitrate is empty ... */
      if(NPlen == 0){
         /* Then we are done, so break from while loop. */
         break;
      }

      /* Assign new subband set to previous set K minus subbands in set NP. */
      nK = K1;
      nKlen = 0;
      for(i = 0; i < Klen; i++){
         if(!NP[K[i]])
            nK[nKlen++] = K[i];
      }

      /* Assign new set as K. */
      K = nK;
      Klen = nKlen;

      /* Bump iteration counter. */
      j++;
   }

   /* Flag subbands that are in set 'K0' (the very first set). */
   nK = K1;
   memset(nK, 0, NUM_SUBBANDS * sizeof(int));
   for(i = 0; i < K0len; i++){
      nK[K0[i]] = TRUE;
   }
   /* Set 'Q' values. */
   for(cnt = 0; cnt < NUM_SUBBANDS; cnt++) {
      if(nK[cnt])
	  {
         quant_vals->qbss[cnt] /= q;
//         quant_vals->qbss[cnt] += 0.0000005;
	  }
      else
         quant_vals->qbss[cnt] = 0.0;
      quant_vals->qzbs[cnt] = 1.2 * quant_vals->qbss[cnt];
   }

   /* Now ready to compute and store bin widths for subbands. */
   for(cnt = 0; cnt < NUM_SUBBANDS; cnt++) {
      fptr = fip + (q_tree[cnt].y * width) + q_tree[cnt].x;

      if(quant_vals->qbss[cnt] != 0.0) {

         zbin = quant_vals->qzbs[cnt] / 2.0;

         for(row = 0;
            row < q_tree[cnt].leny;
            row++, fptr += width - q_tree[cnt].lenx){
            for(col = 0; col < q_tree[cnt].lenx; col++) {
               if(-zbin <= *fptr && *fptr <= zbin)
                  *sptr = 0;
               else if(*fptr > 0.0)
                  *sptr = (short)(((*fptr-zbin)/quant_vals->qbss[cnt]) + 1.0);
               else
                  *sptr = (short)(((*fptr+zbin)/quant_vals->qbss[cnt]) - 1.0);
               sptr++;
               fptr++;
            }
         }
      }
      else if(debug > 0)
         fprintf(stderr, "%d -> %3.6f\n", cnt, quant_vals->qbss[cnt]);
   }

   *osip = sip;
   *ocmp_siz = sptr - sip;
   return(0);
}

/************************************************************************/
/* Compute quantized WSQ subband block sizes.                           */
/************************************************************************/
void Wsq::quant_block_sizes(int *oqsize1, int *oqsize2, int *oqsize3,
                 QUANT_VALS *quant_vals,
                 W_TREE w_tree[], const int w_treelen,
                 Q_TREE q_tree[], const int q_treelen)
{
   int qsize1, qsize2, qsize3;
   int node;

   /* Compute temporary sizes of 3 WSQ subband blocks. */
   qsize1 = w_tree[14].lenx * w_tree[14].leny;
   qsize2 = (w_tree[5].leny * w_tree[1].lenx) +
            (w_tree[4].lenx * w_tree[4].leny);
   qsize3 = (w_tree[2].lenx * w_tree[2].leny) +
            (w_tree[3].lenx * w_tree[3].leny);

   /* Adjust size of quantized WSQ subband blocks. */
   for (node = 0; node < STRT_SUBBAND_2; node++)
      if(quant_vals->qbss[node] == 0.0)
         qsize1 -= (q_tree[node].lenx * q_tree[node].leny);

   for (node = STRT_SUBBAND_2; node < STRT_SUBBAND_3; node++)
      if(quant_vals->qbss[node] == 0.0)
          qsize2 -= (q_tree[node].lenx * q_tree[node].leny);

   for (node = STRT_SUBBAND_3; node < STRT_SUBBAND_DEL; node++)
      if(quant_vals->qbss[node] == 0.0)
         qsize3 -= (q_tree[node].lenx * q_tree[node].leny);

   *oqsize1 = qsize1;
   *oqsize2 = qsize2;
   *oqsize3 = qsize3;
}

/*************************************/
/* Routine to unquantize image data. */
/*************************************/
int Wsq::unquantize(
   float **ofip,         /* floating point image pointer         */
   const DQT_TABLE *dqt_table, /* quantization table structure   */
   Q_TREE q_tree[],      /* quantization table structure         */
   const int q_treelen,  /* size of q_tree                       */
   short *sip,           /* quantized image pointer              */
   const int width,      /* image width                          */
   const int height)     /* image height                         */
{
   float *fip;    /* floating point image */
   int row, col;  /* cover counter and row/column counters */
   float C;       /* quantizer bin center */
   float *fptr;   /* image pointers */
   short *sptr;
   int cnt;       /* subband counter */

   if((fip = (float *) calloc(width*height, sizeof(float))) == NULL) {
      fprintf(stdf,"ERROR : unquantize : calloc : fip\n");
      return(-91);
   }
   if(dqt_table->dqt_def != 1) {
      fprintf(stdf,
      "ERROR: unquantize : quantization table parameters not defined!\n");
      return(-92);
   }

   sptr = sip;
   C = dqt_table->bin_center;
   for(cnt = 0; cnt < NUM_SUBBANDS; cnt++) 
   {
      if(dqt_table->q_bin[cnt] == 0.0)
         continue;
      fptr = fip + (q_tree[cnt].y * width) + q_tree[cnt].x;

      for(row = 0; row < q_tree[cnt].leny; row++, fptr += width - q_tree[cnt].lenx)
	  {
         for(col = 0; col < q_tree[cnt].lenx; col++) 
		 {
            if(*sptr == 0)
               *fptr = 0.0;
            else if(*sptr > 0)
               *fptr = (dqt_table->q_bin[cnt] * ((double)*sptr - C))
                    + (dqt_table->z_bin[cnt] / (double)2.0);
//               *fptr = (dqt_table->q_bin[cnt] * ((float)*sptr - C))
  //                  + (dqt_table->z_bin[cnt] / (float)2.0);
            else if(*sptr < 0)
               *fptr = (dqt_table->q_bin[cnt] * ((float)*sptr + C))
                    - (dqt_table->z_bin[cnt] / (float)2.0);
               //*fptr = (dqt_table->q_bin[cnt] * ((float)*sptr + C))
               //     - (dqt_table->z_bin[cnt] / (float)2.0);
            else {
               fprintf(stdf,
               "ERROR : unquantize : invalid quantization pixel value\n");
               return(-93);
            }
            fptr++;
            sptr++;
         }
      }
   }

   *ofip = fip;
   return(0);
}

/************************************************************************/
/* WSQ decompose the image.  NOTE: this routine modifies and returns    */
/* the results in "fdata".                                              */
/************************************************************************/
int Wsq::wsq_decompose(float *fdata, const int width, const int height,
                  W_TREE w_tree[], const int w_treelen,
                  double *hifilt, const int hisz, double *lofilt, const int losz)
//                  float *hifilt, const int hisz, float *lofilt, const int losz)
{
   int num_pix, node;
   float *fdata1, *fdata_bse;

   num_pix = width * height;
   /* Allocate temporary floating point pixmap. */
   if((fdata1 = (float *) malloc(num_pix*sizeof(float))) == NULL) {
      fprintf(stdf,"ERROR : wsq_decompose : malloc : fdata1\n");
      return(-94);
   }

   /* Compute the Wavelet image decomposition. */
   for(node = 0; node < w_treelen; node++) {
      fdata_bse = fdata + (w_tree[node].y * width) + w_tree[node].x;
      get_lets(fdata1, fdata_bse, w_tree[node].leny, w_tree[node].lenx, width, 1, hifilt, hisz, lofilt, losz, w_tree[node].inv_rw);
      get_lets(fdata_bse, fdata1, w_tree[node].lenx, w_tree[node].leny, 1, width, hifilt, hisz, lofilt, losz, w_tree[node].inv_cl);
   }
   free(fdata1);

   return(0);
}

/************************************************************/
/************************************************************/
void Wsq::get_lets(
   float *new1,     /* image pointers for creating subband splits */
   float *old,
   const int len1,       /* temporary length parameters */
   const int len2,
   const int pitch,      /* pitch gives next row_col to filter */
   const int  stride,    /*           stride gives next pixel to filter */
   double *hi,
//   float *hi,
   const int hsz,   /* NEW */
   double *lo,      /* filter coefficients */
//   float *lo,      /* filter coefficients */
   const int lsz,   /* NEW */
   const int inv)        /* spectral inversion? */
{
   float *lopass, *hipass;	/* pointers of where to put lopass
                                   and hipass filter outputs */
   float *p0,*p1;		/* pointers to image pixels used */
   int pix, rw_cl;		/* pixel counter and row/column counter */
   int i, da_ev;		/* even or odd row/column of pixels */
   int fi_ev;
   int loc, hoc, nstr, pstr;
//   int oloc, ohoc;
   int llen, hlen;
   int lpxstr, lspxstr;
   float *lpx, *lspx; 
   int hpxstr, hspxstr;
   float *hpx, *hspx;
   int olle, ohle;
   int olre, ohre;
   int lle, lle2;
   int lre, lre2;
   int hle, hle2;
   int hre, hre2;


   da_ev = len2 % 2;
   fi_ev = lsz % 2;

   if(fi_ev) {
	   loc = (lsz-1)/2;
	   hoc = (hsz-1)/2 - 1;
	   olle = 0;
	   ohle = 0;
	   olre = 0;
	   ohre = 0;
   }
   else {
	   loc = lsz/2 - 2;
	   hoc = hsz/2 - 2;
	   olle = 1;
	   ohle = 1;
	   olre = 1;
	   ohre = 1;

	   if(loc == -1) {
		   loc = 0;
		   olle = 0;
	   }
	   if(hoc == -1) {
		   hoc = 0;
		   ohle = 0;
	   }

	   for(i = 0; i < hsz; i++)
		   hi[i] *= -1.0;
   }

   pstr = stride;
   nstr = -pstr;

   if(da_ev) {
	   llen = (len2+1)/2;
	   hlen = llen - 1;
   }
   else {
	   llen = len2/2;
	   hlen = llen;
   }


   for(rw_cl = 0; rw_cl < len1; rw_cl++) {

	   if(inv) {
		   hipass = new1 + rw_cl * pitch;
		   lopass = hipass + hlen * stride;
	   }
	   else {
		   lopass = new1 + rw_cl * pitch;
		   hipass = lopass + llen * stride;
	   }

	   p0 = old + rw_cl * pitch;
	   p1 = p0 + (len2-1) * stride;

	   lspx = p0 + (loc * stride);
	   lspxstr = nstr;
	   lle2 = olle;
	   lre2 = olre;
	   hspx = p0 + (hoc * stride);
	   hspxstr = nstr;
	   hle2 = ohle;
	   hre2 = ohre;
	   for(pix = 0; pix < hlen; pix++) {
		   lpxstr = lspxstr;
		   lpx = lspx;
		   lle = lle2;
		   lre = lre2;
		   *lopass = *lpx * lo[0];
		   for(i = 1; i < lsz; i++) {
			   if(lpx == p0) {
				   if(lle) {
					   lpxstr = 0;
					   lle = 0;
				   }
				   else
					   lpxstr = pstr;
			   }
			   if(lpx == p1) {
				   if(lre) {
					   lpxstr = 0;
					   lre = 0;
				   }
				   else
					   lpxstr = nstr;
			   }
			   lpx += lpxstr;
			   *lopass += *lpx * lo[i];
		   }
		   lopass += stride;

		   hpxstr = hspxstr;
		   hpx = hspx;
		   hle = hle2;
		   hre = hre2;
		   *hipass = *hpx * hi[0];
		   for(i = 1; i < hsz; i++) {
			   if(hpx == p0) {
				   if(hle) {
					   hpxstr = 0;
					   hle = 0;
				   }
				   else
					   hpxstr = pstr;
			   }
			   if(hpx == p1) {
				   if(hre) {
					   hpxstr = 0;
					   hre = 0;
				   }
				   else
					   hpxstr = nstr;
			   }
			   hpx += hpxstr;
			   *hipass += *hpx * hi[i];
		   }
		   hipass += stride;

		   for(i = 0; i < 2; i++) {
			   if(lspx == p0) {
				   if(lle2) {
					   lspxstr = 0;
					   lle2 = 0;
				   }
				   else
					   lspxstr = pstr;
			   }
			   lspx += lspxstr;
			   if(hspx == p0) {
				   if(hle2) {
					   hspxstr = 0;
					   hle2 = 0;
				   }
				   else
					   hspxstr = pstr;
			   }
			   hspx += hspxstr;
		   }
	   }
	   if(da_ev) {
		   lpxstr = lspxstr;
		   lpx = lspx;
		   lle = lle2;
		   lre = lre2;
		   *lopass = *lpx * lo[0];
		   for(i = 1; i < lsz; i++) {
			   if(lpx == p0) {
				   if(lle) {
					   lpxstr = 0;
					   lle = 0;
				   }
				   else
					   lpxstr = pstr;
			   }
			   if(lpx == p1) {
				   if(lre) {
					   lpxstr = 0;
					   lre = 0;
				   }
				   else
					   lpxstr = nstr;
			   }
			   lpx += lpxstr;
			   *lopass += *lpx * lo[i];
		   }
		   lopass += stride;
	   }
   }
   if(!fi_ev) {
	   for(i = 0; i < hsz; i++)
		   hi[i] *= -1.0;
   }
}

/************************************************************************/
/* WSQ reconstructs the image.  NOTE: this routine modifies and returns */
/* the results in "fdata".                                              */
/************************************************************************/
int Wsq::wsq_reconstruct(float *fdata, const int width, const int height,
                  W_TREE w_tree[], const int w_treelen,
                  const DTT_TABLE *dtt_table)
{
   int num_pix, node;
   float *fdata1, *fdata_bse;

   if(dtt_table->lodef != 1) {
      fprintf(stdf,
      "ERROR: wsq_reconstruct : Lopass filter coefficients not defined\n");
      return(-95);
   }
   if(dtt_table->hidef != 1) {
      fprintf(stdf,
      "ERROR: wsq_reconstruct : Hipass filter coefficients not defined\n");
      return(-96);
   }

   num_pix = width * height;
   /* Allocate temporary floating point pixmap. */
   if((fdata1 = (float *) malloc(num_pix*sizeof(float))) == NULL) {
      fprintf(stdf,"ERROR : wsq_reconstruct : malloc : fdata1\n");
      return(-97);
   }

   /* Reconstruct floating point pixmap from wavelet subband data. */
   for (node = w_treelen - 1; node >= 0; node--) 
	{
      fdata_bse = fdata + (w_tree[node].y * width) + w_tree[node].x;
      join_lets(fdata1, fdata_bse, w_tree[node].lenx, w_tree[node].leny,
      1, width,dtt_table->hifilt, dtt_table->hisz,dtt_table->lofilt, dtt_table->losz,
                  w_tree[node].inv_cl);
      join_lets(fdata_bse, fdata1, w_tree[node].leny, w_tree[node].lenx,
      width, 1,dtt_table->hifilt, dtt_table->hisz,dtt_table->lofilt, dtt_table->losz,
                  w_tree[node].inv_rw);
   }
   free(fdata1);

   return(0);
}

/****************************************************************/
void  Wsq::join_lets(
   float *new1,    /* image pointers for creating subband splits */
   float *old,
   const int len1,       /* temporary length parameters */
   const int len2,
   const int pitch,      /* pitch gives next row_col to filter */
   const int  stride,    /*           stride gives next pixel to filter */
   double *hi,
//   float *hi,
   const int hsz,   /* NEW */
   double *lo,      /* filter coefficients */
//   float *lo,      /* filter coefficients */
   const int lsz,   /* NEW */
   const int inv)        /* spectral inversion? */
{
   float *lp0, *lp1;
   float *hp0, *hp1;
   float *lopass, *hipass;	/* lo/hi pass image pointers */
   float *limg, *himg;
   int pix, cl_rw;		/* pixel counter and column/row counter */
   int i, da_ev;			/* if "scanline" is even or odd and */
   int loc, hoc;
   int hlen, llen;
   int nstr, pstr;
   int tap;
   int fi_ev;
   int olle, ohle, olre, ohre;
   int lle, lle2, lre, lre2;
   int hle, hle2, hre, hre2;
   float *lpx, *lspx;
   int lpxstr, lspxstr;
   int lstap, lotap;
   float *hpx, *hspx;
   int hpxstr, hspxstr;
   int hstap, hotap;
   int asym, fhre = 0, ofhre;
   float ssfac, osfac, sfac;

   da_ev = len2 % 2;
   fi_ev = lsz % 2;
   pstr = stride;
   nstr = -pstr;
   if (da_ev) {
      llen = (len2+1)/2;
      hlen = llen - 1;
   }
   else {
      llen = len2/2;
      hlen = llen;
   }

   if (fi_ev) {
      asym = 0;
      ssfac = 1.0;
      ofhre = 0;
      loc = (lsz-1)/4;
      hoc = (hsz+1)/4 - 1;
      lotap = ((lsz-1)/2) % 2;
      hotap = ((hsz+1)/2) % 2;
      if(da_ev) {
         olle = 0;
         olre = 0;
         ohle = 1;
         ohre = 1;
      }
      else {
         olle = 0;
         olre = 1;
         ohle = 1;
         ohre = 0;
      }
   }
   else 
	{
      asym = 1;
      ssfac = -1.0;
      ofhre = 2;
      loc = lsz/4 - 1;
      hoc = hsz/4 - 1;
      lotap = (lsz/2) % 2;
      hotap = (hsz/2) % 2;
      if(da_ev) {
         olle = 1;
         olre = 0;
         ohle = 1;
         ohre = 1;
      }
      else {
         olle = 1;
         olre = 1;
         ohle = 1;
         ohre = 1;
      }

      if(loc == -1) {
         loc = 0;
         olle = 0;
      }
      if(hoc == -1) {
         hoc = 0;
         ohle = 0;
      }

      for(i = 0; i < hsz; i++) hi[i] =- hi[i];
//         hi[i] *= -1.0;
   }

   for (cl_rw = 0; cl_rw < len1; cl_rw++)  
	{
      limg = new1 + cl_rw * pitch;
      himg = limg;
      *himg = 0.0;
      *(himg + stride) = 0.0;
      if(inv) {
         hipass = old + cl_rw * pitch;
         lopass = hipass + stride * hlen;
      }
      else {
         lopass = old + cl_rw * pitch;
         hipass = lopass + stride * llen;
      }


      lp0 = lopass;
      lp1 = lp0 + (llen-1) * stride;
      lspx = lp0 + (loc * stride);
      lspxstr = nstr;
      lstap = lotap;
      lle2 = olle;
      lre2 = olre;

      hp0 = hipass;
      hp1 = hp0 + (hlen-1) * stride;
      hspx = hp0 + (hoc * stride);
      hspxstr = nstr;
      hstap = hotap;
      hle2 = ohle;
      hre2 = ohre;
      osfac = ssfac;

      for(pix = 0; pix < hlen; pix++) 
		{ 
         for(tap = lstap; tap >=0; tap--) 
			{
            lle = lle2;
            lre = lre2;
            lpx = lspx;
            lpxstr = lspxstr;

            *limg = *lpx * lo[tap];
            for(i = tap+2; i < lsz; i += 2) 
				{
					if(lpx == lp0) {
                  if(lle) {
                     lpxstr = 0;
                     lle = 0;
                  }
                  else
                     lpxstr = pstr;
					}
               if(lpx == lp1) {
                  if(lre) {
                     lpxstr = 0;
                     lre = 0;
                  }
                  else
                     lpxstr = nstr;
               }
               lpx += lpxstr;
               *limg += *lpx * lo[i];
            }
            limg += stride;
         }
		 if(lspx == lp0) {
            if(lle2) {
               lspxstr = 0;
               lle2 = 0;
            }
            else
               lspxstr = pstr;
			}
         lspx += lspxstr;
         lstap = 1;

         for(tap = hstap; tap >=0; tap--) 
			{
            hle = hle2;
            hre = hre2;
            hpx = hspx;
            hpxstr = hspxstr;
            fhre = ofhre;
            sfac = osfac;

            for(i = tap; i < hsz; i += 2) 
				{
               if(hpx == hp0) {
                  if(hle) {
                     hpxstr = 0;
                     hle = 0;
                  }
                  else {
                     hpxstr = pstr;
                     sfac = 1.0;
                  }
               }
               if(hpx == hp1) {
                  if(hre) {
                     hpxstr = 0;
                     hre = 0;
                     if(asym && da_ev) {
                        hre = 1;
                        fhre--;
                        sfac = (float)fhre;
                        if(sfac == 0.0)
                           hre = 0;
                     }
                  }
                  else {
                     hpxstr = nstr;
                     if(asym)
                        sfac = -1.0;
                  }
               }
               *himg += *hpx * hi[i] * sfac;
               hpx += hpxstr;
            }
            himg += stride;
         }
         if(hspx == hp0) {
            if(hle2) {
               hspxstr = 0;
               hle2 = 0;
            }
            else {
               hspxstr = pstr;
               osfac = 1.0;
            }
         }
         hspx += hspxstr;
         hstap = 1;
      }


      if(da_ev)
         if(lotap)
            lstap = 1;
         else
            lstap = 0;
      else
         if(lotap)
            lstap = 2;
         else
            lstap = 1;

      for(tap = 1; tap >= lstap; tap--) 
		{
         lle = lle2;
         lre = lre2;
         lpx = lspx;
         lpxstr = lspxstr;

         *limg = *lpx * lo[tap];
         for(i = tap+2; i < lsz; i += 2) 
			{
				if(lpx == lp0) {
               if(lle) {
                  lpxstr = 0;
                  lle = 0;
               }
               else
                  lpxstr = pstr;
				}
            if(lpx == lp1) {
               if(lre) {
                  lpxstr = 0;
                  lre = 0;
               }
               else
                  lpxstr = nstr;
            }
            lpx += lpxstr;
            *limg += *lpx * lo[i];
         }
         limg += stride;
      }


      if(da_ev) {
         if(hotap)
            hstap = 1;
         else
            hstap = 0;

         if(hsz == 2) {
            hspx -= hspxstr;
            fhre = 1;
         }
      }
      else
         if(hotap)
            hstap = 2;
         else
            hstap = 1;


      for(tap = 1; tap >= hstap; tap--) 
		{
         hle = hle2;
         hre = hre2;
         hpx = hspx;
         hpxstr = hspxstr;
         sfac = osfac;
         if(hsz != 2)
            fhre = ofhre;

         for(i = tap; i < hsz; i += 2) 
			{
            if(hpx == hp0) {
               if(hle) {
                  hpxstr = 0;
                  hle = 0;
               }
               else {
                  hpxstr = pstr;
                  sfac = 1.0;
               }
            }
            if(hpx == hp1) {
               if(hre) {
                  hpxstr = 0;
                  hre = 0;
                  if(asym && da_ev) {
                     hre = 1;
                     fhre--;
                     sfac = (float)fhre;
                     if(sfac == 0.0)
                        hre = 0;
                  }
               }
               else {
                  hpxstr = nstr;
                  if(asym)
                     sfac = -1.0;
               }
            }
            *himg += *hpx * hi[i] * sfac;
            hpx += hpxstr;
         }
         himg += stride;
      }
   }			

   if(!fi_ev)
      for(i = 0; i < hsz; i++)hi[i] =- hi[i];
//         hi[i] *= -1.0;
}

/*****************************************************/
/* Routine to execute an integer  sign determination */
/*****************************************************/

int Wsq::int_sign(const int power)  /* "sign" power */
{
   int cnt, num = -1;   /* counter and sign return value */

   if(power == 0)
      return 1;

   for(cnt = 1; cnt < power; cnt++)
      num *= -1;

   return num;
}

/*************************************************************/
/* Added by MDG on 02-24-05                                  */
/* Initializes memory used by the WSQ decoder.               */
/*************************************************************/
void Wsq::init_wsq_decoder_resources(DTT_TABLE *dtt_table)
{
   /* Added 02-24-05 by MDG                      */
   /* Init dymanically allocated members to NULL */
   /* for proper memory management in:           */
   /*    read_transform_table()                  */
   /*    getc_transform_table()                  */
   /*    free_wsq_resources()                    */
   dtt_table->lofilt = (double *)NULL;
   dtt_table->hifilt = (double *)NULL;
}

/*************************************************************/
/* Added by MDG on 02-24-05                                  */
/* Deallocates memory used by the WSQ decoder.               */
/*************************************************************/
void Wsq::free_wsq_decoder_resources(DTT_TABLE *dtt_table)
{
   if(dtt_table->lofilt != (double *)NULL){
      free(dtt_table->lofilt);
      dtt_table->lofilt = (double *)NULL;
   }

   if(dtt_table->hifilt != (double *)NULL){
      free(dtt_table->hifilt);
      dtt_table->hifilt = (double *)NULL;
   }
}

//}