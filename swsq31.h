#ifndef WSQ31_H_
#define WSQ31_H_

#ifdef SWSQ_EXPORTS
   #define WSQ_EI __declspec(dllexport)
#else  
   #define WSQ_EI 
#endif

#define  ERROR_MEMORY               -1    // memory error
#define  ERROR_PROCESSING           -2    // processing error 
#define  ERROR_FILE_CONTENT         -3    // wrong content of WSQ file
#define  ERROR_FILTERS              -4    // error in WSQ filters
#define  EMPTY_SOURCE	            -5    // pointer to source buffer is NULL
#define  ERROR_SIZE                 -6    // one of source image size is less than 128
#define  ERROR_EXCEPTION            -7    // unknown exception raised
#define  ERROR_BITRATE				-8    // wrong bitrate parameter
#define  ERROR_CHECK_SUM          -666    // control sum error

typedef void *(*type_init)();
typedef void (*type_cancel)(void *pInstance);
typedef int  (*type_encode)(void *pInstance, unsigned char *idata, unsigned char *odata, int width, int height, float r_bitrate, int chksum);
typedef int  (*type_decode)(void *pInstance, unsigned char *odata, unsigned char *idata, int ilen, int *width, int *height, int chksum);

extern "C" {
void WSQ_EI *init_instance();
void WSQ_EI cancel_instance(void *pInstance);
int  WSQ_EI wsq_encode(void *pInstance, unsigned char *idata, unsigned char *odata, int width, int height, float r_bitrate, int chksum);
int  WSQ_EI wsq_decode(void *pInstance, unsigned char *odata, unsigned char *idata,	int lensrc, int *width, int *height, int chksum );
};
#endif // WSQ31_H_
