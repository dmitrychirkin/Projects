#include "stdafx.h"
#include "swsq31.h"
#include "wsq.h"

void WSQ_EI *init_instance()
{
	Wsq *W = new Wsq();
	return (void *)W;
}

void WSQ_EI cancel_instance(void *pInstance)
{ 
	Wsq *W = (Wsq *)pInstance;
	delete W;
}

int WSQ_EI wsq_encode(void *pInstance,unsigned char *idata, unsigned char *odata, int width, int height, float ratio, int chksum)
{
	Wsq *W = (Wsq *)pInstance;
	return W->encode(idata, odata, width, height, ratio, chksum);
}

int WSQ_EI wsq_decode(void *pInstance, unsigned char *odata, unsigned char *idata, int ilen, int *cols, int *rows, int chksum)
{
	Wsq *W = (Wsq *)pInstance;
	int cls, rws;
	int ans = W->decode(odata, idata, ilen, &cls, &rws, chksum);
	*cols = cls;
	*rows = rws;
	return ans;
}

