/***************************************************************************************
    ipTest is example of using image processing (ip algorithm of IDS Core SDK in Linux
 **************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ip.h"
//#include "ma.h"

unsigned char buffer[1024 * 1024];
//Ip ip;

void usage()
{
    printf ("usage:\n");
    printf ("ipTest bmpFile templFile repeat_count\n");
    printf ("\tbmpFile - input BMP file name for processing\n");
    printf ("\ttemplFile - output template file name\n");
}

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
	usage();
	return 0;
    }
	
    FILE *bmp = 0, *tpl = 0;
    int result = 0;
    CaptureData finger;
    Ip ip;
    TemplateData templ;   
    
     // initialize ti image processor 
    if (ip.init() != IP_OK ||	ip.allocateTemplates(1, &templ) != MA_OK)
    {
	printf ("initialization error\n");
	return -1;
    }
    // read the BMP file
    bmp = fopen(argv[1], "rb");
    if ( bmp == NULL )
    {
	printf("Cannot open file %s\n\r", argv[1]);
	return -1;
    }
    fseek (bmp, 0, SEEK_END);
    unsigned int size = ftell (bmp);
    fseek (bmp, 14, SEEK_SET);
    if (fread(buffer, 1, size - 14, bmp) != size - 14)
    {
	printf("error read image file %s\n\r", argv[1]);
	fclose(bmp);
	return -1;
    }
    fclose(bmp);
    finger.dib = buffer;
    finger.numFinger = FINGPOS_UK;

    // process image
    int count = atoi(argv[3]);
    printf("Start processing %d times...\n",count);
    clock_t start = clock();
    for (int i = 0; i < count; i++) 
        result = ip.process(finger, templ);
    clock_t stop = clock();	
    if (result != IP_OK)
    {
	printf("Processing return error %d\n\r", result);
	return -1;
    }
    printf ("image processing time = %g sec\n", (double)(stop - start) / CLOCKS_PER_SEC / count );
    // write template to the template file 
    tpl = fopen(argv[2],"wb");
    if (!tpl)
    {
	printf("error open file %s for writing\n\r", argv[2]);
	return -1;
    }
    size = ((TemplHeader*)templ.fpTemplate)->size;
    if (fwrite(templ.fpTemplate, 1, size, tpl) != size)
    {
	printf("error write to file %s\n\r", argv[2]);
	fclose(tpl);
	return -1;
    }
    fclose(tpl);
    
    printf("template is written successfully\n\r");
    ip.freeTemplates (1, &templ);
    return 0;
}
