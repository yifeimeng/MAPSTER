//
//  sample.c
//
//  This sample shows how to extract the raw binary data from a dm3 image file. 
//
//  The tag and info array are always wriiten in big-endianness.
//  The raw data's endianness depends on the PC that writes it.
//
//  Check http://www.er-c.org/cbb/info/dmformat/ for more information about dm file format.
// 
//  Created by Yifei Meng on 1/12/17.
//
//

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "parseDMx.h"


#define NUM_TARGET_TAGS 3
#define DM_FILE_NAME "400A_1003.dm3"


int main() {
    FILE *fp, *fp_output;
    char *targetTagsLabels[NUM_TARGET_TAGS];
    int64_t targetTagsPositions[NUM_TARGET_TAGS];

    // define the tags to be found. 
    // In most cases, the raw data is always stored in ":>ImageList->1->ImageData->Data->".
    // the image size x is stored in ":>ImageList->1->ImageData->Dimensions->0->"
    // the image size y is stored in ":>ImageList->1->ImageData->Dimensions->1->"
    // if the data is a 4D dataset, for example, scanning diffraction stack, 
    // dimension0, 1, 2, 3 will store scanning size x and y, DP size x and y respectively 		
    // ImageList->0 often contains a thumbnail of the original image
    targetTagsLabels[0] = ":>ImageList->1->ImageData->Data->"; // the raw data tag
    targetTagsLabels[1] = ":>ImageList->1->ImageData->Dimensions->0->"; // image size x 
    targetTagsLabels[2] = ":>ImageList->1->ImageData->Dimensions->1->"; // image size y


    targetTagsPositions[0] = -1;
    targetTagsPositions[1] = -1;
    targetTagsPositions[2] = -1;

    // find the starting positions of the interested tages. The start position sits at the starting byte of the info array
    fp = fopen(DM_FILE_NAME, "rb");
    if (fp == NULL)
        printf("No such file!!\n");
    else{
        parseDMx(fp, targetTagsPositions, targetTagsLabels, NUM_TARGET_TAGS);
        fclose(fp);
    }

    for (uint16_t targetTagIndex = 0; targetTagIndex < NUM_TARGET_TAGS; targetTagIndex++) {
        if (targetTagsPositions[targetTagIndex] == -1) {
            printf("Cannot find this tag!\n");
        }
        else{
            printf("Data of targe tag %d starts at %lld\n", targetTagIndex, targetTagsPositions[targetTagIndex]);
        }
    }
    printf("DM file parsing finished.\n");
    printf("\n");

    
    // print the information acquired and determine the raw data type/position
    fp = fopen(DM_FILE_NAME, "rb"); // be sure to open the file in the binary mode
    uint32_t sX, sY;//the index is stored in uint32_t
    // read the scanning pixel number X
    fseek(fp, targetTagsPositions[1] + 12, SEEK_SET);// shift 12 bytes for dm3, 20 bytes for dm4
    fread(&sX, 4, 1, fp);
    printf("Number of pixels at X: %u\n", sX);	

    //read the scanning pixel number Y
    fseek(fp, targetTagsPositions[2] + 12, SEEK_SET);// shift 12 bytes for dm3, 20 bytes for dm4
    fread(&sY, 4, 1, fp);
    printf("Number of pixels at Y: %u\n", sY);	

    printf("\n"); 	
    // move the file pointer to the starting position of the rawData tag, read the info array
    // for dm3 file, info array element is 4 bytes integer(BE), for dm4 file, it is 8 bytes 
    uint32_t fourBytes;
    fseek(fp, targetTagsPositions[0] + 4, SEEK_SET);// skip the %%%%, shift 4 bytes
    fread(&fourBytes, 4, 1, fp); // read the size of the info array, should be 0x3
    fread(&fourBytes, 4, 1, fp); // should be 0x14, which means its an array
    fread(&fourBytes, 4, 1, fp); // read the element type

    /* Check this table for the data type of the image	
     02h =  2  i2 signed    (short)
     03h =  3  i4 signed    (long)
     04h =  4  i2 unsigned  (ushort) or unicode string
     05h =  5  i4 unsigned  (ulong)
     06h =  6  f4           (float)
     07h =  7  f8           (double)
     08h =  8  i1            (boolean)
     09h =  9  a1            (char)
     0ah = 10  i1
     0bh = 11  i8  ?         (long long) not sure if signed or unsigned

     In most cases, DM stores the data in unsigned 16-bit or 32-bit integer.
    */


    printf("The image datatype is: %u\n", swapBytesUnsigned32(fourBytes));
    fread(&fourBytes, 4, 1, fp);
    printf("The number of pixels in the image is: %u\n", swapBytesUnsigned32(fourBytes));	

    printf("The binary data starts at %u.\n", targetTagsPositions[0] + 20);	

    // read the raw data and perform the processing, here we simply store the raw data in a binary format
    // The endianness of the raw data depends on the computer that writes it. 
    // Most modern PC use Intel x86 CPU, which means little endianness.
    fseek(fp, targetTagsPositions[0] + 20, SEEK_SET); // shift 20 bytes for dm3, 36 bytes for dm4
    uint32_t *image = (uint32_t *)malloc(sX*sY*sizeof(uint32_t));
    fread(image, 4, sX*sY, fp);
    	
    fp_output = fopen("image.bin", "wb");
    fwrite(image, 4, sX*sY, fp_output);
	

    fclose(fp);	 
    fclose(fp_output);
    free(image);
	

}
