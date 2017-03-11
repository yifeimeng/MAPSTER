//
//  main.c
//
//
//  Created by Yifei Meng on 1/12/17.
//
//

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "hdf5.h"
#include "parseDMx.h"


#define NUM_TARGET_TAGS 5
#define DM4_FILE_NAME "/home/yifei/Documents/test_data/STOPTO-test.dm4"
#define HDF5_FILE_NAME_LENGTH 30
#define HDF5_FILE_NAME "hdf5_test"
/*
 * The buffer scale means the number of image pixels being read every time.
 * The buffer scale must be a divisor of the image (diffraction pattern) total pixel number.
 * The totol buffer size will be BUFFER_SCALE*total scanning pixels.
 * The I/O may be faster when the buffer size is power of 2
 * Choose the buffer size properly based on the system.
 * Stack memory is used for the reader buffer. Do not use a large number to avoid stack overflow.
 * 8 is a 512KB buffer. Everytime we call fread, a 512KB data is read.
 *
 */
#define READER_BUFFER_SCALE 8
/*
 *
 * The buffer scale means the number of pixels being written into every image every time.
 * The buffer scale must be a multiplier of the READER_BUFFER_SCALE.
 * The total size of the writer buffer is NUM_WRITER*numberImagePerWriter*WRITER_BUFFER_SCALE*4.
 * Dynamical memory is used for the writer buffer. Avoid allocating a heap larger than the system allowed.
 * 4096 is around 300MB here. Everytime we call H5Dwrite, around 16KB data is written.
 *
 */
#define WRITER_BUFFER_SCALE 4096
/*
 * The number of writer must be a divisior of the scanning pixel number.
 * More writer, smaller size for each file.
 */
#define NUM_WRITER 16



void raw2HDF5(int64_t *targetTagsPositions);

void raw2HDF5(int64_t *targetTagsPositions) {
    FILE *fp;
    fp = fopen(DM4_FILE_NAME, "rb");

    uint64_t sX, sY, dX, dY;//the index is stored in uint64_t
    uint16_t fourBytes;
    //read the scanning pixel number X
    fseek(fp, targetTagsPositions[1] + 20, SEEK_SET);//magic number here
    fread(&fourBytes, 4, 1, fp);
    sX = fourBytes;//implicit typecast
    printf("Number of scanning step at X: %llu\n", sX);
    //read the scanning pixel number Y
    fseek(fp, targetTagsPositions[2] + 20, SEEK_SET);//magic number here
    fread(&fourBytes, 4, 1, fp);
    sY = fourBytes;
    printf("Number of scanning step at Y: %llu\n", sY);
    //read the diffraction pattern pixel X
    fseek(fp, targetTagsPositions[3] + 20, SEEK_SET);//magic number here
    fread(&fourBytes, 4, 1, fp);
    dX = fourBytes;
    printf("Diffraciton pattern length at X: %llu\n", dX);
    //read the diffraction pattern pixel Y
    fseek(fp, targetTagsPositions[4] + 20, SEEK_SET);//magic number here
    fread(&fourBytes, 4, 1, fp);
    dY = fourBytes;
    printf("Diffraction pattern length at Y: %llu\n", dY);

    //construct the reader buffer. better to set the size as mutiplies of sX*sY
    uint64_t totalNumberImages = sX*sY;
    uint16_t readerBuffer[READER_BUFFER_SCALE][totalNumberImages];

    //construct the write buffer, the number of writer is same as the number of file generated
    uint64_t numberImagePerWriter = totalNumberImages/NUM_WRITER;
    uint16_t *writerBuffer = (uint16_t *)malloc(NUM_WRITER*numberImagePerWriter*WRITER_BUFFER_SCALE*sizeof(uint16_t));

    //move the file pointer to the starting position of the rawData
    uint64_t eightBytes;
    fseek(fp, targetTagsPositions[0] + 4, SEEK_SET);//another magic number,representing %%%%
    fread(&eightBytes, 8, 1, fp);
    printf("The number of elements in info array is: %llu\n",swapBytesUnsigned64(eightBytes));
    fread(&eightBytes, 8, 1, fp);
    fread(&eightBytes, 8, 1, fp);
    printf("The array datatype is: %llu\n", swapBytesUnsigned64(eightBytes));
    fread(&eightBytes, 8, 1, fp);
    printf("The number of elements in the array is: %llu\n", swapBytesUnsigned64(eightBytes));

    //perform the conversion

    //create the dataspace and the datatype
    hid_t datatype, dataspace;
    hsize_t dimsf[2];
    herr_t status;

    dimsf[0] = numberImagePerWriter;
    dimsf[1] = dX*dY;
    dataspace = H5Screate_simple(2, dimsf, NULL);

    datatype = H5Tcopy(H5T_STD_U16LE);
    status = H5Tset_order(datatype, H5T_ORDER_LE);

    //create the same number of hdf5 files and datasets as the writer
    hid_t hdf5FileSet[NUM_WRITER];
    hid_t hdf5DataSet[NUM_WRITER];
    char hdf5FileName[NUM_WRITER][HDF5_FILE_NAME_LENGTH];
    char indexString[5];
    for (uint64_t fileIndex = 0; fileIndex < NUM_WRITER; fileIndex++) {
        sprintf(indexString, "%04llu", fileIndex);
        snprintf(hdf5FileName[fileIndex], HDF5_FILE_NAME_LENGTH, "%s%c%s", HDF5_FILE_NAME, '_', indexString);
        hdf5FileSet[fileIndex] = H5Fcreate(hdf5FileName[fileIndex], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        char *datasetName = "DPStack";
        hdf5DataSet[fileIndex] = H5Dcreate2(hdf5FileSet[fileIndex], datasetName, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    //read the whole dm4 file
    uint64_t numberReads = (dX*dY)/READER_BUFFER_SCALE;
    uint64_t numberClearWriter = WRITER_BUFFER_SCALE/READER_BUFFER_SCALE;
    uint64_t writerFiller, fileFiller;

    for (uint64_t readIndex = 0; readIndex < numberReads; readIndex ++) {//mods
        //read a chunk of the data into the buffer
        fread(readerBuffer, sizeof(uint16_t), totalNumberImages*READER_BUFFER_SCALE, fp);
        printf("read segment %llu out of %llu\n", readIndex, numberReads);
        writerFiller = readIndex % numberClearWriter;
        fileFiller = readIndex / numberClearWriter;

        //transfer all data from the reader buffer into the writer buffer, fix any reshape/transpose
        uint64_t writerIndex, imageIndex;
        for (uint64_t bi = 0; bi < READER_BUFFER_SCALE; bi ++){
            for (uint64_t si = 0; si < totalNumberImages; si ++) {


                    uint64_t ptrPosition = si * WRITER_BUFFER_SCALE + writerFiller * READER_BUFFER_SCALE + bi;
                    writerBuffer[ptrPosition] = readerBuffer[bi][si];//perform the reshape
            }
        }

        //write the data in the writer buffer into the hdf files when the buffer is full
        if ((writerFiller + 1) == numberClearWriter) {

            for (uint64_t writerIndex = 0; writerIndex < 1; writerIndex ++) {//mods
                for (uint64_t imageIndex = 0; imageIndex < numberImagePerWriter; imageIndex ++) {

                    //set the writing block parameters
                    hsize_t offset[2] = {imageIndex, fileFiller*WRITER_BUFFER_SCALE};
                    hsize_t count[2] = {1, WRITER_BUFFER_SCALE};
                    hsize_t stride[2] = {1, 1};
                    hsize_t block[2] = {1, 1};

                    //create the memory sapce for the subset data
                    hsize_t dimsm[1];
                    dimsm[0] = WRITER_BUFFER_SCALE;
                    hid_t memspace;
                    memspace = H5Screate_simple(1, dimsm, NULL);
                    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, stride, count, block);

                    //write the data into the hdf5 file
                    uint16_t sdata[WRITER_BUFFER_SCALE];
                    uint64_t ptrShift = (writerIndex * numberImagePerWriter + imageIndex) * WRITER_BUFFER_SCALE;
                    memcpy(sdata, (writerBuffer + ptrShift), WRITER_BUFFER_SCALE*sizeof(uint16_t));
                    status = H5Dwrite(hdf5DataSet[writerIndex], datatype, memspace, dataspace, H5P_DEFAULT, sdata);

                    //close the memory space used
                    status = H5Sclose(memspace);

                }
            }

        }


    }

    //free the writer buffer
    free(writerBuffer);

    //close all the hdf5 files and corresponding datasets
    for (uint64_t fileIndex = 0; fileIndex < NUM_WRITER; fileIndex++) {

        status = H5Dclose(hdf5DataSet[fileIndex]);
        status = H5Fclose(hdf5FileSet[fileIndex]);
    }


    status = H5Sclose(dataspace);
    status = H5Tclose(datatype);



    //close the dm4 file
    fclose(fp);

}


int main() {
    FILE *fpDM4;
    char *targetTagsLabels[NUM_TARGET_TAGS];
    int64_t targetTagsPositions[NUM_TARGET_TAGS];

    //define the tages to be found
    targetTagsLabels[0] = ":>ImageList->1->ImageData->Data->";//the raw data tag
    targetTagsLabels[1] = ":>ImageList->1->ImageData->Dimensions->0->";
    targetTagsLabels[2] = ":>ImageList->1->ImageData->Dimensions->1->";
    targetTagsLabels[3] = ":>ImageList->1->ImageData->Dimensions->2->";
    targetTagsLabels[4] = ":>ImageList->1->ImageData->Dimensions->3->";

    targetTagsPositions[0] = -1;
    targetTagsPositions[1] = -1;
    targetTagsPositions[2] = -1;
    targetTagsPositions[3] = -1;
    targetTagsPositions[4] = -1;

    //find the data starting positions of the interested tages. The start posiition sits at the starting byte of the info array
    fpDM4 = fopen(DM4_FILE_NAME, "rb");
    if (fpDM4 == NULL)
        printf("No such file!!\n");
    else{
        parseDMx(fpDM4, targetTagsPositions, targetTagsLabels, NUM_TARGET_TAGS);
        fclose(fpDM4);
    }

    for (uint16_t targetTagIndex = 0; targetTagIndex < NUM_TARGET_TAGS; targetTagIndex++) {
        if (targetTagsPositions[targetTagIndex] == -1) {
            printf("Cannot find this tag!\n");
        }
        else{
            printf("Data of targe tag %d starts at %lld\n", targetTagIndex, targetTagsPositions[targetTagIndex]);
        }
    }

    //cut the raw data into separate hdf5 files. Fix the byte order.
    raw2HDF5(targetTagsPositions);



}
