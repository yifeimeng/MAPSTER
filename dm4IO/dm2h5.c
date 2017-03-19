#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "parseDMx.h"
#include "dm2h5.h"
#include "hdf5.h"

#define MAX_FILENAME_LENGTH 50
#define MAX_DIRECTORY_LENGTH 200
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
#define NUM_TAG 5


void dm2h5(char *input_dir, char *dmFileName, char *output_dir, uint32_t numWriter) {

    int64_t targetTagsPositions[NUM_TAG];
    char *targetTagsLabels[NUM_TAG];

    targetTagsLabels[0] = ":>ImageList->1->ImageData->Data->"; // the raw data tag
    targetTagsLabels[1] = ":>ImageList->1->ImageData->Dimensions->0->"; // scan x
    targetTagsLabels[2] = ":>ImageList->1->ImageData->Dimensions->1->"; // scan y
    targetTagsLabels[3] = ":>ImageList->1->ImageData->Dimensions->2->"; // dp size x
    targetTagsLabels[4] = ":>ImageList->1->ImageData->Dimensions->3->"; // dp size y

    targetTagsPositions[0] = -1;
    targetTagsPositions[1] = -1;
    targetTagsPositions[2] = -1;
    targetTagsPositions[3] = -1;
    targetTagsPositions[4] = -1;

    // find the starting positions of the interested tages. The start position sits at the starting byte of the info array
    FILE *fp;
    char openFileName[MAX_DIRECTORY_LENGTH + MAX_FILENAME_LENGTH];
    strcpy(openFileName, input_dir);
    strcat(openFileName, dmFileName);
    strcat(openFileName, ".dm4");
    printf("File to be opened: %s\n", openFileName);

    fp = fopen(openFileName, "rb");
    if (fp == NULL) {
        printf("No such file!!\n");
        return;
    }
    else{
        parseDMx(fp, targetTagsPositions, targetTagsLabels, NUM_TAG);
        fclose(fp);
    }

    for (uint16_t targetTagIndex = 0; targetTagIndex < NUM_TAG; targetTagIndex++) {
        if (targetTagsPositions[targetTagIndex] == -1) {
            printf("Cannot find this tag!\n");
        }
        else{
            printf("Data of target tag %d starts at %ld\n", targetTagIndex, targetTagsPositions[targetTagIndex]);
        }
    }
    printf("DM file parsing finished.\n");
    printf("\n");




    fp = fopen(openFileName, "rb");

    uint64_t sX, sY, dX, dY;//the index is stored in uint64_t
    uint16_t fourBytes;
    //read the scanning pixel number X
    fseek(fp, targetTagsPositions[1] + 20, SEEK_SET);//magic number here
    fread(&fourBytes, 4, 1, fp);
    sX = fourBytes;//implicit typecast
    printf("Number of scanning step at X: %lu\n", sX);
    //read the scanning pixel number Y
    fseek(fp, targetTagsPositions[2] + 20, SEEK_SET);//magic number here
    fread(&fourBytes, 4, 1, fp);
    sY = fourBytes;
    printf("Number of scanning step at Y: %lu\n", sY);
    //read the diffraction pattern pixel X
    fseek(fp, targetTagsPositions[3] + 20, SEEK_SET);//magic number here
    fread(&fourBytes, 4, 1, fp);
    dX = fourBytes;
    printf("Diffraciton pattern length at X: %lu\n", dX);
    //read the diffraction pattern pixel Y
    fseek(fp, targetTagsPositions[4] + 20, SEEK_SET);//magic number here
    fread(&fourBytes, 4, 1, fp);
    dY = fourBytes;
    printf("Diffraction pattern length at Y: %lu\n", dY);

    //construct the reader buffer. better to set the size as mutiplies of sX*sY
    uint64_t totalNumberImages = sX*sY;
    uint16_t readerBuffer[READER_BUFFER_SCALE][totalNumberImages];

    //construct the write buffer, the number of writer is same as the number of file generated
    uint64_t numberImagePerWriter = totalNumberImages/numWriter;
    uint16_t *writerBuffer = (uint16_t *)malloc(numWriter*numberImagePerWriter*WRITER_BUFFER_SCALE*sizeof(uint16_t));

    //move the file pointer to the starting position of the rawData
    uint64_t eightBytes;
    fseek(fp, targetTagsPositions[0] + 4, SEEK_SET);// skip the %%%%
    fread(&eightBytes, 8, 1, fp);
    printf("The number of elements in info array is: %lu\n",swapBytesUnsigned64(eightBytes));
    fread(&eightBytes, 8, 1, fp);
    fread(&eightBytes, 8, 1, fp);
    printf("The array datatype is: %lu\n", swapBytesUnsigned64(eightBytes));
    fread(&eightBytes, 8, 1, fp);
    printf("The number of elements in the array is: %lu\n", swapBytesUnsigned64(eightBytes));


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
    hid_t hdf5FileSet[numWriter];
    hid_t hdf5DataSet[numWriter];
    char hdf5FileName[numWriter][MAX_FILENAME_LENGTH + MAX_DIRECTORY_LENGTH];
    char indexString[5];
    for (uint64_t fileIndex = 0; fileIndex < numWriter; fileIndex++) {
        sprintf(indexString, "%04lu", fileIndex);
        snprintf(hdf5FileName[fileIndex], MAX_FILENAME_LENGTH + MAX_DIRECTORY_LENGTH, "%s%s%c%s%s", output_dir, dmFileName, '_', indexString, ".h5");
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
        printf("read segment %lu out of %lu\n", readIndex, numberReads);
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
    for (uint64_t fileIndex = 0; fileIndex < numWriter; fileIndex++) {

        status = H5Dclose(hdf5DataSet[fileIndex]);
        status = H5Fclose(hdf5FileSet[fileIndex]);
    }


    status = H5Sclose(dataspace);
    status = H5Tclose(datatype);


    //close the dm4 file
    fclose(fp);

}

int main(int argc, char *argv[]) {

    char dmFileName[MAX_FILENAME_LENGTH];
    char input_dir[MAX_DIRECTORY_LENGTH];
    char output_dir[MAX_DIRECTORY_LENGTH];

    if (argc == 1) {
        printf("At least one argument (input file name) is expected.\n");
        return 0;
    }

    //default number of writer is 1
    uint32_t numWriter = 1, numSwitch = 0;
    output_dir[0] = '\0';

    // examine the switch
    for (uint32_t i = 1; i < argc; i++) {
        /* Check for a switch (leading "-"). */
        if (argv[i][0] == '-') {
            /* Use the next character to decide what to do. */
            switch (argv[i][1]) {
                case 'o':
                        printf("flag o was called\n");
                        numSwitch ++;
                        strcpy(output_dir, argv[i + 1]);
                        break;
                case 'w':
                        printf("flag w was called\n");
                        numSwitch ++;
                        sscanf(argv[i + 1], "%u", &numWriter);
                        break;
                default:
                        printf("unrecognized flag %c\n", argv[i][1]);
                        break;
            }
        }
    }

    if (argc > (2 + numSwitch * 2)) {
        printf("Unknown arguments.\n");
        return 0;
    }

    // examine the input file argument
    char *pch;
    pch = strrchr(argv[1], '.');
    if (strcmp(pch, ".dm4") != 0) {
        printf("This is not a dm4 file.\n");
        return 0;
    }

    // separate the file name and directory
    char *pchDot, *pchSlash;
    pchDot = strrchr(argv[1], '.');
    pchSlash = strrchr(argv[1], '/');

    if (pchSlash != NULL) {
        // extract the file name
        int fileNameLength = pchDot - pchSlash - 1;
        for (uint32_t i = 0; i < fileNameLength; i ++) {
            dmFileName[i] = pchSlash[i + 1];
        }

        dmFileName[fileNameLength] = '\0';
        printf("The input DM file name is: %s\n", dmFileName);

        // extract the directory
        int directoryLength = pchSlash - argv[1] + 1;
        for (uint32_t i = 0; i < directoryLength; i ++) {
            input_dir[i] = argv[1][i];

        }
        printf("The input directory is: %s\n", input_dir);
    }
    else {

        strcpy(dmFileName, argv[1]);
        strcpy(input_dir, "./");

    }

    if (strcmp(output_dir, "\0") == 0) {
        strcpy(output_dir, input_dir);

    }

    dm2h5(input_dir, dmFileName, output_dir, numWriter);

}
