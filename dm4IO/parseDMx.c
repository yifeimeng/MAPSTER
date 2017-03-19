//
//  parseDMx.c
//
//  This subroutine parses the dm3 or dm4 files and finds the starting position of raw data in the file.
//
//  Created by Yifei Meng on 1/6/17.
//
//

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "parseDMx.h"


uint16_t swapBytesUnsigned16(uint16_t v){
    return (v << 8) | (v >> 8);
}

uint32_t swapBytesUnsigned32(uint32_t v){
    uint32_t s;
    s = ((v>>24)&0xff) | ((v>>8)&0xff00) | ((v<<8)&0xff0000) | ((v<<24)&0xff000000);
    return s;
}

uint64_t swapBytesUnsigned64(uint64_t v){
    uint64_t s;
    s = ((v << 8) & 0xff00ff00ff00ff00ULL ) | ((v >> 8) & 0x00ff00ff00ff00ffULL );
    s = ((s << 16) & 0xffff0000ffff0000ULL ) | ((s >> 16) & 0x0000ffff0000ffffULL );
    return (s << 32) | (s >> 32);
}

//read specific type of data. Define the infoArray when use this function recursively.
void readTagData(FILE *fp, uint64_t *infoArray, uint64_t infoArraySize) {

    //determine the data type of the tag
    uint64_t tagType = infoArray[0];
    //printf("data type is 0x%llx:",tagType); // used for debug

    //process the tag based on different data type
    switch (tagType) {
        case 0x0f: {//the tag is a structure(group)
            printf("struct-");
            uint64_t numberEntry = infoArray[2];
            uint64_t entryIndex, entryDataType;
            for (entryIndex = 0; entryIndex < numberEntry; entryIndex++) {
                entryDataType = infoArray[2*entryIndex+4];
                printf("(0x%lx):",entryDataType);
                uint64_t tempInfoArray[1] = {infoArray[2*entryIndex+4]};
                readTagData(fp, tempInfoArray, 1);

            }

            break;
        }

        case 0x14: {//the tag is an array
            if (infoArraySize == 3) {//array of data
                uint64_t elementDataType = infoArray[1];
                uint64_t numberElement = infoArray[2];
                printf("array of 0x%lx, 0x%lx elements",elementDataType, numberElement);
                uint64_t tempInfoArray[1] = {infoArray[1]};
                for (uint64_t elementIndex = 0; elementIndex < numberElement; elementIndex++){
                    readTagData(fp, tempInfoArray, 1);
                }
            }
            else {//arrary of structure (AoS)
                uint64_t elementDataType = infoArray[1];
                if (elementDataType == 0x0f) {//comfirm this is an AoS
                    uint64_t numberEntry = infoArray[3];
                    uint64_t numberElement = infoArray[infoArraySize-1];
                    uint64_t *tempInfoArray = &infoArray[1];
                    for (uint64_t elementIndex = 0; elementIndex < numberElement; elementIndex++) {
                        readTagData(fp, tempInfoArray, infoArraySize - 2);
                    }
                }
            }
            break;
        }
        case 0x12:{//the tag is a string
            //under construction
            uint32_t stringLength;
            fread(&stringLength, 4, 1, fp);
            char currString[stringLength];
            fread(&currString, 1, stringLength, fp);
            puts(currString);
            break;
        }
        case 0x02: {//the tag is a 2 bytes signed integer
            signed short currValue;
            fread(&currValue, 2, 1, fp);
            //printf("%d",currValue);
            break;
        }
        case 0x03: {//the tag is a 4 bytes signed integer
            signed int currValue;
            fread(&currValue, 4, 1, fp);
            //printf("%d",currValue);
            break;
        }
        case 0x04: {//the tag is a 2 bytes uint16 integer or an Unicode character (UTF-16)
            uint16_t currValue;
            fread(&currValue, 2, 1, fp);
            // printf("%d",currValue);
            break;
        }
        case 0x05: {//the tag is a 4 bytes uint32 integer
            uint32_t currValue;
            fread(&currValue, 4, 1, fp);
            //printf("%d",currValue);
            break;
        }
        case 0x06: {//the tag is a 4 bytes float number
            float currValue;
            fread(&currValue, 4, 1, fp);
            //printf("%f",currValue);
            break;
        }
        case 0x07: {//the tag is a 8 bytes double-precesion float number
            double currValue;
            fread(&currValue, 8, 1, fp);
            //printf("%f",currValue);
            break;
        }
        case 0x08: {//the tag is a 1 byte boolean number
            char currValue;
            fread(&currValue, 1, 1, fp);
            //printf("%d",currValue);
            break;
        }
        case 0x09: {//the tag is a 1 byte char
            char currValue;
            fread(&currValue, 1, 1, fp);
            //printf("%c",currValue);
            break;
        }
        case 0x0a:{//the tag is a 1 byte integer,
            char currValue;
            fread(&currValue, 1, 1, fp);
            //printf("%d", currValue);
            break;
        }
        case 0x0b: {//the tag is a 8 bytes integer, not clear signed or unsigned
            uint64_t currValue;
            fread(&currValue, 8, 1, fp);
            //printf("%llu",currValue);
            break;
        }
        case 0x0c: {//the tag is a 8 bytes unknown data type
            uint64_t currValue;
            fread(&currValue, 8, 1, fp);
            //printf("%llu",currValue);
            break;
        }
        default:
            //printf("Unknown data type!");
            break;
    }
    //printf("\n");
}

void readTag(FILE *fp, uint32_t versionNumber) {
    //read the %%%%
    char sepString[5];
    fread(sepString, 1, 4, fp);
    sepString[4] = '\0';

    //read the info array
    uint64_t eightBytes, infoArraySize;
    uint32_t fourBytes;

    if (versionNumber == 4) {
        // read the info array size
        fread(&eightBytes, 8, 1, fp);
        infoArraySize = swapBytesUnsigned64(eightBytes);
        //read the whole info array
        uint64_t infoArray[infoArraySize];
        fread(infoArray, 8, infoArraySize, fp);
        //swap the bytes to little endian for elements in the info array
        for (uint64_t i = 0; i < infoArraySize; i ++) {
            infoArray[i] = swapBytesUnsigned64(infoArray[i]);
        }

	//read the data stored in the tag
        readTagData(fp, infoArray, infoArraySize);

    }
    else if (versionNumber == 3) {
        fread(&fourBytes, 4, 1, fp);
        infoArraySize = (uint64_t)swapBytesUnsigned32(fourBytes);
	printf("info array size is %lu.\n", infoArraySize);
        uint64_t infoArray[infoArraySize];
        // read the 4-byte integer one by one, and convert it into 8-byte
        for (uint64_t i = 0; i < infoArraySize; i ++) {
            fread(&fourBytes, 4, 1, fp);
            infoArray[i] = (uint64_t)swapBytesUnsigned32(fourBytes);
        }

	//read the data stored in the tag
        readTagData(fp, infoArray, infoArraySize);

    }




    printf("\n");

}


void readTagEntry(FILE *fp, uint32_t versionNumber, uint64_t tagIndex, uint32_t tagLevel, char *tagLabel, char **targetTagsLabels, uint16_t numberTargetTags, int64_t *targetTagsPositions) {
    char tagType;
    uint16_t twoBytes, currLevelLabelLength;
    uint64_t eightBytes, childrenTagSize;
    uint32_t tagLabelLength;

    tagLabelLength = strlen(tagLabel);
    tagLevel++;
    //printf("tag level is %d\n", tagLevel);
    fread(&tagType, 1, 1, fp);
    printf("tag type is 0x%x.\n",tagType);
    fread(&twoBytes, 2, 1, fp);
    currLevelLabelLength = swapBytesUnsigned16(twoBytes);
    printf("tag label length is %d.\n",currLevelLabelLength);
    if (currLevelLabelLength>0) {
        char currLevelLabel[currLevelLabelLength+1];
        fread(currLevelLabel, 1, currLevelLabelLength, fp);//copy the tag of current level
        currLevelLabel[currLevelLabelLength] = '\0';//put the end mark
        strcat(tagLabel,currLevelLabel);//extend the total tag label
        strcat(tagLabel,"->");//add the indication of the next level
        printf("%s\n",tagLabel);
    }
    else {
        //tag name doesn't exist.
        char tagIndexLabel[10];//use the current tag index as the name of this unnamed tag
        sprintf(tagIndexLabel, "%lu", tagIndex);
        strcat(tagLabel,tagIndexLabel);
        strcat(tagLabel,"->");//add the indication of the next level
        printf("%s\n",tagLabel);
    }

    //read the total bytes in the current tag directory including all children directories
    //this is new for dm4 files, no such information in dm3 files
    if (versionNumber == 4) {
        fread(&eightBytes, 8, 1, fp);
        childrenTagSize = swapBytesUnsigned64(eightBytes);
    }

    //0x14 indicates this is a tag directory, 0x15 indicates this is a tag
    if (tagType == 0x14) {
        readTagGroup(fp, versionNumber, tagLevel, tagLabel, targetTagsLabels, numberTargetTags, targetTagsPositions);
    }
    else if (tagType == 0x15) {
        //check if this tag is one of the target tags
        char isTargetTag = 0;
        for (uint16_t targetTagIndex = 0; targetTagIndex < numberTargetTags; targetTagIndex++) {
            if (strcmp(tagLabel, targetTagsLabels[targetTagIndex]) == 0) {
                targetTagsPositions[targetTagIndex] = ftell(fp);
                isTargetTag = 1;
                //break;
            }
        }

        if (isTargetTag == 1) {

            if (versionNumber == 4) {
                // skip the whole target tag for dm4 files
                // avoid reading all data into the memory for dm4 files
                fseek(fp, childrenTagSize, SEEK_CUR);
            }
            else if (versionNumber == 3) {
                readTag(fp, versionNumber);

            }


        }
        else {

            readTag(fp, versionNumber);
        }
    }
    else {

    }

    tagLabel[tagLabelLength] = '\0';//convert it to the original tag label
}

void readTagGroup(FILE *fp, uint32_t versionNumber, uint32_t tagLevel, char *tagLabel, char **targetTagsLabels, uint16_t numberTargetTags, int64_t *targetTagsPositions) {
    char oneByte, isSort, isOpen;
    uint64_t eightBytes;
    uint32_t fourBytes;

    fread(&oneByte, 1, 1 ,fp);
    isSort = oneByte;
    if (isSort == 1)
        printf("This tag directory is sorted.\n");
    else
        printf("This tag directory is not sorted.\n");

    fread(&oneByte, 1, 1, fp);
    isOpen = oneByte;
    if (isOpen == 1)
        printf("This tag directory is open.\n");
    else
        printf("This tag directory is closed.\n");

    //read number of tags in the current directory
    uint64_t numberTags;
    if (versionNumber == 4) {
        fread(&eightBytes, 8, 1, fp);
        numberTags = swapBytesUnsigned64(eightBytes);
    }
    else if (versionNumber == 3) {
        fread(&fourBytes, 4, 1, fp);
        numberTags = (uint64_t)swapBytesUnsigned32(fourBytes);
    }
    printf("%lu tags in this directory.\n",numberTags);

    //read all tags
    for (uint64_t tagIndex = 0;tagIndex < numberTags;tagIndex++) {
        //printf("tag level is %d.\n",tagLevel);
        readTagEntry(fp, versionNumber, tagIndex, tagLevel, tagLabel, targetTagsLabels, numberTargetTags, targetTagsPositions);
    }
}

void parseDMx(FILE *fp, int64_t *targetTagsPositions, char **targetTagsLabels, uint16_t numberTargetTags) {
    uint32_t fourBytes;
    uint64_t eightBytes;

    //read the version number, 3 or 4
    uint32_t versionNumber;
    fread(&fourBytes, 4, 1, fp);
    versionNumber = swapBytesUnsigned32(fourBytes);
    printf("The version number is %d\n", versionNumber);

    //read the total file size in bytes, it is the file length - 24 bytes (dm4) / -16 bytes (dm3)
    uint64_t fileSize;
    if (versionNumber == 4) {
    	fread(&eightBytes, 8, 1, fp);
    	fileSize = swapBytesUnsigned64(eightBytes);
    }
    else if (versionNumber == 3) {
        fread(&fourBytes, 4, 1, fp);
        fileSize = (uint64_t)swapBytesUnsigned32(fourBytes);
    }
    else {
        printf("This file is not supported!\n");
        return;
    }
    printf("The file size is %lu bytes\n", fileSize);

    //read the order of bytes, 1 is little endian (LE), 0 is big endian (BE)
    uint32_t byteOrder;
    fread(&fourBytes, 4, 1, fp);
    byteOrder = swapBytesUnsigned32(fourBytes);
    if (byteOrder == 1) {
        printf("The dataset byte order is little endian.\n");}
    else {
        printf("The dataset byte order is big endian.\n");}

    //read all tags, starting from the root directory
    char tagLabel[200] = ":>";//initiate the tag label
    puts(tagLabel);
    uint32_t initialTagLevel = 0;
    readTagGroup(fp, versionNumber, initialTagLevel, tagLabel, targetTagsLabels, numberTargetTags, targetTagsPositions);

}
