//
//  convertDM4ToHDF5.c
//  
//
//  Created by Yifei Meng on 1/6/17.
//
//

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "convertDM4ToHDF5.h"

//function declarations
void readTagGroup(FILE *fp, uint32_t tagLevel, char *tagLabel);
void readTagData(FILE *fp, uint64_t *infoArray, uint64_t infoArraySize);

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
    uint64_t tagType = swapBytesUnsigned64(infoArray[0]);
    printf("data type is 0x%llx:",tagType);
    
    //process the tag based on different data type
    switch (tagType) {
        case 0x0f: {//the tag is a structure(group)
            printf("struct-");
            uint64_t numberEntry = swapBytesUnsigned64(infoArray[2]);
            uint64_t entryIndex, entryDataType;
            for (entryIndex = 0; entryIndex < numberEntry; entryIndex++) {
                entryDataType = swapBytesUnsigned64(infoArray[2*entryIndex+4]);
                printf("(0x%llx):",entryDataType);
                uint64_t tempInfoArray[1] = {infoArray[2*entryIndex+4]};
                readTagData(fp, tempInfoArray, 1);
                
            }
            
            break;
        }
            
        case 0x14: {//the tag is an array of data or an array of structure (AoS)
            if (infoArraySize == 3) {//array of data
                uint64_t elementDataType = swapBytesUnsigned64(infoArray[1]);
                uint64_t numberElement = swapBytesUnsigned64(infoArray[2]);
                printf("array of 0x%llx, 0x%llx elements",elementDataType, numberElement);
                uint64_t tempInfoArray[1] = {infoArray[1]};
                for (uint64_t elementIndex = 0; elementIndex < numberElement; elementIndex++){
                    readTagData(fp, tempInfoArray, 1);
                }
            }
            else {//AoS
                uint64_t elementDataType = swapBytesUnsigned64(infoArray[1]);
                if (elementDataType == 0x0f) {//comfirm this is an AoS
                    uint64_t numberEntry = swapBytesUnsigned64(infoArray[3]);
                    uint64_t numberElement = swapBytesUnsigned64(infoArray[infoArraySize-1]);
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
            printf("%d",currValue);
            break;
        }
        case 0x03: {//the tag is a 4 bytes signed integer
            signed int currValue;
            fread(&currValue, 4, 1, fp);
            printf("%d",currValue);
            break;
        }
        case 0x04: {//the tag is a 2 bytes uint32_teger or an Unicode character (UTF-16)
            uint16_t currValue;
            fread(&currValue, 2, 1, fp);
            printf("%d",currValue);
            break;
        }
        case 0x05: {//the tag is a 4 bytes uint32_teger
            uint32_t currValue;
            fread(&currValue, 4, 1, fp);
            printf("%d",currValue);
            break;
        }
        case 0x06: {//the tag is a 4 bytes float number
            float currValue;
            fread(&currValue, 4, 1, fp);
            printf("%f",currValue);
            break;
        }
        case 0x07: {//the tag is a 8 bytes double-precesion float number
            double currValue;
            fread(&currValue, 8, 1, fp);
            printf("%f",currValue);
            break;
        }
        case 0x08: {//the tag is a 1 byte boolean number
            char currValue;
            fread(&currValue, 1, 1, fp);
            printf("%d",currValue);
            break;
        }
        case 0x09: {//the tag is a 1 byte char
            char currValue;
            fread(&currValue, 1, 1, fp);
            printf("%c",currValue);
            break;
        }
        case 0x0a:{//the tag is a 1 byte integer,
            char currValue;
            fread(&currValue, 1, 1, fp);
            printf("%d", currValue);
            break;
        }
        case 0x0b: {//the tag is a 8 bytes integer, not clear signed or unsigned
            uint64_t currValue;
            fread(&currValue, 8, 1, fp);
            printf("%llu",currValue);
            break;
        }
        case 0x0c: {//the tag is a 8 bytes unknown data type
            uint64_t currValue;
            fread(&currValue, 8, 1, fp);
            printf("%llu",currValue);
            break;
        }
        default:
            printf("Unknown data type!");
            break;
    }
    printf("\n");
}

void readTag(FILE *fp) {
    //read the %%%%
    char sepString[5];
    fread(sepString, 1, 4, fp);
    sepString[4] = '\0';
    
    //read the size of the info array
    uint64_t eightBytes, infoArraySize;
    fread(&eightBytes, 8, 1, fp);
    infoArraySize = swapBytesUnsigned64(eightBytes);
    //printf("info array size is %llu\n",infoArraySize);
    
    //read the whole info array
    uint64_t infoArray[infoArraySize];
    fread(infoArray, 8, infoArraySize, fp);
    
    //read the data stored in the tag
    readTagData(fp, infoArray, infoArraySize);
    
    
    printf("\n");
}


void readTagEntry(FILE *fp, uint32_t tagLevel, char *tagLabel) {
    char tagType;
    uint16_t twoBytes, currLevelLabelLength;
    uint64_t eightBytes, chillurenTagSize;
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
        strcat(tagLabel,"->");//add the indication of the next level
        printf("%s\n",tagLabel);
    }
    
    //read the total bytes in the current tag directory including all chilluren directories
    //this is new for dm4 files, no such fiellu in dm3 files
    fread(&eightBytes, 8, 1, fp);
    chillurenTagSize = swapBytesUnsigned64(eightBytes);
    
    //0x14 indicates this is a tag directory, 0x15 indicates this is a tag
    if (tagType == 0x14) {
        readTagGroup(fp, tagLevel, tagLabel);
    }
    else if (tagType == 0x15) {
        readTag(fp);
    }
    else {
        
    }
    
    tagLabel[tagLabelLength] = '\0';//convert it to the original tag label
}

void readTagGroup(FILE *fp, uint32_t tagLevel, char *tagLabel) {
    char oneByte, isSort, isOpen;
    uint64_t eightBytes;
    
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
    fread(&eightBytes, 8, 1, fp);
    numberTags = swapBytesUnsigned64(eightBytes);
    printf("%llu tags in this directory.\n",numberTags);
    
    //read all tags
    uint64_t tagIndex;
    for (tagIndex = 0;tagIndex < numberTags;tagIndex++) {
        //printf("tag level is %d.\n",tagLevel);
        readTagEntry(fp, tagLevel, tagLabel);
    }
}

void readDM4(FILE *fp) {
    uint32_t fourBytes;
    uint64_t eightBytes;
    
    //read the version number, 3 or 4
    uint32_t versionNumber;
    fread(&fourBytes, 4, 1, fp);
    versionNumber = swapBytesUnsigned32(fourBytes);
    printf("The version number is %d\n", versionNumber);
    
    //read the total file size in bytes, it is the file length - 24 bytes
    uint64_t fileSize;
    fread(&eightBytes, 8, 1, fp);
    fileSize = swapBytesUnsigned64(eightBytes);
    printf("The file size is %llu bytes\n", fileSize);
    
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
    readTagGroup(fp, initialTagLevel, tagLabel);
    
}

int main() {
    FILE *targetFP;
    
    targetFP = fopen("./imageTest.dm4","rb");
    if (targetFP == NULL)
        printf("No such file!!\n");
    else{
        readDM4(targetFP);
        fclose(targetFP);
    }
}