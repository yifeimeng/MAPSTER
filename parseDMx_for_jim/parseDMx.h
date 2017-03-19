//
//  parseDMx.h
//
//
//  Created by Yifei Meng on 1/6/17.
//
//

#ifndef parseDMx_h
#define parseDMx_h

#include <stdio.h>
//function declarations
uint16_t swapBytesUnsigned16(uint16_t v);
uint32_t swapBytesUnsigned32(uint32_t v);
uint64_t swapBytesUnsigned64(uint64_t v);
void readTagData(FILE *fp, uint64_t *infoArray, uint64_t infoArraySize);
void readTag(FILE *fp, uint32_t versionNumber);
void readTagEntry(FILE *fp, uint32_t versionNumber, uint64_t tagIndex, uint32_t tagLevel, char *tagLabel, char **targetTagsLabels, uint16_t numberTargetTags, int64_t *targetTagsPositions);
void readTagGroup(FILE *fp, uint32_t versionNumber, uint32_t tagLevel, char *tagLabel, char **targetTagsLabels, uint16_t numberTargetTags, int64_t *targetTagsPositions);
void parseDMx(FILE *fp, int64_t *targetTagsPositions, char **targetTagsLabels, uint16_t numberTargetTags);

#endif /* parseDMx_h */
