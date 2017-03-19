Subroutine:
===============================
void parseDMx(FILE *fp, int64_t *targetTagsPositions, char **targetTagsLabels, uint16_t numberTargetTags)
===============================

This subroutine parses both dm3 and dm4 files. Define the tags you are interested in the "targetTagsLabels". 
Also assign the number of tags to be searched to "numberTargetTags". 
The search results are stored in the "targetTagsPositions". See sample.c for an example. 

Parameters:
fp: The DigitalMicrograph file.
targetTagsPositions: Output results are stored here. The start positions of interested tags. They are not the starting positions of raw data. The info array must be skipped to reach the real raw data. See the sample for details.
targetTagsLabel: The name of interested tags. In most cases, the raw image data is stored in ":>ImageList->1->ImageData->Data->".
numberTargetTags: Number of tags to be searched.

Check http://www.er-c.org/cbb/info/dmformat/ for more information about dm file format.

