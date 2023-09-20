# giri_ct_alignment
 Exploratory codes to try to align 3D image stacks captured separately through CT scanning and serial grinding and imaging

## General instructions and ideas so far
This core sample was CT scanned and then serially ground and imaged with GIRI at princeton. In order to understand the relationship between density returns from the CT stack and visible-light imagery, and use the two for a segmentation of the embedded fossils, we first need to align the two data sources so their voxels correspond. Although registration marks (scores) were made on the sample to help with alignment, simply using the outlines of the core with these scores for alignment does not yeild a unique, precise solution. I currenly opt for identifying several pyrite grains that are visible in both image stacks, getting their outlines through threshold segmentation, and using those points to align the stacks. This workflow takes the following steps, conducted in giri_ct_alignment_wrapper.m:

1. 
