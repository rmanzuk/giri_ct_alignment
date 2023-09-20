# giri_ct_alignment
 Exploratory codes to try to align 3D image stacks captured separately through CT scanning and serial grinding and imaging

## General instructions and ideas so far
This core sample was CT scanned and then serially ground and imaged with GIRI at princeton. In order to understand the relationship between density returns from the CT stack and visible-light imagery, and use the two for a segmentation of the embedded fossils, we first need to align the two data sources so their voxels correspond. Although registration marks (scores) were made on the sample to help with alignment, simply using the outlines of the core with these scores for alignment does not yeild a unique, precise solution. I currenly opt for identifying several pyrite grains that are visible in both image stacks, getting their outlines through threshold segmentation, and using those points to align the stacks. This workflow takes the following steps, conducted in giri_ct_alignment_wrapper.m:

1. Set up basic directory structure for files and set the scales of pixels for each stack. Note the CT stack has perfectly cubic voxels, so it only has one voxel scale term, where the giri Z and pixel scales to not match, so those are placed in separate terms. Lines 1-34.
2. Becasue the core was in 2 pieces, we can't guarantee that we can perfectly align the whole thing at once in case the pieces were slightly offset for the two different imaging methods. So, I separate the ct and giri image directories at the approximate split between the two pieces to allow for separate alignments. Lines 35-41.
3. Use custom function crop_3d.m which asks for user inputs to get cropped volumes around several pyrite grains in each stack. Lines 46-57.
4. Use custom function outline_point_cloud.m to use thresholds and user inputs to get outline point clouds of each pyrite grain. Lines 58-69.
5. Use the position of the cropped volumes in the core to put the pyrite grain outlines into a single point cloud in absolute space. Lines 70-87.
6. Solve for the transform between the GIRI and CT pyrite point clouds. Lines 88-91.
7. Load in the CT imagery and apply the transform to get an aligned stack. Lines 92-100

Following that procedure, I've never gotten a perfect alignment between the two, so the lines after that in giri_ct_slignment_wrapper.m are some experiments to improve the fit after the initial transform or try other methods altogether. Those lines are worth perusing for ideas but might not yield much. 
