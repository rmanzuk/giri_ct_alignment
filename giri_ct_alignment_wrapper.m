% Wrapper for to take a CT stack and resample it to create a stack aligned
% with GIRI's outputs
%
% Ryan A. Manzuk
% Created: 08/01/2022
% Last edited: 09/12/2022
%% Setup
% paths for images
ct_folder = '/Users/ryan/Dropbox (Princeton)/giri_ct_alignment_project/siberia_1ab/ct_resample1';
giri_folder = '/Users/ryan/Dropbox (Princeton)/giri_ct_alignment_project/siberia_1ab/giri_downsampled/625nm';

% make those paths into nice, sorted directories
ct_ext = '*.tif';
ct_dir = dir(fullfile(ct_folder,ct_ext));
[~,ct_idx,~] = natsortfiles({ct_dir.name});
ct_dir = ct_dir(ct_idx);

giri_ext = '*.jpg';
giri_dir = dir(fullfile(giri_folder,giri_ext));
sorted_idx = im_num_sort(giri_folder, giri_ext);
giri_dir = giri_dir(sorted_idx);

% we'll discard some of the ct images at the start and end
%ct_usable_range = [75:3021];
%ct_dir = ct_dir(ct_usable_range);

% scales in um
ct_voxel_scale = 34.245;
[~,giri_pixel_scale] = dot_target_scale(fullfile(giri_folder,giri_dir(end-2).name));
giri_z_scale = 34.29;

% cut off the last few from the giri_dir
giri_dir = giri_dir(1:end-4);

%% The core is also broken and should be aligned as 2 separate chunks
giri_upper_dir = giri_dir(1:1263);
giri_lower_dir = giri_dir(1324:end);

ct_upper_dir = ct_dir(1:1277);
ct_lower_dir = ct_dir(1317:end);

%% we will use the pyrites to align in 3d. 
% first we need to get some 3d crops around the pyrites in both stacks
vis_3d(giri_upper_dir)
vis_3d(ct_upper_dir)
%% 
[giri_p1_ulc,giri_p1_vol] = crop_3d(giri_upper_dir);
[ct_p1_ulc,ct_p1_vol] = crop_3d(ct_upper_dir);
[giri_p2_ulc,giri_p2_vol] = crop_3d(giri_upper_dir);
[ct_p2_ulc,ct_p2_vol] = crop_3d(ct_upper_dir);
[giri_p3_ulc,giri_p3_vol] = crop_3d(giri_upper_dir);
[ct_p3_ulc,ct_p3_vol] = crop_3d(ct_upper_dir);
[giri_p4_ulc,giri_p4_vol] = crop_3d(giri_upper_dir);
[ct_p4_ulc,ct_p4_vol] = crop_3d(ct_upper_dir);
[giri_p5_ulc,giri_p5_vol] = crop_3d(giri_upper_dir);
[ct_p5_ulc,ct_p5_vol] = crop_3d(ct_upper_dir);

%% get the pyrite outlines
[giri_p1_pc] = outline_point_cloud(giri_p1_vol, giri_pixel_scale, giri_z_scale);
[ct_p1_pc] = outline_point_cloud(ct_p1_vol, ct_voxel_scale, ct_voxel_scale);
[giri_p2_pc] = outline_point_cloud(giri_p2_vol, giri_pixel_scale, giri_z_scale);
[ct_p2_pc] = outline_point_cloud(ct_p2_vol, ct_voxel_scale, ct_voxel_scale);
[giri_p3_pc] = outline_point_cloud(giri_p3_vol, giri_pixel_scale, giri_z_scale);
[ct_p3_pc] = outline_point_cloud(ct_p3_vol, ct_voxel_scale, ct_voxel_scale);
[giri_p4_pc] = outline_point_cloud(giri_p4_vol, giri_pixel_scale, giri_z_scale);
[ct_p4_pc] = outline_point_cloud(ct_p4_vol, ct_voxel_scale, ct_voxel_scale);
[giri_p5_pc] = outline_point_cloud(giri_p5_vol, giri_pixel_scale, giri_z_scale);
[ct_p5_pc] = outline_point_cloud(ct_p5_vol, ct_voxel_scale, ct_voxel_scale);

%% put the pyrite outine coordinates in absolute space and make them into point clouds
giri_top_z = numel(giri_upper_dir) * giri_z_scale;
ct_top_z = numel(ct_upper_dir) * ct_voxel_scale;

giri_p1_pc  = giri_p1_pc + [giri_p1_ulc(2) * giri_pixel_scale,-giri_p1_ulc(1) * giri_pixel_scale, giri_top_z - (giri_p1_ulc(3) * giri_z_scale)];
ct_p1_pc  = ct_p1_pc + [ct_p1_ulc(2) * ct_voxel_scale,-ct_p1_ulc(1) * ct_voxel_scale, ct_top_z - (ct_p1_ulc(3) * ct_voxel_scale)];
giri_p2_pc  = giri_p2_pc + [giri_p2_ulc(2) * giri_pixel_scale,-giri_p2_ulc(1) * giri_pixel_scale, giri_top_z - (giri_p2_ulc(3) * giri_z_scale)];
ct_p2_pc  = ct_p2_pc + [ct_p2_ulc(2) * ct_voxel_scale,-ct_p2_ulc(1) * ct_voxel_scale, ct_top_z - (ct_p2_ulc(3) * ct_voxel_scale)];
giri_p3_pc  = giri_p3_pc + [giri_p3_ulc(2) * giri_pixel_scale,-giri_p3_ulc(1) * giri_pixel_scale, giri_top_z - (giri_p3_ulc(3) * giri_z_scale)];
ct_p3_pc  = ct_p3_pc + [ct_p3_ulc(2) * ct_voxel_scale,-ct_p3_ulc(1) * ct_voxel_scale, ct_top_z - (ct_p3_ulc(3) * ct_voxel_scale)];
giri_p4_pc  = giri_p4_pc + [giri_p4_ulc(2) * giri_pixel_scale,-giri_p4_ulc(1) * giri_pixel_scale, giri_top_z - (giri_p4_ulc(3) * giri_z_scale)];
ct_p4_pc  = ct_p4_pc + [ct_p4_ulc(2) * ct_voxel_scale,-ct_p4_ulc(1) * ct_voxel_scale, ct_top_z - (ct_p4_ulc(3) * ct_voxel_scale)];
giri_p5_pc  = giri_p5_pc + [giri_p5_ulc(2) * giri_pixel_scale,-giri_p5_ulc(1) * giri_pixel_scale, giri_top_z - (giri_p5_ulc(3) * giri_z_scale)];
ct_p5_pc  = ct_p5_pc + [ct_p5_ulc(2) * ct_voxel_scale,-ct_p5_ulc(1) * ct_voxel_scale, ct_top_z - (ct_p5_ulc(3) * ct_voxel_scale)];

giri_pyrite_pc = pointCloud([giri_p1_pc;giri_p2_pc;giri_p3_pc;giri_p4_pc;giri_p5_pc]);
ct_pyrite_pc = pointCloud([ct_p1_pc;ct_p2_pc;ct_p3_pc;ct_p4_pc;ct_p5_pc]);

%% and solve for the transform between the point clouds
threed_transform = pcregistericp(ct_pyrite_pc, giri_pyrite_pc,'Metric','pointToPlane');

eul = rotm2eul(threed_transform.Rotation);

%% load in the full ct stack. It's small, so should fit.

upper_ct_stack = load_ct_stack(ct_upper_dir);

%% transform the ct stack
ct_stack_transformed = upper_ct_stack;

ct_stack_transformed = imwarp(ct_stack_transformed,threed_transform);

%% % after the operations above, the position of the sample in the images
% might not be aligned between giri and ct. So we can just calculate that 
% out. 
%
% looks like we don't know from all of the translations where exactly the z
% alignment (or x-y for that matter). but we can manually estimate it
% roughly and massage it into place. 
%
giri_test_im = im2double(imread(fullfile(giri_folder,giri_dir(300).name)));
% select a part of the ct stack that roughly alignns with that image
ct_test_portion = ct_stack_transformed(:,:,300:360);

% go through the test portion and see what the best x-y shift is on the
% best match to the giri test im
max_corrs = zeros(size(ct_test_portion,3),1);
for i = 1:size(ct_test_portion,3)
    [max_corrs(i),~] = nxc_offset(giri_test_im,imresize(ct_test_portion(:,:,i), ct_voxel_scale/giri_pixel_scale));
end

[~,best_match] = max(max_corrs);

% and use the best correlation image to get the offset
[~,image_offset] = nxc_offset(giri_test_im,imresize(ct_test_portion(:,:,best_match), ct_voxel_scale/giri_pixel_scale));

%%

% get the outlines
giri_slice_outline = outline_point_cloud(giri_test_im,giri_pixel_scale,giri_z_scale,giri_poly);
ct_chunk_outline = outline_point_cloud(ct_test_portion,ct_voxel_scale,ct_voxel_scale);

% fit circles to get centers
[giri_xc,giri_yc,~,~] = circfit(giri_slice_outline(:,1),giri_slice_outline(:,2));
[ct_xc,ct_yc,~,~] = circfit(ct_chunk_outline(:,1),ct_chunk_outline(:,2));

giri_centermean = round([giri_xc,giri_yc]/giri_pixel_scale);
ct_centermean = round([ct_xc,ct_yc]/ct_voxel_scale);

% now we inform the crop of the chunk outline based upon the centers
% get the distances from the center of the giri sample to the edges of the
% image in microns
giri_test_im = im2double(imread(fullfile(giri_folder,giri_dir(300).name)));
giri_dist_above = abs(giri_centermean(2))*giri_pixel_scale;
giri_dist_below = (size(giri_test_im,1) - abs(giri_centermean(2))) * giri_pixel_scale;
giri_dist_left = abs(giri_centermean(1))*giri_pixel_scale;
giri_dist_right = (size(giri_test_im,2) - abs(giri_centermean(1))) * giri_pixel_scale;

% and use those to set indices for grabbing the right portion of the ct
% test area
ct_above_ind = abs(ct_centermean(2)) - round(giri_dist_above/ct_voxel_scale);
ct_beolow_ind = abs(ct_centermean(2)) + round(giri_dist_below/ct_voxel_scale);
ct_left_ind = abs(ct_centermean(1)) - round(giri_dist_left/ct_voxel_scale);
ct_right_ind = abs(ct_centermean(1)) + round(giri_dist_right/ct_voxel_scale);

%% 
ct_test_portion = ct_test_portion([ct_above_ind:ct_beolow_ind],[ct_left_ind:941],:);

%% We'll want to get a rough outline for the sample in the giri images for smart cropping
[giri_poly] = user_polygon(giri_dir);
%% load in the full ct stack. It's small, so should fit.

ct_stack = load_ct_stack(ct_dir);

%% first thing to do is get the ct stack approximately aligned with giri
% starting with judging the initial x-y rotation and translation
% load in two images that should be approximately in the same place in the
% stack
[ct_test_outline] = outline_point_cloud(ct_dir(420), ct_voxel_scale, ct_voxel_scale);
[giri_test_outline] = outline_point_cloud(giri_dir(415), giri_pixel_scale, giri_z_scale, giri_poly);

%% and estimate the initial transformation
%  get the center of a circle fit to either center
[giri_xc,giri_yc,~,~] = circfit(giri_test_outline(:,1),giri_test_outline(:,2));
[ct_xc,ct_yc,~,~] = circfit(ct_test_outline(:,1),ct_test_outline(:,2));
% the translation is the difference
xy_translation = [giri_xc - ct_xc, giri_yc - ct_yc];
ct_translated = ct_test_outline(:,[1,2]) + xy_translation;

%% now we'll estimate the rotation by manually selecting potential corresponding points. 

% manually estimating alpha for now.
alpha = 14*pi/13; 
% perform the rotation on the ct points to check
rot_mat = [cos(-alpha) -sin(-alpha); sin(-alpha) cos(-alpha)];
ct_centered = ct_translated - [giri_xc,giri_yc];
ct_rotated = rot_mat * ct_centered';
ct_2d_registered = ct_rotated' +[giri_xc,giri_yc];
%% Get the points of the external geometries with fiducial markers
[ct_outline_cloud] = outline_point_cloud(ct_dir, ct_voxel_scale, ct_voxel_scale);

[giri_outline_cloud] = outline_point_cloud(giri_dir(1:end), giri_pixel_scale, giri_z_scale, giri_poly);

%% make the initial translations we found for the point cloud from the 2d exercise
ct_centered = ct_outline_cloud - mean(ct_outline_cloud);
alpha = 14*pi/13;
rot_mat3d = [cos(alpha), sin(alpha), 0 ; -sin(alpha), cos(alpha), 0; 0, 0, 1];
ct_rotated = (rot_mat3d * ct_centered')' + mean(ct_outline_cloud);

ct_initial_translate = ct_rotated +[xy_translation(1), xy_translation(2), 0];

%% actually make those into point cloud objects
ct_pt_cloud = pointCloud(ct_initial_translate);

giri_pt_cloud = pointCloud(giri_outline_cloud);

%% denoise and downsample the point clouds
ct_pt_cleaned = pcdenoise(pcdownsample(ct_pt_cloud,'random',0.005));
giri_pt_cleaned = pcdenoise(pcdownsample(giri_pt_cloud,'random',0.005));

%% and register the ct cloud to the giri one

threed_transform = pcregistericp(ct_pt_cleaned, giri_pt_cleaned,'Metric','pointToPlane');


%% 
ct_test_portion = ct_test_portion([ct_above_ind:ct_beolow_ind],[ct_left_ind:ct_right_ind],:);
%%
% run through and use similarity to see where the best match is 
sims = [];
for i = 1:size(ct_test_portion,3)
    sims(i) = ssim(giri_test_im,imresize(ct_test_portion(:,:,i),size(giri_test_im)));
end

%%

ct_reg = imresize(ct_test_portion(:,:,1),size(giri_test_im));

false_color_im = cat(3,giri_test_im,ct_reg,ct_reg);
imshow(false_color_im)

%% 

ct_stack_transformed2 = ct_stack_transformed([ct_above_ind:ct_beolow_ind],[ct_left_ind:ct_right_ind],:);
for i = 1:size(ct_stack_transformed2,3)
    this_im = ct_stack_transformed2(:,:,i);
    imwrite(this_im,['/Users/ryan/Dropbox (Princeton)/giri_ct_alignment_project/siberia_1ab/ct_resample1/' num2str(i) '.tif']);
end