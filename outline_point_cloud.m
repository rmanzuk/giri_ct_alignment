function [point_cloud] = outline_point_cloud(input_ims, im_scale, z_scale,cropping_polygon)
    
    % we may be dealing with either a directory of images or a stack, so
    % check
    is_dir = false;
    num_ims = 0;
    if isa(input_ims,'struct')
        is_dir = true;
        num_ims = numel(input_ims);
    else
        is_dir = false;
        num_ims = size(input_ims,3);
    end
    % set up the empty point cloud
    point_cloud = [];

    % and set up the z maximum
    z_max = num_ims * z_scale;
    
 
    % we need to let the user pick some metrics for binarization
    % if the user put in an outlining polygon, we'll need to black out
    % the stuff outside the polygon
    % start with an image in the middle of the stack
    if is_dir
        initial_im = im2double(imread(fullfile(input_ims(ceil(end/2)).folder,input_ims(ceil(end/2)).name)));
    else
        initial_im = input_ims(:,:,ceil(end/2));
    end
    
    if nargin == 4
        mask = poly2mask(cropping_polygon(:,1),cropping_polygon(:,2), size(initial_im,1), size(initial_im,2));
        initial_im = mask .* initial_im;
    end
    
    figure()
    title('threshold binarizations')
    threshes = [0.05:0.05:0.9];
    for j = 1:numel(threshes)
        subplot(ceil(numel(threshes)/5),5,j)
        imshow(imbinarize(initial_im,threshes(j)))
        title(sprintf('%0.2f', threshes(j)))
        drawnow
    end
    
    % let the user set the binarization threshold
    binarize_thresh = input('Which threshold would you like to use for binarization?\n');
    
    % and make a binary for further testing
    test_binary = imbinarize(initial_im,binarize_thresh);

    % now binarize with that threshold and test island sizes for opening
    island_sizes = [10, 100, 1000, 5000, 10000, 50000, 100000];
    % do conn comps
    CC = bwconncomp(test_binary,4);
    S = regionprops(CC, 'Area');
    L = labelmatrix(CC);
    figure()
    title('island sizes for closing')
    for j = 1:numel(island_sizes)
        subplot(2,4,j)
        % and remove small islands
        imshow(ismember(L, find([S.Area] >= island_sizes(j))));
        title(sprintf('%0.0f', island_sizes(j)))
        drawnow
    end

    % ask the user for the island area they like best
    island_area = input('Which island area would you like to use for opening?\n');
    
    % and update the binary
    test_binary2 = ismember(L, find([S.Area] >= island_area));

    % and invert the binary to do the same thing and close it at several
    % hole sizes
    hole_sizes = [10, 100, 1000, 5000, 10000, 50000, 100000];
    CC = bwconncomp(~test_binary2,4);
    S = regionprops(CC, 'Area');
    L = labelmatrix(CC);
    figure()
    title('hole sizes for closing')
    for j = 1:numel(hole_sizes)
        subplot(2,4,j)
        % and remove small islands
        imshow(~ismember(L, find([S.Area] >= hole_sizes(j))));
        title(sprintf('%0.0f', hole_sizes(j)))
        drawnow
    end
    
    % ask the user for the hole area they like best
    hole_area = input('Which hole area would you like to use for closing?\n');
    
    % and double check the user is happy with the anticipated output
    user_happy = input('Are you happy with the anticipated outputs (Y/N)\n', 's');

    if strcmpi(user_happy,'N')
        disp('Sorry, try again.')
        return
    end

    % now we iterate through all of the images and extract their outlines
    for i = 1:num_ims
        if is_dir
            this_im = im2double(imread(fullfile(input_ims(i).folder,input_ims(i).name)));
        else 
            this_im = input_ims(:,:,i);
        end
        % if the user put in an outlining polygon, we'll need to black out
        % the stuff outside the polygon
        if nargin == 4
            mask = poly2mask(cropping_polygon(:,1),cropping_polygon(:,2), size(this_im,1), size(this_im,2));
            this_im = mask .* this_im;
        end

        
        initial_binary = imbinarize(this_im,binarize_thresh);

        % and we want to open and close that image a bit so we just have the major
        % outline
        % start with conn comp
        CC = bwconncomp(initial_binary,4);
        S = regionprops(CC, 'Area');
        L = labelmatrix(CC);
        % and remove small islands
        binary_2 = ismember(L, find([S.Area] >= island_area));

        % and invert the binary to do the same thing and close it
        CC = bwconncomp(~binary_2,4);
        S = regionprops(CC, 'Area');
        L = labelmatrix(CC);
        % and remove small islands
        binary_3 = ~ismember(L, find([S.Area] >= hole_area));

        % and trace the boundaries
        bound = bwboundaries(binary_3);

        % just in case there is a weird extra boundary, we'll just take the
        % one with the max number of elements
        [~,max_ind] = max(cellfun(@numel,bound));

        % we also need a z coordinate for everything
        z_coord = z_max - (i-1) * z_scale;

        % put the 3 coordinates together in 1 matrix with cols x,y,z
        if numel(bound) > 0
            to_append = [bound{max_ind}(:,2).*im_scale, -bound{max_ind}(:,1).*im_scale, repmat(z_coord,size(bound{max_ind},1),1)];
        else
            to_append = [];
        end

        % and append
        point_cloud = [point_cloud; to_append];
    end
end