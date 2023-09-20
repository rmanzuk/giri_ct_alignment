function [ct_stack] = load_ct_stack(input_dir)
    
    % First we'll get the size information of the stack so we can
    % preallocate
    test_im = im2double(imread(fullfile(input_dir(1).folder,input_dir(1).name)));

    ct_stack = zeros(size(test_im,1), size(test_im,2), numel(input_dir));
    
    % then just iterate through, read the images, and place them
    for i = 1:numel(input_dir)
        ct_stack(:,:,i) = im2double(imread(fullfile(input_dir(i).folder,input_dir(i).name))); 
    end
end