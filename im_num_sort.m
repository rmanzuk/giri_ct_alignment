function [sorted_inds] = im_num_sort(im_folder, im_ext)
    
    % first make the input folder a directory
    im_dir = dir(fullfile(im_folder,im_ext));

    % and get just the names 
    im_names = {im_dir(:).name};

    % split the strings of the name with the underscore
    split_names = cellfun(@(x) strsplit(x, '_'), im_names, 'UniformOutput', false);
    split_names = vertcat(split_names{:}); % To remove nesting of cell array

    % remove the extension from the names and make into double array
    [~, im_nums, ~] = cellfun(@fileparts, split_names(:,2), 'UniformOutput', false);
    im_nums = str2double(im_nums);

    % and finally sort
    [~, sorted_inds] = sort(im_nums);
end