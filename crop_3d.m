function [upper_left_coord,cropped_vol] = crop_3d(input_dir)
    % propmt the user to id the start and end image for their deisred
    % volume
    top_im_num = input('Which image number should be the top?\n');
    bot_im_num = input('Which image number should be the bottom?\n');

    % and then we'll want to display 1 image for the user to click the row
    % and column coordinates 
    mean_im_num = round((top_im_num +bot_im_num)/2);

    disp('Please click the upper left and lower right corners of the crop.')
    figure()
    imshow(imread(fullfile(input_dir(mean_im_num).folder, input_dir(mean_im_num).name)))

    [col_coords,row_coords] = ginput(2);

    row_coords = round(row_coords);
    col_coords = round(col_coords);

    upper_left_coord = [row_coords(1),col_coords(1),top_im_num];

    cropped_vol = zeros(diff(row_coords)+1, diff(col_coords)+1, bot_im_num-top_im_num+1);
    
    for jj = 1:size(cropped_vol,3)
        for_crop = im2double(imread(fullfile(input_dir(jj).folder, input_dir(top_im_num + jj).name)));
        cropped_vol(:,:,jj) = for_crop(row_coords(1):row_coords(2),col_coords(1):col_coords(2));
    end
    
end