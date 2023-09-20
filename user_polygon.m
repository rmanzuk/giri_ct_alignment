function [outer_poly] = user_polygon(input_dir)
    
    % set up the empty outer polygon coordinates
    outer_poly = [];
    
    % we'll want to sample a few images from the stack so we know the
    % polygon works throughout the stack
    sample_inds = randsample(numel(input_dir),5);

    % read in all of the images and take their mean
    comp_im = [];
    for i = 1:numel(sample_inds)
        comp_im(:,:,i) = im2double(imread(fullfile(input_dir(sample_inds(i)).folder,input_dir(sample_inds(i)).name)));
    end
    
    mean_im = mean(comp_im,3);
    
    % and then we'll show that mean images and ask the user to input the
    % polygon, displaying the clicking as we go.
    disp('please click the outline of a rough polygon that contains the sample.')
    imshow(mean_im)
    ax = gca;
    ax.Toolbar.Visible = 'off';
    hold on

    title('Polygon tracing'); 
    y_coords = [];
    x_coords = [];
    n = 0;
    while true
        [x_i,y_i] = ginput(1);
        if isempty(x_i) ; break; end
        n = n+1;
        x_coords(n) = x_i(1);
        y_coords(n) = y_i(1);
        plot(x_coords,y_coords,'r','LineWidth',1)
        drawnow
    end

    % and set the polygon for output. x_coords are the columns, y_coords
    % are the rows
    outer_poly = [round(x_coords)', round(y_coords)'];
end