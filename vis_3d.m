function [] = vis_3d(input_dir)

    %% begin the function

    n_slices = numel(input_dir);
    slide_step = 1/(n_slices-1);
    
    f = figure('Visible','off');
    c = uicontrol(f,'Style','slider');
    c.Min = 1;
    c.Max = n_slices;
    c.Value = 1;
    c.SliderStep = [slide_step slide_step];
    c.Position = [270 10 60 20];
    c.Callback = @selection;
    f.Visible='on';

    function selection(src,event)
        val = c.Value;
        i = round(val);
        this_im = imread(fullfile(input_dir(i).folder, input_dir(i).name));
        imshow(this_im,'InitialMagnification',400)
        title("image #" + (i));
        drawnow
    end

end