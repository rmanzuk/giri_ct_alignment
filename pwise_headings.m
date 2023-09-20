function [headings] = pwise_headings(X, Y);

    % start with the vectors X->Y
    vecs = Y-X;

    % and convert them to polar coordinates
    [headings,~] = cart2pol(vecs(:,1),vecs(:,2));

    % and get rid of negatives because that'll be easiest
    headings(headings < 0) = (2*pi) + headings(headings < 0);

end
