function [max_corr,calculated_offset] = nxc_offset(im1,im2)
    nxc = normxcorr2(edge(im1,'Canny',0.4),edge(im2,'Canny',0.4));
    [max_corr,index_max] = max(abs(nxc(:)));
    [ypeak,xpeak] = ind2sub(size(nxc),index_max(1));
    
    calculated_offset = [(ypeak-size(im1,1)),(xpeak-size(im1,2))];
end