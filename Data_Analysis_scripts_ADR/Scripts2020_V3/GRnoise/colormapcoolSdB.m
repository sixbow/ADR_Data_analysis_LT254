function cmap_out = colormapcoolSdB(points)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
cmap = colormap(cool(points));
cmap_out(:,:)=cmap(:,:).*0.8;

end

