function [Ixx, Iyy, Ixy, Iyx , C, D] = imcurlop( image,pad )

G = rescale(double(im2gray(image)));

% store the values in our output variables, for clarity
[ Ix,Iy] = imgradientxy(G,'intermediate');
[Ixx,Ixy] = imgradientxy(pad*Ix,'intermediate');
[Iyx,Iyy] = imgradientxy(pad*Iy,'intermediate');
D = abs(pad*(Ixx + Iyy));
C = abs(pad*(Ixy-Iyx));


end

