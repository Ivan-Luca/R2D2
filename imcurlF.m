function [Ix, Iy, Ixy, Iyx , C, D] = imcurlF( image,dir1, dir2, pad )

G = im2gray(image);

% convert to double
G2 = im2double(G);
G2 = rescale(G2);
% create X and Y filters
horizontal_filter = [[0 0 0]; [-1 0 1];[0 0 0]];
vertical_filter = [[0 -1 0]; [0 0 0];[0 1 0]];

fusion_filter = [[1 0 -1]; [0 0 0];[-1 0 1]];

if dir1 == 0
horizontal_filter = [-1 1];
end
if dir2 == 0
vertical_filter = [-1 1]';
end
% using imfilter to get our gradient in each direction
filtered_x = pad*conv2(G2, horizontal_filter,'same');
filtered_y = pad*conv2(G2, vertical_filter,'same');
%refilter
refiltered_xy = pad*conv2(filtered_x, vertical_filter,'same');
refiltered_yx = pad*conv2(filtered_y, horizontal_filter,'same');

refiltered_xx = pad*conv2(filtered_x, horizontal_filter,'same');
refiltered_yy = pad*conv2(filtered_y, vertical_filter,'same');

fusion1 = pad*conv2(G2, fusion_filter,'same');
fusion2 = pad*conv2(G2, fusion_filter,'same');

ffilterh = fft2(horizontal_filter)
ffilterv = fft2(vertical_filter)
fim = fft2(G2)

% store the values in our output variables, for clarity
Ix = filtered_x;
Iy = filtered_y;
Ixy = refiltered_xy;
Iyx = refiltered_yx;
Ixx = refiltered_xx;
Iyy = refiltered_yy;

Iyx = ifft2((fim .* ffilterh) .* ffilterv);
Ixy = ifft2((fim .* ffilterv) .* ffilterh);
D = abs(Ixx + Iyy);
C = abs(Iyx-Ixy);


end

