function [Ix, Iy, Ixy, Iyx , C, D] = imcurlT( image,dir1, dir2, pad )

G = im2gray(image);

% convert to double
G2 = im2double(G);
G2 = rescale(G2);
% create X and Y filters
horizontal_filter = [[0 0 0]; [-1 0 1];[0 0 0]];
vertical_filter = [[0 -1 0]; [0 0 0];[0 1 0]];
% horizontal_filter = 1/9*[1 0 -1; 2 0 -2; 1 0 -1];
% vertical_filter = 1/9*[1 2 1; 0 0 0 ; -1 -2 -1];

if dir1 == 0
horizontal_filter = [0 0 1; -1 0 1;-1 0 0];
end
if dir2 == 0
vertical_filter = [-1 -1 0;0 0 0;0 1 1];
end
% using imfilter to get our gradient in each direction
filtered_x = pad*conv2(G2, horizontal_filter,'same');
filtered_y = pad*conv2(G2, vertical_filter,'same');
%refilter
refiltered_xy = pad*conv2(filtered_x, vertical_filter,'same');
refiltered_yx = pad*conv2(filtered_y, horizontal_filter,'same');

refiltered_xx = pad*conv2(filtered_x, horizontal_filter,'same');
refiltered_yy = pad*conv2(filtered_y, vertical_filter,'same');

% store the values in our output variables, for clarity
Ix = filtered_x;
Iy = filtered_y;
Ixy = refiltered_xy;
Iyx = refiltered_yx;
Ixx = refiltered_xx;
Iyy = refiltered_yy;
D = abs(Ixx + Iyy);
C = abs(refiltered_xy-refiltered_yx);


end

