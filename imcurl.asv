function [Ixx, Iyy, Ixy, Iyx , C, D] = imcurl( image,dir1, dir2, pad )

G = im2gray(image);

% convert to double
G2 = im2double(G);
G2 = rescale(G2);
% create X and Y filters
horizontal_filter = [-1 2 -1];
vertical_filter = [-1 2 -1]';
% horizontal_filter = 1/9*[1 0 -1; 2 0 -2; 1 0 -1];
% vertical_filter = 1/9*[1 2 1; 0 0 0 ; -1 -2 -1];

if dir1 == 0
horizontal_filter = [-1 1];
end
if dir2 == 0
vertical_filter = [-1 1]';
end
% using imfilter to get our gradient in each direction
filtered_x = pad*conv2(G2, [-1 0 1],'same');
filtered_y = pad*conv2(G2, [-1 0 1]','same');
%refilter
refiltered_xy = pad*conv2(filtered_x, [-1 0 1]','same');
refiltered_yx = pad*conv2(filtered_y, [-1 0 1],'same');

refiltered_xx = pad*conv2(filtered_x, [-1 0 1],'same');
refiltered_yy = pad*conv2(filtered_y, [-1 0 1]','same');

Ixxy = pad*conv2(refiltered_xx, [-1 0 1]','same');
Iyyx = pad*conv2(refiltered_yy, [-1 0 1],'same');
% store the values in our output variables, for clarity
Ix = filtered_x;
Iy = filtered_y;
Ixy = refiltered_xy;
Iyx = refiltered_yx;
Ixx = refiltered_xx;
Iyy = refiltered_yy;
D = abs(Ixx + Iyy);
%D =  (Ix.^2 + Iy.^2)./2;
%C = abs(Ixx.*Iy - Iyy.*Ix);
C = abs(Ixy - Iyx);
%D = C;
% 'Ix et Ixx et C'
% max(Ix,[],'all')
% min(Ix,[],'all')
% max(Ixy,[],'all')
% min(Ixy,[],'all')
% max(C,[],'all')
% min(C,[],'all')
% max(D,[],'all')
% min(D,[],'all')
% figure, imshow(rescale(Ix)),colormap('Turbo'),title('Ix');
% figure, imshow(rescale(Iy)),colormap('Turbo'),title('Iy');
% figure, imshow(rescale(Ixx)),colormap('Turbo'),title('Ixx');
% figure, imshow(rescale(Iyy)),colormap('Turbo'),title('Iyy');
% figure, imshow(rescale(Ixy)),colormap('Turbo'),title('Ixy');
% figure, imshow(rescale(Iyx)),colormap('Turbo'),title('Iyx');
% figure, imshow(rescale(C)),colormap('Turbo'),title('C');
% figure, imshow(rescale(D)),colormap('Turbo'),title('D');

end

