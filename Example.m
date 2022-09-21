% Polar/Rectangular Conversion 
% V0.1 16 Dec 2007 (Created) Prakash Manandhar, pmanandhar@umassd.edu
im = rgb2gray(imread('TestIm.PNG'));
im = double(im)/255.0;
figure, imshow(im),, title('im');
imP = ImToPolar(im, 0.6, 1, 40, 200);
figure, imshow(imP), title('imp');

imR = PolarToIm(imP, 0.6, 1, 250, 250);
figure, imshow(imR), title('imr');

rMin = 0.25; rMax = 0.8;

im2 = imread('TestIm2.jpg');
figure, imshow(im2), title('im2');
imR2 = PolarToIm(im2, rMin, rMax, 300, 300);
figure; imshow(imR2, [0 255]), title('im2R');