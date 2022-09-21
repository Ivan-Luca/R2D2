
clear all, close all;
folder = {'descriptor','external_code','ssim'};
for it=1:length(folder)
    p = genpath(folder{it});
    addpath(p);
end
descriptor = 'LSSD';  %  'LGHD'|'EHD'|'PCEHD' | 'ELGHD' | 'HOG' | 'LSSD' | 'EMD'
% use classic = 1 to use descriptors above (not MIND)
%folder = 'test_images/'
folderTh = 'M3FD2/Ir/';
folderVis = 'M3FD2/Vis/';
% folderTh = 'test_images/Ir/';
% folderVis = 'test_images/Vis/';
imagesRgb = dir(folderVis);
imagesLwir = dir(folderTh);
numbs = randi(400,80,1);
step = 20
miss = 0; i = 0; show = 1; classic = 0, sharpen = 0, sc =1;
missTab = int16.empty(200,0);
missAmount = single.empty(200,0);
totError = single.empty(200,0);
timeExec = single.empty(200,0);
pointsStat = zeros([200,4,1],"uint16");
pointsStatMiss = zeros([200,4,1],"uint16");
nbImages = 0;
k = 2201;
maxIndex = min(nbImages*step+4+k,length(imagesLwir)+1)
ws = 47;
idexes = uint8.empty(200,0);

k
nameRgb = imagesRgb(k).name
nameLwir = imagesLwir(k).name
im_rgb = im2gray(imread(strcat(folderVis,nameRgb)));
im_lwirori = im2gray(imread(strcat(folderTh,nameLwir)));
im_lwirori = imresize(im_lwirori,[round(size(im_lwirori)./sc)]);
