clear all, close all;
folder = {'descriptor','external_code','ssim'};
for it=1:length(folder)
    p = genpath(folder{it});
    addpath(p);
end
folderTh = 'LGHDDB/Ir/';
folderVis = 'LGHDDB/Vis/';

% folderTh = 'test_images/Ir/';
% folderVis = 'test_images/Vis/';

imagesRgb = dir(folderVis);
imagesLwir = dir(folderTh);
numbs = randi(400,80,1);
step = 1
miss = 0; i = 0; sc =1;
show = 0; 
missTab = int16.empty(200,0);
Amount = single.empty(200,0);
totError = single.empty(200,0);
timeExec = single.empty(200,0);
pointsStat = zeros([200,4,1],"uint16");
pointsStatMiss = zeros([200,4,1],"uint16");
nbImages = 100;
k =3;

maxIndex = min(nbImages*step+4+k,length(imagesLwir)+1)
ws = 47;
idexes = uint8.empty(200,0);
id = 0;
while k < maxIndex
    k
    nameRgb = imagesRgb(k).name;
    nameLwir = imagesLwir(k).name;
    im_rgb = im2gray(imread(strcat(folderVis,nameRgb)));
    im_lwirori = im2gray(imread(strcat(folderTh,nameLwir)));

    im_lwir = im_lwirori;

%     rgb_points = detectMinEigenFeatures(im_rgb,"MinQuality",0.000001);
%     lwir_points = detectMinEigenFeatures(im_lwir,"MinQuality",0.000001);
%     rgb_points = detectHarrisFeatures(im_rgb,"MinQuality",0.000001);
%     lwir_points = detectHarrisFeatures(im_lwir,"MinQuality",0.000001);
%     rgb_points = detectHessianFeatures('Shi',im_rgb,0.000001);
%     lwir_points = detectHessianFeatures('Shi',im_lwir,0.000001);
    rgb_points = detectSURFFeatures(im_rgb,'MetricThreshold',0.01);
    lwir_points = detectSURFFeatures(im_lwir,'MetricThreshold',0.01);
   % rgb_points = detectSIFTFeatures(im_rgb,'ContrastThreshold',0.01); %%
   % version 2022
    %lwir_points = detectSIFTFeatures(im_lwir,'ContrastThreshold',0.01);
    
    
    
    rgb_points = rgb_points.selectStrongest(5000);
    lwir_points = lwir_points.selectStrongest(5000);
    repet = 0;
    mat = pdist2(rgb_points.Location,lwir_points.Location);
    repet = (mat <= 2);
    corrs = sum(repet);
    tot = sum(corrs > 0);

    id = id+1;
    Amount(id) = tot;
    k = k+step;
    
    if show ==1
        figure , imshowpair(im_rgb,im_lwirori)
        figure ,imshow(im_rgb);
        figure, imshow(im_lwirori);
        figure, imshowpair(im_lwirori,im_lwir);
        figure, imshow(im_rgb),hold on, plot(rgb_points.Location(:,1),rgb_points.Location(:,2),'o'),hold off;
        figure, imshow(im_lwir),hold on, plot(lwir_points.Location(:,1),lwir_points.Location(:,2),'o');
    end
end


meanRep = sum(Amount,'all')/id




