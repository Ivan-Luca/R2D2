clear all, close all;
folder = {'descriptor','external_code','ssim'};
for it=1:length(folder)
    p = genpath(folder{it});
    addpath(p);
end
folderTh = 'M3FD/Ir/';
folderVis = 'M3FD/Vis/';

imagesRgb = dir(folderVis);
imagesLwir = dir(folderTh);
numbs = randi(400,80,1);
step = 3
miss = 0; i = 0; sc =1;
show = 0; 
missTab = int16.empty(200,0);
missAmount = single.empty(200,0);
totError = single.empty(200,0);
timeExec = single.empty(200,0);
pointsStat = zeros([200,4,1],"uint16");
pointsStatMiss = zeros([200,4,1],"uint16");
nbImages = 100;
k = 3;
maxIndex = min(nbImages*step+4+k,length(imagesLwir)+1)
idexes = uint8.empty(200,0);
while k < maxIndex
    k
    nameRgb = imagesRgb(k).name
    nameLwir = imagesLwir(k).name
    im_rgb = im2gray(imread(strcat(folderVis,nameRgb)));
    im_lwirori = im2gray(imread(strcat(folderTh,nameLwir)));
    %im_lwirori = imresize(im_lwirori,[round(size(im_lwirori)./sc)]);

    i = i+1;

    im_lwir = im_lwirori;

    lwir_pointsOri = detectMinEigenFeatures(im_lwirori,"MinQuality",0.01, 'ROI',[50 50 size(im_lwir,2)-50 size(im_lwir,1)-50]);
    lwir_pointsControl = lwir_pointsOri.Location;
    
    start = tic;

   [des_m1,des_m2] = RIFT_no_rotation_invariance(im_rgb,im_lwir,4,6,96);

   
   
    disp('nearest matching')
    % nearest matching
    [indexPairs,matchmetric] = matchFeatures(des_m1.des,des_m2.des,'MaxRatio',1,'MatchThreshold', 100);
    matchedPoints1 = des_m1.kps(indexPairs(:, 1), :);
    matchedPoints2 = des_m2.kps(indexPairs(:, 2), :);
    [matchedPoints2,IA]=unique(matchedPoints2,'rows');
    matchedPoints1=matchedPoints1(IA,:);
    
    disp('outlier removal')
    %outlier removal
    H=FSC(matchedPoints1,matchedPoints2,'affine',2);
    Y_=H*[matchedPoints1';ones(1,size(matchedPoints1,1))];
    Y_(1,:)=Y_(1,:)./Y_(3,:);
    Y_(2,:)=Y_(2,:)./Y_(3,:);
    E=sqrt(sum((Y_(1:2,:)-matchedPoints2').^2));
    inliersIndex=E<3;
    cleanedPoints1 = matchedPoints1(inliersIndex, :);
    cleanedPoints2 = matchedPoints2(inliersIndex, :);
        
    
    stop = toc(start);

     tform = estimateGeometricTransform2D(cleanedPoints1,cleanedPoints2,'projective');
    
     lwir_pointProj = transformPointsForward(tform,lwir_pointsControl);
     Reproj_error = pdist([0,0; sum(sqrt(lwir_pointProj-lwir_pointsOri.Location).^2)/size(lwir_pointProj,1)]);
    
    if Reproj_error > 50 
        miss = miss +1; i = i- 1;
        missTab(miss) = k;
        missAmount(miss) = Reproj_error;
        pointsStatMiss(miss,:) = [length(des_m1.kps),length(des_m2.kps),length(matchedPoints1),length(cleanedPoints1)];  
    else
        totError(i) = Reproj_error;
        timeExec(i) = stop;
        pointsStat(i,:) = [length(des_m1.kps),length(des_m2.kps),length(matchedPoints1),length(cleanedPoints1)];
    end
    
     if show ==1
        figure , imshowpair(im_rgb,im_lwirori)
        figure ,imshow(im_rgb);
        figure, imshow(im_lwirori);
        figure, imshowpair(im_lwirori,im_lwir);
        figure, imshow(im_rgb),hold on, plot(des_m1.kps(:,1),des_m1.kps(:,2),'o'),hold off;
        figure, imshow(im_lwir),hold on, plot(des_m2.kps(:,1),des_m2.kps(:,2),'o');
    end
    
    if   Reproj_error > 50 || show == 1
        % Show images
        figure; showMatchedFeatures(im_rgb, im_lwir, matchedPoints1, matchedPoints2,'method', 'montage');
        figure; showMatchedFeatures(im_rgb, im_lwir, cleanedPoints1, cleanedPoints2,'method', 'montage'),title(nameRgb);
        %figure; showMatchedFeatures(im_rgb, im_lwir, matchedPoints12, matchedPoints22,'method', 'montage'),title('LTFC');
        %followOutput = affineOutputView(size(im_rgb),tform,'BoundsStyle','SameasInput');
        %followOutput2 = affineOutputView(size(im_rgb),tform2,'BoundsStyle','SameasInput');
        %warped = imwarp(im_lwir,tform,'OutputView',followOutput);
        Rfixed = imref2d(size(im_rgb));
        warped = imwarp(im_lwir,tform,'OutputView',Rfixed);
        %warped2 = imwarp(im_lwir,tform2,'OutputView',followOutput2);
        figure, imshow(im_lwirori),hold on, plot(lwir_pointsOri.Location(:,1),lwir_pointsOri.Location(:,2),'o');
         plot( lwir_pointProj(:,1), lwir_pointProj(:,2),'x'), title('TRE'), hold off;    
        
        figure , imshowpair(im_rgb,im_lwir), title('displaced');
        figure , imshowpair(im_rgb,warped),title('registered');
        %figure , imshowpair(im_rgb,warped2),title('registered2');
        figure , imshowpair(warped,im_lwirori),title('TRE');
        idx = uint8(warped > 50);
        figure , imshow(((warped .* idx)+im_rgb)/2),title('fusion');
    
    end
    k = k+step;
end

TRE = sum(totError)/i;
meanTime = sum(timeExec)/i;
% = 
% mirgb = MIND_descriptor2D(im_rgb,1);
% milwir = MIND_descriptor2D(im_lwir,1);
% imshow( mirgb(:,:,1:3) );
% imshow( milwir(:,:,1:3) );


