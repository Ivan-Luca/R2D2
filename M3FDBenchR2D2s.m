clear all, close all;
folder = {'descriptor','external_code','ssim'};
for it=1:length(folder)
    p = genpath(folder{it});
    addpath(p);
end
folderTh = 'M3FD/Ir/';
folderVis = 'M3FD/Vis/';
folderTh = 'LGHDDB/Ir/';
folderVis = 'LGHDDB/Vis/';
% folderTh = 'test_images/Ir/';
% folderVis = 'test_images/Vis/';

imagesRgb = dir(folderVis);
imagesLwir = dir(folderTh);
numbs = randi(400,80,1);
step = 5
miss = 0; i = 0; sc =1;
show = 1; 
missTab = int16.empty(200,0);
missAmount = single.empty(200,0);
totError = single.empty(200,0);
timeExec = single.empty(200,0);
pointsStat = zeros([200,4,1],"uint16");
pointsStatMiss = zeros([200,4,1],"uint16");
nbImages = 00;
k = 10;

maxIndex = min(nbImages*step+4+k,length(imagesLwir)+1)
ws = 64;
idexes = uint8.empty(200,0);
while k < maxIndex
    k
    nameRgb = imagesRgb(k).name;
    nameLwir = imagesLwir(k).name;
    im_rgb = im2gray(imread(strcat(folderVis,nameRgb)));
    im_lwirori = im2gray(imread(strcat(folderTh,nameLwir)));
    %im_lwirori = imresize(im_lwirori,[round(size(im_lwirori)./sc)]);
    %im_lwirori = imrotate(im2gray(imread(strcat(folderTh,nameLwir))),270);
    i = i+1;

    im_lwir = im_lwirori;

    lwir_pointsOri = detectMinEigenFeatures(im_lwirori,"MinQuality",0.01, 'ROI',[50 50 size(im_lwir,2)-50 size(im_lwir,1)-50]);
    lwir_pointsControl = lwir_pointsOri.Location;
    
    start = tic;
% 
        rgb_points = detectMinEigenFeatures(im_rgb,"MinQuality",0.001);
        lwir_points = detectMinEigenFeatures(im_lwir,"MinQuality",0.001);
%         [m1,~,~,~,~,eo1,~] = phasecong3(im_rgb,s,6,3,'mult',1.6,'sigmaOnf',0.75,'g', 3, 'k',1);
%         [m2,~,~,~,~,eo2,~] = phasecong3(im_lwir,s,o,3,'mult',1.6,'sigmaOnf',0.75,'g', 3, 'k',1);
%         
%         a=max(m1(:)); b=min(m1(:)); m1=(m1-b)/(a-b);
%         a=max(m2(:)); b=min(m2(:)); m2=(m2-b)/(a-b); 
%         rgb_points = detectFASTFeatures(im_rgb,'MinContrast',0.05);
%         lwir_points = detectFASTFeatures(im_lwir,'MinContrast',0.05);
        
        
        rgb_points = rgb_points.selectStrongest(5000);
        lwir_points = lwir_points.selectStrongest(5000);
   
        [descrgb, pointsrgb, rgbrapports,cv,dv] = R2D2v1(im_rgb,rgb_points.Location,ws,pattern);
        [desclwir, pointslwir,lwirrapports,ct,dt] = R2D2v1(im_lwir,lwir_points.Location,ws,pattern);
        
%             [descrgb2, pointsrgb2, rgbrapports2,cv2,dv2] = R2D2v1(im_rgb,rgb_points.Location,16,pattern);
%             [desclwir2, pointslwir2,lwirrapports2,ct2,dt2] = R2D2v1(im_lwir,lwir_points.Location,16,pattern);
%             
%             [descrgb3, pointsrgb3, rgbrapports3,cv3,dv3] = R2D2v1(im_rgb,rgb_points.Location,32,pattern);
%             [desclwir3, pointslwir3,lwirrapports3,ct3,dt3] = R2D2v1(im_lwir,lwir_points.Location,32,pattern);
%             
%             [indexPairs,matchmetric] = matchFeatures(descrgb2,desclwir2,'MaxRatio',1,'MatchThreshold', 20,'Unique',true,"Method","Exhaustive"); %100 in paper
%             matchedPoints12 = pointsrgb2(indexPairs(:, 1), :);
%             matchedPoints22 = pointslwir2(indexPairs(:, 2), :);
%         
%             [indexPairs,matchmetric] = matchFeatures(descrgb3,desclwir3,'MaxRatio',1,'MatchThreshold', 20,'Unique',true,"Method","Exhaustive"); %100 in paper
%             matchedPoints13 = pointsrgb3(indexPairs(:, 1), :);
%             matchedPoints23 = pointslwir3(indexPairs(:, 2), :);

    [indexPairs,matchmetric] = matchFeatures(descrgb,desclwir,'MaxRatio',1,'MatchThreshold', 99,'Unique',true,"Method","Exhaustive"); %100 in paper
    matchedPoints1 = pointsrgb(indexPairs(:, 1), :);
    matchedPoints2 = pointslwir(indexPairs(:, 2), :);
    
%             matchedPoints1 = [matchedPoints1 ; matchedPoints12; matchedPoints13];
%             matchedPoints2 = [matchedPoints2 ; matchedPoints22 ; matchedPoints23];
%             allmatches = unique([matchedPoints1 matchedPoints2],'rows');
%             matchedPoints1 = allmatches(:,[1 2]);
%             matchedPoints2 = allmatches(:,[3 4]);
    [tform,inliersIndex] = estimateGeometricTransform2D(matchedPoints2, matchedPoints1, 'similarity',"MaxNumTrials",3000,"MaxDistance",5);
    %[tform,inliersIndex] = estimateGeometricTransform2D(matchedPoints2, matchedPoints1, 'projective',"MaxNumTrials",3000,"MaxDistance",5);
    matchedPoints11 = matchedPoints1(inliersIndex, :);
    matchedPoints21 = matchedPoints2(inliersIndex, :);
    matchedBadPoints11 = matchedPoints1(inliersIndex==0, :);
    matchedBadPoints21 = matchedPoints2(inliersIndex==0, :);
        
    
    stop = toc(start);
%     
%     outs = LTFC(im_rgb,im_lwir,matchedPoints1,matchedPoints2);
%     matchedPoints12 = matchedPoints1(outs==0, :);
%     matchedPoints22 = matchedPoints2(outs==0, :);
%     tform2 = fitgeotrans(matchedPoints22, matchedPoints12,'affine');
    
    lwir_pointProj = transformPointsForward(tform,lwir_pointsControl);
    Reproj_error = pdist([0,0; sum(sqrt(lwir_pointProj-lwir_pointsOri.Location).^2)/size(lwir_pointProj,1)]);
    
    if Reproj_error > 50 
        miss = miss +1; i = i- 1;
        missTab(miss) = k;
        missAmount(miss) = Reproj_error;
        pointsStatMiss(miss,:) = [length(pointsrgb),length(pointslwir),length(matchedPoints1),length(matchedPoints11)];  
    else
        totError(i) = Reproj_error;
        timeExec(i) = stop;
        pointsStat(i,:) = [length(pointsrgb),length(pointslwir),length(matchedPoints1),length(matchedPoints11)];
    end
    
     if show ==1
        figure , imshowpair(im_rgb,im_lwirori)
        figure ,imshow(im_rgb);
        figure, imshow(im_lwirori);
        figure, imshowpair(im_lwirori,im_lwir);
        figure, imshow(im_rgb),hold on, plot(rgb_points.Location(:,1),rgb_points.Location(:,2),'o'),hold off;
        figure, imshow(im_lwir),hold on, plot(lwir_points.Location(:,1),lwir_points.Location(:,2),'o');
    end
    
    if    show == 1
        % Show images
        figure; showMatchedFeatures(im_rgb, im_lwir, matchedPoints1, matchedPoints2,'method', 'montage');
        figure; showMatchedFeatures(im_rgb, im_lwir, matchedPoints11, matchedPoints21,'method', 'montage'),title(nameRgb);
        %figure; showMatchedFeatures(im_rgb, im_lwir, matchedPoints12, matchedPoints22,'method', 'montage'),title('LTFC');
        %followOutput = affineOutputView(size(im_rgb),tform,'BoundsStyle','SameasInput');
        %followOutput2 = affineOutputView(size(im_rgb),tform2,'BoundsStyle','SameasInput');
        %warped = imwarp(im_lwir,tform,'OutputView',followOutput);
        Rfixed = imref2d(size(im_rgb));
        warped = imwarp(im_lwir,tform,'OutputView',Rfixed);
        %warped2 = imwarp(im_lwir,tform2,'OutputView',followOutput2);
        figure, imshow(im_lwirori),hold on, plot(lwir_pointsOri.Location(:,1),lwir_pointsOri.Location(:,2),'o');
         plot( lwir_pointProj(:,1), lwir_pointProj(:,2),'x'), title('TRE'), hold off;    
        
        figure , imshowpair(im_rgb,im_lwir), title(strcat(k,' displaced'));
        figure , imshowpair(im_rgb,warped),title(strcat(k,' registered'));
        %figure , imshowpair(im_rgb,warped2),title('registered2');
        figure , imshowpair(warped,im_lwirori),title('TRE');
        idx = uint8(warped > 50);
        figure , imshow(((warped .* idx)+im_rgb)/2),title('fusion');
        descrgb1 = descrgb(indexPairs(:, 1), :);
        desclwir1 = desclwir(indexPairs(:, 2), :);
        descrgb1 = descrgb1(inliersIndex, :);
        desclwir1 = desclwir1(inliersIndex, :);
        figure, subplot(1,2,1), imshow(descrgb1), subplot(1,2,2), imshow(desclwir1), colormap('Turbo');
        [l,c] = size(im_rgb);
        [Y,X] = meshgrid(1:4:l,1:4:c);
        cv = cv(1:4:l,1:4:c);
        dv = dv(1:4:l,1:4:c);
        mv = sqrt(cv.^2 + dv.^2);
        ct = ct(1:4:l,1:4:c);
        dt = dt(1:4:l,1:4:c);
        mt = sqrt(ct.^2 + dt.^2);
        figure, imshow(im_rgb),hold on, quiver(X',Y',cv,dv), axis;
        figure, subplot(1,3,1), imshow(cv),subplot(1,3,2), imshow(ct),subplot(1,3,3), imshow(cv-ct);
        figure, subplot(1,3,1), imshow(dv),subplot(1,3,2), imshow(dt),subplot(1,3,3), imshow(dv-dt);
        figure, subplot(1,3,1), imshow(mv),subplot(1,3,2), imshow(mt),subplot(1,3,3), imshow(mv-mt);
    end
    k = k+step;
end

TRE = sum(totError)/i;
meanTime = sum(timeExec)/i;
length(matchedPoints11)
% = 
% mirgb = MIND_descriptor2D(im_rgb,1);
% milwir = MIND_descriptor2D(im_lwir,1);
% imshow( mirgb(:,:,1:3) );
% imshow( milwir(:,:,1:3) );



