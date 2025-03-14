function [vec_of_vecs, points,rapports] = testDescriptorRI(im,kp,area_size,patch_size)
if nargin < 1
    area_size = 10;
    patch_size = 3;%3;
end
nkp = size(kp,1);
% hp = floor(patch_size/2);
% ha = floor(area_size)/2;
% %prepare results and variables
% patch = zeros(patch_size);
% patch2 = zeros(patch_size);
% sq = zeros(2*area_size+1,2*area_size+1,8,'double');
% 
% vec_of_vecs = zeros(nkp,((2*area_size+1)^2)*8,'double');
% vec_of_vecs2 = zeros(nkp,640,'double');
% varmin = var(double(im),0,'all');
[L,C] = size(im);
im = rescale(double(im));
mindimage = zeros(L,C,1,'double');
% se = strel('disk',3);


% totp = 0;
% points = zeros(nkp,2);
% for py = 3:L-2
%     for px = 3:C-2
%         patch = im(round(py-hp:py+hp), round(px-hp:px+hp));
%         index = 0;
%         for l = -1:1
%             for c = -1:1
%                 if l ~= 0 || c ~= 0
%                    index = index+1;
%                    patch2 = im(round(py+l-hp:py+l+hp), round(px+c-hp:px+c+hp));
%                    val = (sum(sum((patch-patch2).^2)))/(patch_size^2);
%                    
%                    mindimage(py,px,index) = val;          
% 
%                 end
%             end
%         end
%         vari = var(mindimage(py,px,:));
%         mindimage(py,px,:) = exp(-mindimage(py,px,:)/vari);
%         
%         
% 
%         
% %         pfort = sum(mindimage(py,px,:));
% %         if pfort < 4 
% %             mindimage(px,py,:) = 0;
% %         end
%     end
% end


% mindimage = MIND_descriptor2D(im,1);
% mindimage(:,:,:) = rescale(mindimage(:,:,:),0,1);
%mindimage(:,:,i) = imerode(mindimage(:,:,i),se);
kernel = [ 1/16 1/8 1/16; 1/8 1/4 1/8; 1/16 1/8 1/16];
kernel2 = (1/273)*[1 4 7 4 1; 1 16 26 16 4; 7 26 41 26 7; 1 16 26 16 4; 1 4 7 4 1];
% im2 = conv2(im,kernel2,'same');
% im3 = conv2(im2,kernel2,'same');
% im4 = conv2(im3,kernel2,'same');

im2 = imgaussfilt(im,0.75);
%im3 = imgaussfilt(im2,0.75);
%im4 = imgaussfilt(im3,0.75);
%im5 = imgaussfilt(im4,0.75);

%figure, subplot(2,2,1) , imshow(im), subplot(2,2,2) , imshow(im2),subplot(2,2,3) , imshow(im3),subplot(2,2,4) , imshow(im4);
%[Ix, Iy, Ixy, Iyx, D] = imcurl(im,1,1);
[Ix, Iy, Ixy, Iyx, C2, D2] = imcurl(im2,0,0);
%[Ix, Iy, Ixy, Iyx, D3] = imcurl(im3,1,1);
%[Ix, Iy, Ixy, Iyx, D4] = imcurl(im4,1,1);
%[Ix, Iy, Ixy, Iyx, D5] = imcurl(im4,1,1);
%[Ix, Iy, Ixy, Iyx, D2d] = imcurl(imrotate(im2,45,'crop'),1,1);
%figure, imshow(imrotate(im2,45))
% [Ix, Iy, Ixy, Iyx, D2] = imcurl(im,1,-1);
% [Ix, Iy, Ixy, Iyx, D3] = imcurl(im,-1,1);
% [Ix, Iy, Ixy, Iyx, D4] = imcurl(im,1,1);

% [m1,~,~,~,~,eo1,~] = phasecong3(im,4,6,3,'mult',1.6,'sigmaOnf',0.75,'g', 3, 'k',1);
% 'eo'
% size(eo1)
% size(m1)
% eo1
% m1
% for i = 1:4
%     for j = 1:6
%         figure,imshow(abs(eo1{i,j})), hold off;
%     end
% end
% for c = t+1:C-t-1
%    for l = t+1:L-t-1
%         Dv2(l,c) = D(l,c)/max(var(D(l-t:l+t,c+t:c-t),0,'all'),0);
%         %Dv2(l,c) = D(l,c)/max(mean(D(l-t:l+t,c+t:c-t),'all'),0);
%    end
% end

%mindimage(:,:,1) = atan(Dp-mean(Dp,'all'));
% mindimage(:,:,1) =  D/max(D,[],'all');%uint8(255*rescale(Dp));
% mindimage(:,:,2) =  D/min(D,[],'all');
% A = [ 1 0 1; 0 1 0; 1 0 1];
% B = [0 1 0; 1 0 1; 0 1 0];
% max(A,B)

mindimage(:,:,1) =C2/max(C2,[],'all');   %(abs(D2)/max(abs(D2),[],'all')).^(1/3);
mindimage(:,:,2) =  D2/max(D2,[],'all');
% mindimage(:,:,3) =  Dp4;
% mindimage(:,:,4) =  Dp;
%mindimage(:,:,1) = rescale(Dv);
figure, histogram(mindimage(:,:,1)), ylim([0 10000]);
%figure, histogram(mindimage(:,:,2)), ylim([0 10000]);
% figure, histogram(mindimage(:,:,3));
% figure, histogram(mindimage(:,:,4));

%mindimage(:,:,1) = 1./(1+exp(-2*(Dp-mean(Dp,'all'))));
% mindimage(:,:,2) = 1./(1+exp(2*(Dm-mean(Dm,'all'))));
%[min(mindimage(:,:,1),[],'all'),max(mindimage(:,:,1),[],'all')]
%max(mindimage(:,:,1),[],'all')
% mindimage(:,:,1) = exp(-mindimage(:,:,1));
% max(mindimage(:,:,1),[],'all')
%mindimage(:,:,2) = abs(Dtp);
% max(mindimage(:,:,1),[],'all')
% mindimage(:,:,3) = D3/mean(D3,'all');
% max(mindimage(:,:,1),[],'all')
% mindimage(:,:,4) = D4/mean(D4,'all');
% max(mindimage(:,:,4),[],'all')
% figure,imshow(mindimage(:,:,1)),hold off;
% figure,imshow(D2/max(D2,[],'all')),hold off;
% figure,imshow(mindimage(:,:,2));
% figure,imshow(mindimage(:,:,3));
% figure,imshow(mindimage(:,:,4));
%iterate through kps
totp = 0;
ns = area_size;
o = 6;
sq = zeros(2*area_size+1,2*area_size+1,size(mindimage,3),'double');
rapports =zeros([nkp 1]);
for m = 1:nkp %nkp
    
    y = kp(m,2);
    x = kp(m,1);

    if (x > (area_size+1))  && (x < (C-area_size))
        if (y > (area_size+1)) && (y < (L-area_size))
            totp=totp+1;
            %des = [];
            sqt = mindimage(round(y-area_size:y+area_size),round(x-area_size:x+area_size),:);
            [ys,xs,zs] = size(sqt);
            covM = cov(sqt(:,:,1),sqt(:,:,2));
            f1 = covM(1,1)/(covM(1,1)+covM(2,2));
            f2 = 1-f1;
            [covM(1,1) , var(sqt(:,:,1),0,'all')];
            rapports(m) = f1;
            sqt2 = sqt(:,:,1)*f1 + sqt(:,:,2)*f2;
            %sqt2 = (sqt(:,:,1) + sqt(:,:,2))/2;
            %sqt2 = sqt2*covM(1,2);
%             for j = 1:ns
%                 for i = 1:ns
%                      clip = sqt2(round((j-1)*ys/ns+1):round(j*ys/ns),round((i-1)*xs/ns+1):round(i*xs/ns),:);
%                      size(clip);
%                      %%for k = 1:zs
%                         %clip2 = clip(:,:,k);
%                         %des(k,j,i,:) = permute(hist(clip2(:), 20), [1 3 2]);
%                      des(i,j) = mean(clip,'all');
%                      
%                         %des(i,j,:,:) = logsample(clip,1,round(xs/ns),round(xs/ns),round(ys/ns),8,32);
%                     %% end
%                      
%                end 
%             end
            %imshow(des./max(des(:),[],'all'))
%             pos1 = 0; pos2 = 0; i = 1;
%             p = 6;
%             ws = area_size;
%             while (pos1 + floor(p/2)) < ws
%                 ps = floor(p/2);
%                 clip1 = sqt2((ws-ps-pos1):(ws+ps-pos1),(ws-ps-pos2):(ws+ps-pos2),:);
%                 clip2 = sqt2((ws-ps+pos1):(ws+ps+pos1),(ws-ps+pos2):(ws+ps+pos2),:);
%                 clip3 = sqt2((ws-ps+pos1):(ws+ps+pos1),(ws-ps-pos2):(ws+ps-pos2),:);
%                 clip4 = sqt2((ws-ps-pos1):(ws+ps-pos1),(ws-ps+pos2):(ws+ps+pos2),:);
%                 clip5 = sqt2((ws-ps):(ws+ps),(ws-ps-pos2):(ws+ps-pos2),:);
%                 clip6 = sqt2((ws-ps):(ws+ps),(ws-ps+pos2):(ws+ps+pos2),:);
%                 clip7 = sqt2((ws-ps-pos1):(ws+ps-pos1),(ws-ps):(ws+ps),:);
%                 clip8 = sqt2((ws-ps+pos1):(ws+ps+pos1),(ws-ps):(ws+ps),:);
%                 
%                 clip9 = sqt2((ws-ps-floor(pos1/2)):(ws+ps-floor(pos1/2)),(ws-ps-pos2):(ws+ps-pos2),:);
%                 clip10 = sqt2((ws-ps+floor(pos1/2)):(ws+ps+floor(pos1/2)),(ws-ps+pos2):(ws+ps+pos2),:);
%                 clip11 = sqt2((ws-ps+floor(pos1/2)):(ws+ps+floor(pos1/2)),(ws-ps-pos2):(ws+ps-pos2),:);
%                 clip12 = sqt2((ws-ps-floor(pos1/2)):(ws+ps-floor(pos1/2)),(ws-ps+pos2):(ws+ps+pos2),:);
%                 
%                 clip13 = sqt2((ws-ps-pos1):(ws+ps-pos1),(ws-ps-floor(pos2/2)):(ws+ps-floor(pos2/2)),:);
%                 clip14 = sqt2((ws-ps+pos1):(ws+ps+pos1),(ws-ps+floor(pos2/2)):(ws+ps+floor(pos2/2)),:);
%                 clip15 = sqt2((ws-ps+pos1):(ws+ps+pos1),(ws-ps-floor(pos2/2)):(ws+ps-floor(pos2/2)),:);
%                 clip16 = sqt2((ws-ps-pos1):(ws+ps-pos1),(ws-ps+floor(pos2/2)):(ws+ps+floor(pos2/2)),:);
%                 
%                 
%                 des(1,i) = mean(clip1,'all');
%                 des(2,i) = mean(clip2,'all');
%                 des(3,i) = mean(clip3,'all');
%                 des(4,i) = mean(clip4,'all');
%                 des(5,i) = mean(clip5,'all');
%                 des(6,i) = mean(clip6,'all');
%                 des(7,i) = mean(clip7,'all');
%                 des(8,i) = mean(clip8,'all');
%                 des(9,i) = mean(clip9,'all');
%                 des(10,i) = mean(clip10,'all');
%                 des(11,i) = mean(clip11,'all');
%                 des(12,i) = mean(clip12,'all');
%                 des(13,i) = mean(clip13,'all');
%                 des(14,i) = mean(clip14,'all');
%                 des(15,i) = mean(clip15,'all');
%                 des(16,i) = mean(clip16,'all');
%                 
%                 
%                 %p = p +2;
%                 pos1 = pos1 + ps;
%                 pos2 = pos2 + ps;
%                 i = i +1;
%             end
            
            %des(:,:) =  logsample(sqt2,1,area_size,area_size+1,area_size+1,16,16);
            %des = des/max(des,'all');
%             des = des(:,2:end);
%              [ma, mid] = max(sum(des,2));
%              des = circshift(des,[-(mid-1) 0]);
            %V = mean(sqt,'all');
            %sq = des(:);
            %sq = des(:)./max(des(:),[],'all');
            %sq = uint8(255*(sqt2(:)./max(sqt2(:),[],'all')));%sqt2
            
%              sq = uint8(255*(des(:)./max(des(:),[],'all')));
             [sq,visu] = extractHOGFeatures(sqt2./max(sqt2,[],'all'),'CellSize',[95 95],'BlockSize',[1 1],'NumBins',72,'UseSignedOrientation',true);
             if mod(m,300) == 0
                figure, imshow(sqt2), hold on, plot(visu), hold off;
             end
%             if norm(sq) ~= 0
%                 sq = sq/norm(sq);
%             end
%             for n = 1:*area_size+1
%                 for k = 1:*area_size+1
%                     for p = 1:8
%                    sq(n,k,p) = quantityMIND(sqt(n,k,p))
%                 end
%             end
            
            %figure
            %imshow(sq(:,:,1:3))
            %arr(:) = sq;
            %arr2 = logsample(sq,4,area_size,area_size+1,area_size+1,4,20);
%             if( mod(m,110) == 100)
%             figure
%             imshow(sq(:,:,1:3))
%             figure
%             %imshow(arr2(:,:,1:3))
%             end
            %%arr2 = reshape(arr2,1,[]);
            %arr = rescale(arr,0,1);
%             arr = rescale(arr,0,1);
%             arr2 = rescale(arr2,0,1);
            %sq = histogram(sq,256);
            vec_of_vecs(totp,:) = sq;%sq.Values;
            %vec_of_vecs2(totp,:) = arr2;
            points(totp,:) = [x,y];
        end
    end
end



vec_of_vecs = vec_of_vecs(1:totp,:);
figure, histogram(vec_of_vecs), ylim([0,20000]), title('total descriptors'), hold off;
%vec_of_vecs2 = vec_of_vecs2(1:totp,:);
points = points(1:totp,:);
final_image = (mindimage(:,:,1)*mean(rapports,'all') + mindimage(:,:,2)*(1 - mean(rapports,'all')));
final_image = final_image/max(final_image,[],'all');
figure, imshow( uint8(255*(final_image)));
% figure
% imshow(mindimage(:,:,1:1))
% figure
% imshow(mindimage(:,:,4:4))
% figure
% imshow(mindimage(:,:,6:8))
end


