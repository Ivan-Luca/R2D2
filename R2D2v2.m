function [vec_of_vecs, points,rapports] = R2D2v1(im,kp,area_size,pattern)
if nargin < 1
    area_size = 10;
    patch_size = 3;%3;
end
nkp = size(kp,1);

[L,C] = size(im);
im = rescale(double(im));
mindimage = zeros(L,C,1,'double');

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
[Ix, Iy, Ixy, Iyx, C2, D2] = imcurl(im2,0,0,100);
%[Ix, Iy, Ixy, Iyx, D3] = imcurl(im3,1,1);
%[Ix, Iy, Ixy, Iyx, D4] = imcurl(im4,1,1);
%[Ix, Iy, Ixy, Iyx, D5] = imcurl(im4,1,1);
%[Ix, Iy, Ixy, Iyx, D2d] = imcurl(imrotate(im2,45,'crop'),1,1);
%figure, imshow(imrotate(im2,45))
% [Ix, Iy, Ixy, Iyx, D2] = imcurl(im,1,-1);
% [Ix, Iy, Ixy, Iyx, D3] = imcurl(im,-1,1);
% [Ix, Iy, Ixy, Iyx, D4] = imcurl(im,1,1);


mindimage(:,:,1) = C2/max(C2,[],'all');   %(abs(D2)/max(abs(D2),[],'all')).^(1/3);
mindimage(:,:,2) = D2/max(D2,[],'all');

figure, histogram(mindimage(:,:,1)), ylim([0 10000]);
%figure, histogram(mindimage(:,:,2)), ylim([0 10000]);
% figure, histogram(mindimage(:,:,3));
% figure, histogram(mindimage(:,:,4));

%iterate through kps
totp = 0;
ns = 16;
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
            for j = 1:ns
                for i = 1:ns
                     clip = sqt2(round((j-1)*ys/ns+1):round(j*ys/ns),round((i-1)*xs/ns+1):round(i*xs/ns),:);
                     %%for k = 1:zs
                        %clip2 = clip(:,:,k);
                        %des(k,j,i,:) = permute(hist(clip2(:), 20), [1 3 2]);
                     des(i,j) = mean(clip,'all');
                     
                        %des(i,j,:,:) = logsample(clip,1,round(xs/ns),round(xs/ns),round(ys/ns),8,32);
                    %% end
                     
               end 
           end


            
            sq = uint8(255*(des(:)./max(des(:),[],'all')));
            des2 = uint8(255*(des./max(des(:),[],'all')));
%             for ki = 1:length(pattern)
%                sq(ki) = des2(9+pattern(ki,2), 9+pattern(ki,1)) < des2(9+pattern(ki,4), 9+pattern(ki,3));
%             end

            for ki = 1:length(sq)
                if ki == length(sq)
                    sq2(ki) = sq(ki) < sq(1);
                else
                    sq2(ki) = sq(ki) < sq(ki+1);
                end
%                 if mod(ki,8) == 0
%                     skb(floor(ki/8)) = uint8(0);
%                     for p = 1:8
%                         skb(floor(ki/8)) = skb(floor(ki/8)) + sq2(ki-8+p)*2^(p-1);
%                     end
%                 end            
            end
            sq = uint8(sq2);
%             sq = sq(:);
            
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
            %sq
            vec_of_vecs(totp,:) = sq;%sq.Values;
            %vec_of_vecs2(totp,:) = arr2;
            points(totp,:) = [x,y];
        end
    end
end



vec_of_vecs = binaryFeatures(vec_of_vecs);
%figure, histogram(vec_of_vecs), ylim([0,30000]), title('total descriptors'), hold off;
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


