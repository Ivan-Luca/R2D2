function [vec_of_vecs, points,rapports,m1,m2] = R2D2vBi(im,kp,area_size,patch_size)

nkp = size(kp,1);

[L,C] = size(im);
im = rescale(double(im));


% kernel = [ 1/16 1/8 1/16; 1/8 1/4 1/8; 1/16 1/8 1/16];
% kernel2 = (1/273)*[1 4 7 4 1; 1 16 26 16 4; 7 26 41 26 7; 1 16 26 16 4; 1 4 7 4 1];


    im2 = imgaussfilt(im,1.6);


 [Ixx, Iyy, Ixy, Iyx, C2, D2] = imcurl(im,1,1,100);
% [Ix, Iy, Ixy, Iyx, C2, D2] = imcurlFr(im2);
%C2 = abs(Ixx .* Iyy - Ixy .* Iyx);
%D2 = abs( Ixx + Iyy ) ;
m1 = C2./max(C2,[],'all');
m2 = D2./max(D2,[],'all');
CD = cat(3,m1,m2);
% figure,imshow(CD(:,:,1)),colormap('turbo'),title('CCD');
% figure,imshow(CD(:,:,2)),colormap('turbo'),title('DCD');
% figure,imshow(CD(:,:,1)*0.3 + CD(:,:,2)*0.7),colormap('turbo');
%iterate through kps
totp = 0;
ns = 16;
ps = floor((2*area_size)/ns);
kernel = genGausKernel(ps,(ps-1)/6);
des = zeros(ns,ns,'uint8');
vec_of_vecs = zeros(nkp,ns*ns,'uint8');
points = zeros(nkp,2);
rapports =zeros([nkp 1]);

rng('default');
rng(1); % reapetable random
pairs = randi(16,256,4);
for m = 1:nkp %nkp
    
    y = kp(m,2);
    x = kp(m,1);

    if (x > (area_size+1))  && (x < (C-area_size))
        if (y > (area_size+1)) && (y < (L-area_size))
            totp=totp+1;
            sqt = CD(round(y-area_size:y+area_size),round(x-area_size:x+area_size),:);

            [ys,xs,zs] = size(sqt);
            covM = cov(sqt(:,:,1),sqt(:,:,2));
            f1 = covM(1,1)/(covM(1,1)+covM(2,2));
            f2 = 1-f1;
            rapports(m) = f1;
            sqt2 = sqt(:,:,1)*f1 + sqt(:,:,2)*f2;
            %sqt2 = CD(round(y-area_size:y+area_size),round(x-area_size:x+area_size),2);
            %size(sqt2)
            %patch = sqt2(round(((ns/2)-1)*ys/ns+1):round(((ns/2)-1)*ys/ns+(ns/2)),round(((ns/2)-1)*xs/ns+1):round(((ns/2)-1)*xs/ns+(ns/2)),:);
            for j = 1:ns
                for i = 1:ns
                    q = pairs((j-1)*ns+i,1);s = pairs((j-1)*ns+i,2);d = pairs((j-1)*ns+i,3);f = pairs((j-1)*ns+i,4);
                    patch = sqt2(floor((q-1)*(ys/ns))+1:floor((q-1)*(ys/ns))+ps,floor((s-1)*(xs/ns))+1:floor((s-1)*(xs/ns))+ps,:);
                    clip = sqt2(floor((d-1)*(ys/ns))+1:floor((d-1)*(ys/ns))+ps,floor((f-1)*(xs/ns))+1:floor((f-1)*(xs/ns))+ps,:);
                    %m = mean((clip-patch),[1 2]);
                    
                    %m = sum(patch(:,:).*kernel,'all') - sum(clip(:,:).*kernel,'all');
                    m = max(patch(:,:),[],'all') - max(clip(:,:),[],'all');
                    if  m(:,:) <= 0
                        des(i,j) = 0;
                    else
                        des(i,j) = 1;
                    end
%                    m = exp(-3*(sum(sum(patch(:,:).*kernel) - sum(clip(:,:).*kernel))));
%                     if  m(:,:) <= 0.8
%                         des(i,j) = 0;
%                     else
%                         des(i,j) = 1;
%                    end
%                     m1 = sum(sum(patch(:,:,1).*kernel) - sum(clip(:,:,1).*kernel));
%                     m2 = sum(sum(patch(:,:,2).*kernel) - sum(clip(:,:,2).*kernel));
%                     m= cat(3,m1,m2);
%                     for mit = 1:size(m,3)
%                          if  m(:,:,mit) <= 0
%                             des(i,j,mit) = 0;
%                         else
%                             des(i,j,mit) = 1;
%                          end
%                     end
%                      
                     %maxi, des(i,j) = max(hist(clip, 1:4),'all');
                     %des(i,j,:) = [mean(cumsum(clip,1),'all'), mean(cumsum(clip,2),'all')] ;
                     %des(i,j) = max(clip,[],'all');
               end 
            end

            %sq = des(:)./max(des(:),[],'all');
            %sq = compress(des(:));
            %sq = dec2bin(compress(des(:)));
            vec_of_vecs(totp,:) =  des(:);
            %vec_of_vecs(totp,:) = uint8(255*(des(:)./max(des(:),[],'all')));
            %vec_of_vecs(totp,:) = uint16(65535*(des(:)./max(des(:),[],'all')));
            points(totp,:) = [x,y];
        end
    end
end

%vec_of_vecs = vec_of_vecs(1:totp,:);

vec_of_vecs = binaryFeatures(vec_of_vecs(1:totp,:));
%level = graythresh(vec_of_vecs);
%vec_of_vecs =  imbinarize(vec_of_vecs);
%histogram(vec_of_vecs)

end

function bin = compress(raw)
bin = zeros(512/8,1,'uint8');
for i=1:8:length(raw)  
    bin(ceil(i/8)) = raw(i) + 2*raw(i+1) + 4*raw(i+2) + 8*raw(i+3) + 16*raw(i+4) + 32*raw(i+5) + 64*raw(i+6) + 128*raw(i+7);
end
end

function gFilter = genGausKernel(size,sigma)
gFilter = zeros(size,size);
c = ceil(size/2);
for col = 1 : size
  for row = 1 : size
    gFilter(row, col) = exp(-((col-c)^2+(row-c)^2)/(2*sigma)^2);
  end
end
end
