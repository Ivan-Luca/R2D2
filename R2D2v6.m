function [vec_of_vecs, points,rapports,m1,m2] = R2D2v6(im,kp,area_size,patch_size)
if nargin < 1
    area_size = 96;
    patch_size = 3;%3;
end
base = 20;
nkp = size(kp,1);

[L,C] = size(im);
im = rescale(double(im));


% kernel = [ 1/16 1/8 1/16; 1/8 1/4 1/8; 1/16 1/8 1/16];
% kernel2 = (1/273)*[1 4 7 4 1; 1 16 26 16 4; 7 26 41 26 7; 1 16 26 16 4; 1 4 7 4 1];


    im2 = imgaussfilt(im,0.75);


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
des = zeros(ns,ns);
vec_of_vecs = zeros(nkp,256);
points = zeros(nkp,2);
rapports =zeros([nkp 1]);
for m = 1:nkp %nkp
    
    y = kp(m,2);
    x = kp(m,1);

    if (x > (area_size+1))  && (x < (C-area_size))
        if (y > (area_size+1)) && (y < (L-area_size))
            totp=totp+1;
            sqt = CD(round(y-area_size:y+area_size),round(x-area_size:x+area_size),:);
            [ys,xs,zs] = size(sqt);
            %sqt2 = CD(round(y-area_size:y+area_size),round(x-area_size:x+area_size),2);
            for j = 1:ns
                for i = 1:ns
                     clip = sqt(round((j-1)*ys/ns+1):round(j*ys/ns),round((i-1)*xs/ns+1):round(i*xs/ns),:);
                     
                     des(i,j) = atan(mean(clip(:,:,2),'all')/mean(clip(:,:,1),'all'));
                     %des(i,j,:) = [mean(cumsum(clip,1),'all'), mean(cumsum(clip,2),'all')] ;
                     %des(i,j) = max(clip,[],'all');
               end 
            end
            
            %sq = des(:)./max(des(:),[],'all');
            vec_of_vecs(totp,:) = des(:)./max(des(:),[],'all');% sq;
            %vec_of_vecs(totp,:) = uint8(255*(des(:)./max(des(:),[],'all')));
            points(totp,:) = [x,y];
        end
    end
end

vec_of_vecs = vec_of_vecs(1:totp,:);
%level = graythresh(vec_of_vecs);
%vec_of_vecs =  imbinarize(vec_of_vecs);
histogram(vec_of_vecs)

end


