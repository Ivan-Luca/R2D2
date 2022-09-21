function [vec_of_vecs, points,rapports,C2,D2] = R2D2v1(im,kp,area_size,patch_size)

nkp = size(kp,1);

[L,C] = size(im);
im = rescale(double(im));


% kernel = [ 1/16 1/8 1/16; 1/8 1/4 1/8; 1/16 1/8 1/16];
% kernel2 = (1/273)*[1 4 7 4 1; 1 16 26 16 4; 7 26 41 26 7; 1 16 26 16 4; 1 4 7 4 1];


    im2 = imgaussfilt(im,0.75);


 [Ixx, Iyy, Ixy, Iyx, C2, D2] = imcurl(im2,1,1,100);
% [Ix, Iy, Ixy, Iyx, C2, D2] = imcurlFr(im2);
%C2 = abs(Ixx .* Iyy - Ixy .* Iyx);
%D2 = abs( Ixx + Iyy ) ;
m1 = C2./max(C2,[],'all');
m2 = D2./max(D2,[],'all');
%CD = cat(3,m1,m2);
% figure,imshow(CD(:,:,1)),colormap('turbo'),title('CCD');
% figure,imshow(CD(:,:,2)),colormap('turbo'),title('DCD');
% figure,imshow(CD(:,:,1)*0.3 + CD(:,:,2)*0.7),colormap('turbo');
%iterate through kps
totp = 0;
ns = 16;
des = zeros(ns,ns,1);
vec_of_vecs = zeros(nkp,ns*ns*1,'uint8');
points = zeros(nkp,2);
rapports =zeros([nkp 1]);

 angles = atan2(m1,m2);

 %figure, histogram(angles,[-pi/4,pi/4])
for m = 1:nkp %nkp
    
    y = kp(m,2);
    x = kp(m,1);

    if (x > (area_size+1))  && (x < (C-area_size))
        if (y > (area_size+1)) && (y < (L-area_size))
            totp=totp+1;
%             sqt = CD(round(y-area_size:y+area_size),round(x-area_size:x+area_size),:);
            sqt3 = angles(round(y-area_size:y+area_size),round(x-area_size:x+area_size),:);
             [ys,xs,zs] = size(sqt3);
%             covM = cov(sqt(:,:,1),sqt(:,:,2));
%             f1 = covM(1,1)/(covM(1,1)+covM(2,2));
%             f2 = 1-f1;
%             rapports(m) = f1;
%             sqt2 = sqt(:,:,1)*f1 + sqt(:,:,2)*f2;
            %sqt2 = CD(round(y-area_size:y+area_size),round(x-area_size:x+area_size),2);
            
            for j = 1:ns
                for i = 1:ns
                     clip = sqt3(round((j-1)*ys/ns+1):round(j*ys/ns),round((i-1)*xs/ns+1):round(i*xs/ns),:);
                     %clip2 = sqt2(round((j-1)*ys/ns+1):round(j*ys/ns),round((i-1)*xs/ns+1):round(i*xs/ns),:);

                     zerg = sum(clip==0,'all');

                     norma = (2*area_size/ns)^2;
                     %des(i,j,:) = (norma-zerg);
                     if j == 1 && i == 1
                        des(i,j,:) = (norma-zerg);
                     elseif i ==1
                        des(i,j-1,:) = (norma-zerg);
                     else
                         des(i,j,:) =des(i-1,j,:) - (norma-zerg);
                     end
               end 
            end

            sq = des(:);%./max(des(:),[],'all');
            vec_of_vecs(totp,:) =  sq;
            %vec_of_vecs(totp,:) = uint8(255*(des(:)./max(des(:),[],'all')));
            %vec_of_vecs(totp,:) = uint16(65535*(des(:)./max(des(:),[],'all')));
            points(totp,:) = [x,y];
        end
    end
end

vec_of_vecs = vec_of_vecs(1:totp,:);
%level = graythresh(vec_of_vecs);
%vec_of_vecs =  imbinarize(vec_of_vecs);
%histogram(vec_of_vecs)

end


