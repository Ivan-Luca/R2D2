function [vec_of_vecs, points,rapports,m1,m2] = R2D2v1(im,kp,area_size,patch_size)

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
CD = cat(3,m1,m2);
% figure,imshow(CD(:,:,1)),colormap('turbo'),title('CCD');
% figure,imshow(CD(:,:,2)),colormap('turbo'),title('DCD');
% figure,imshow(CD(:,:,1)*0.3 + CD(:,:,2)*0.7),colormap('turbo');
%iterate through kps
totp = 0;
ns = 16;
ps = floor((2*area_size)/ns);
kernel = genGausKernel(ps,(ps-1)/6);
kernel2 = genGausKernel(2*area_size+1,(2*area_size)/6);
nb = sum(kernel ~= 0,'all');
des = zeros(ns,ns);
vec_of_vecs = zeros(nkp,ns*ns);
points = zeros(nkp,2);
rapports =zeros([nkp 1]);

% figure,imshow(kernel)
% figure,imshow(kernel2)
for m = 1:nkp %nkp
    
    y = kp(m,2);
    x = kp(m,1);

    if (x > (area_size+1))  && (x < (C-area_size))
        if (y > (area_size+1)) && (y < (L-area_size))
            totp=totp+1;
            sqt = CD(round(y-area_size:y+area_size),round(x-area_size:x+area_size),:);
            %dom = atan2(sum(Iyy(round(y-area_size:y+area_size),round(x-area_size:x+area_size),:),'all'),sum(Ixx(round(y-area_size:y+area_size),round(x-area_size:x+area_size),:),'all'));
            [ys,xs,zs] = size(sqt);
            covM = cov(sqt(:,:,1),sqt(:,:,2));
            f1 = covM(1,1)/(covM(1,1)+covM(2,2));
            f2 = 1-f1;
            rapports(m) = f1;
            sqt2 = sqt(:,:,1)*f2+ sqt(:,:,2)*f1;
            %sqt2 = imresize(sqt2,[16 16]);
            %sqt2 = sqt(:,:,1)./sqt(:,:,2);
            %sqt2 = CD(round(y-area_size:y+area_size),round(x-area_size:x+area_size),2);
            %size(sqt2)
            %patch = sqt2(round(((ns/2)-1)*ys/ns+1):round(((ns/2)-1)*ys/ns+(ns/2)),round(((ns/2)-1)*xs/ns+1):round(((ns/2)-1)*xs/ns+(ns/2)),:);
            for j = 1:ns
                for i = 1:ns
                    clip = sqt2(floor((j-1)*(ys/ns))+1:floor((j-1)*(ys/ns))+ps,floor((i-1)*(xs/ns))+1:floor((i-1)*(xs/ns))+ps,:);
                    clip = clip.*kernel;
                    
                    des(i,j) = sum(clip.*kernel,'all');
               end 
            end
            %figure, subplot(1,2,1),imshow(des), hold on
            %des = imrotate(des,-(dom*180)/pi,'crop');
           % subplot(1,2,2),imshow(des);
            %sq = des(:)./max(des(:),[],'all');
            %sq = compress(des(:));
            %sq = dec2bin(compress(des(:)));
            vec_of_vecs(totp,:) =  des(:)./max(des(:),[],'all');
            %vec_of_vecs(totp,:) = uint8(255*(des(:)./max(des(:),[],'all')));
            %vec_of_vecs(totp,:) = uint16(65535*(des(:)./max(des(:),[],'all')));
            points(totp,:) = [x,y];
        end
    end
end

%vec_of_vecs = vec_of_vecs(1:totp,:);

vec_of_vecs = vec_of_vecs(1:totp,:);
%level = graythresh(vec_of_vecs);
%vec_of_vecs =  imbinarize(vec_of_vecs);
%histogram(vec_of_vecs)

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


