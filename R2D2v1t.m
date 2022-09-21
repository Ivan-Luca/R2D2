function [vec_of_vecs, points,rapports] = R2D2v1(im,kp,area_size,patch_size,type)
if nargin < 1
    area_size = 96;
    patch_size = 3;%3;
end
base = 20;
nkp = size(kp,1);

[L,C] = size(im);
im = rescale(double(im));
mindimage = zeros(L,C,1,'double');

kernel = [ 1/16 1/8 1/16; 1/8 1/4 1/8; 1/16 1/8 1/16];
kernel2 = (1/273)*[1 4 7 4 1; 1 16 26 16 4; 7 26 41 26 7; 1 16 26 16 4; 1 4 7 4 1];


im2 = imgaussfilt(im,0.75);

[Ix, Iy, Ixy, Iyx, C2, D2] = imcurl(im2,0,0,100);



mindimage(:,:,1) = C2/max(C2,[],'all');
mindimage(:,:,2) = D2/max(D2,[],'all');



%iterate through kps
totp = 0;
ns = 16;
kernel = fspecial('gaussian', [95 95], 40);
kernel = (kernel/max(kernel,[],'all'))

rapports =zeros([nkp 1]);
for m = 1:nkp %nkp
    
    y = kp(m,2);
    x = kp(m,1);

    if (x > (area_size+1))  && (x < (C-area_size))
        if (y > (area_size+1)) && (y < (L-area_size))
            totp=totp+1;
            sqt = mindimage(round(y-area_size:y+area_size),round(x-area_size:x+area_size),:) .* kernel;
            [ys,xs,zs] = size(sqt);
            covM = cov(sqt(:,:,1),sqt(:,:,2));
            f1 = covM(1,1)/(covM(1,1)+covM(2,2));
            f2 = 1-f1;
            rapports(m) = f1;
            sqt2 = sqt(:,:,1)*f1 + sqt(:,:,2)*f2;

            for j = 1:ns
                for i = 1:ns
                     clip = sqt2(round((j-1)*ys/ns+1):round(j*ys/ns),round((i-1)*xs/ns+1):round(i*xs/ns),:);
                     if type == 1
                        clip = (log(clip)/log(base))+1;
                        clip = max(clip,0);
                     end
                     
                     des(i,j) = mean(clip,'all');

               end 
            end
            
            sq = uint8(255*(des(:)./max(des(:),[],'all')));
            vec_of_vecs(totp,:) = sq;

            points(totp,:) = [x,y];
        end
    end
end

vec_of_vecs = vec_of_vecs(1:totp,:);
if type == 1
    final_image = (log(mindimage(:,:,1)*mean(rapports,'all') + mindimage(:,:,2)*(1 - mean(rapports,'all')))/log(base))+1;
    final_image = max(final_image,0);
else
    final_image = (mindimage(:,:,1)*mean(rapports,'all') + mindimage(:,:,2)*(1 - mean(rapports,'all')));
end
final_image = final_image/max(final_image,[],'all');
figure, imshow(kernel);
figure, imshow( uint8(255*(final_image)));

end


