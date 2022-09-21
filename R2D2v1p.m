function [vec_of_vecs, points,rawpoints,rapports] = R2D2v1p(im,area_size,patch_size,amount,border)
if nargin < 1
    area_size = 96;
    patch_size = 3;%3;
    border = 10;
    amount = 6000;
end
base = 20;


[L,C] = size(im);
im = rescale(double(im));
mindimage = zeros(L,C,1,'double');

im2 = imgaussfilt(im,0.75);

[Ix, Iy, Ixy, Iyx, C2, D2] = imcurl(im2,0,0,100);



mindimage(:,:,1) = C2/max(C2,[],'all');
mindimage(:,:,2) = D2/max(D2,[],'all');
CD = (mindimage(:,:,1) + mindimage(:,:,2));
figure,subplot(1,2,1),imshow(mindimage(:,:,1))
subplot(1,2,2),imshow(mindimage(:,:,2))
Cid = mindimage(:,:,1)>0;
Did = mindimage(:,:,2)>0.1;
CD = C2+D2;
CDid = (Cid+Did) > 0.1;
figure,imshow(Cid);
zero = zeros(size(Cid));
figure,imshow(zero);
zero(20:(end-20), 20:(end-20)) = 1;
figure,imshow(zero);
Cid = Cid.*zero;
CDid = CDid.*zero;
figure,imshow(CDid);
[row, col] = find(CDid);
ind = sub2ind(size(C2), row, col);
scores = CD(ind);
rawpoints = cornerPoints([col row],'Metric',scores);
%kp = [rawpoints.Location(:,1) rawpoints.Location(:,2)]
kp= rawpoints.selectStrongest(amount);
kp = [kp.Location(:,1) kp.Location(:,2)];
figure, imshow(CDid),hold on, plot(col,row,'o'),hold off;
figure, imshow(im),hold on, plot(rawpoints.Location(:,1),rawpoints.Location(:,2),'o'),hold off;
figure, imshow(im),hold on, plot(kp(:,1),kp(:,2),'o'),hold off;
nkp = size(kp,1);
%iterate through kps
totp = 0;
ns = 16;


rapports =zeros([nkp 1]);
for m = 1:nkp %nkp
    
    y = kp(m,2);
    x = kp(m,1);

    if (x > (area_size+1))  && (x < (C-area_size))
        if (y > (area_size+1)) && (y < (L-area_size))
            totp=totp+1;
            sqt = mindimage(round(y-area_size:y+area_size),round(x-area_size:x+area_size),:);
            [ys,xs,zs] = size(sqt);
            covM = cov(sqt(:,:,1),sqt(:,:,2));
            f1 = covM(1,1)/(covM(1,1)+covM(2,2));
            f2 = 1-f1;
            rapports(m) = f1;
            sqt2 = sqt(:,:,1)*f1 + sqt(:,:,2)*f2;

            for j = 1:ns
                for i = 1:ns
                     clip = sqt2(round((j-1)*ys/ns+1):round(j*ys/ns),round((i-1)*xs/ns+1):round(i*xs/ns),:);

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


end


