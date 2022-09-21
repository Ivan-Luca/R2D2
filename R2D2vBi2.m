function [vec_of_vecs, points,rapports,m1,m2] = R2D2vBi2(im,kp,area_size,patch_size)

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

des = zeros(248,1,'uint8');
vec_of_vecs = zeros(nkp,248,'uint8');
points = zeros(nkp,2);
rapports =zeros([nkp 1]);
pas = ((area_size-15)-3)/4

radius = [round(3:pas:(area_size-7))]

rng('default');
rng(1); % reapetable random
it=1;
X = [];
Y = [];
for m = 1:nkp %nkp
    
    y = kp(m,2);
    x = kp(m,1);

    if (x > (area_size+1))  && (x < (C-area_size))
        if (y > (area_size+1)) && (y < (L-area_size))
            totp=totp+1;
            it=1;
            noyau = 5;
            step = pi/4;
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
            for j = 1:5
                angles = [0:step:2*pi];
                pairs = round([radius(j)*cos(angles)+(area_size+1); radius(j)*sin(angles)+(area_size+1)]');
                X = [X pairs(:,1)']; Y = [Y pairs(:,2)'];
                kernel = genGausKernel(noyau,(ps-1)/6);
                samples = 2^(j+2);
                for i = 1:samples
                    x1 = pairs(i,1); x2 = pairs(i+1,1); y1 = pairs(i,2); y2 = pairs(i+1,2);
                    patch = sqt2(y1-floor(noyau/2):floor(noyau/2)+y1,x1-floor(noyau/2):x1-floor(noyau/2),:);
                    %clip = sqt2(y2-floor(noyau/2):floor(noyau/2)+y2,x2-floor(noyau/2):x2-floor(noyau/2),:);
                    clip = sqt2(area_size+1-floor(noyau/2):floor(noyau/2)+area_size+1,area_size+1-floor(noyau/2):area_size+1-floor(noyau/2),:);

                    
                    m = sum(patch(:,:).*kernel,'all') - sum(clip(:,:).*kernel,'all');
                    if  m(:,:) <= 0
                        des(it) = 0;
                    else
                        des(it) = 1;
                    end
                    it = it+1;
                end 
                
                step = step/2;
                noyau = noyau+4;
               
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
figure, plot(X,Y);
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
