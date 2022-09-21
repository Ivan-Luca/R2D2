function [vec_of_vecs, points] = MINDescriptorLQ2(im,kp,area_size,patch_size)
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
mindimage = zeros(L,C,8,'single');
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


mindimage = MIND_descriptor2D(im,1);
mindimage(:,:,:) = rescale(mindimage(:,:,:),0,1);
% 'before'
% figure
% imshow(mindimage(:,:,1:1))
% for p = 1:8
%     mindimage(:,:,p) = imbinarize(mindimage(:,:,p),0.5);
% end
% 'after'
% figure
% imshow(mindimage(:,:,1:1))
%mindimage(:,:,i) = imerode(mindimage(:,:,i),se);
sq = zeros(2*area_size+1,2*area_size+1,8,'uint8');
%iterate through kps
totp = 0;
for m = 1:nkp %nkp
    
    y = kp(m,2);
    x = kp(m,1);
    if (x > 2*area_size)  && (x < (C-2*area_size))
        if (y > 2*area_size) && (y < (L-2*area_size))
            totp=totp+1;
            
            sq = imbinarize(mindimage(round(y-area_size:y+area_size),round(x-area_size:x+area_size),:),0.5);
%             for n = 1:*area_size+1
%                 for k = 1:*area_size+1
%                     for p = 1:8
%                    sq(n,k,p) = quantityMIND(sqt(n,k,p))
%                 end
%             end
            
            %figure
            %imshow(sq(:,:,1:3))
            arr = uint8(reshape(sq,1,[]));
            %arr = pow2(7:-1:0)*reshape(sq,8,[]);
%             arr2 = zeros(1,size(arr,2)/8,'uint8');
%             step = 1;step2 = 1;
%             while step < size(arr2,2)
%                 arr2(1,step2) = bi2de(arr(step:step+7));
%                 step=step+8;step2 = step+1;
%             end
            
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
            vec_of_vecs(totp,:) = arr;
            %vec_of_vecs2(totp,:) = arr2;
            points(totp,:) = [x,y];
        end
    end
end

vec_of_vecs = binaryFeatures(vec_of_vecs(1:totp,:));
%vec_of_vecs2 = vec_of_vecs2(1:totp,:);
points = points(1:totp,:);
% figure
% imshow(mindimage(:,:,1:1))
% figure
% imshow(mindimage(:,:,4:4))
% figure
% imshow(mindimage(:,:,6:8))

end


