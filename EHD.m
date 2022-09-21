function res = EHD(im, kps) 
            eh = zeros(80, size(kps,2));
            kps_to_ignore = zeros(1,size(kps,2));
            patch_size = 100;
            for i = 1: size(kps,2)
                % Patch location
                x = round(kps(1, i));
                y = round(kps(2, i));
                % Top-left point of patch
                x1 = max(1,x-floor( patch_size/2));
                y1 = max(1,y-floor( patch_size/2));
                x2 = min(x+floor( patch_size/2),size(im,2));
                y2 = min(y+floor( patch_size/2),size(im,1));
                % Ignore incomplete patches
                if y2-y1 ~=  patch_size || x2-x1 ~=  patch_size
                    kps_to_ignore(i)=1;
                    continue;
                end  
                %[eh1]=ehd(edge_map(y1:y2, x1:x2),[],3,0);  
                [eh1]=ehd(im(y1:y2, x1:x2),[],3,0);  
                eh(:,i)= eh1;
            end           
            res = struct('kps', kps(:,kps_to_ignore ==0)', 'des', eh(:,kps_to_ignore==0)');
        end