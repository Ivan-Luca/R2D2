function cornerPoints = textFilter(im,points,thresh)

img2 = imgaussfilt(normalize(single(im)),5/3);
img2 = imgaussfilt(img2,5/3);

cornerPoints = [points(1)];
[L,C] = size(img2);
length = size(points);
for kp = 1:length
    i = round(points.Location(kp,2)); 
    j = round(points.Location(kp,1));
    if (5 < i) && (i < L-5) 
        if (5 < j) && (j < C-5)
            dp1 = img2((i-3)-1:(i-3)+2,j-1:j+2);
            dp2 = img2((i+3)-1:(i+3)+2,j-1:j+2);
            dp3 = img2(i-1:i+2,(j-3)-1:(j-3)+2);
            dp4 = img2(i-1:i+2,(j+3)-1:(j+3)+2);
            dpx = abs(sum(sum(dp1))-sum(sum(dp2)));
            dpy = abs(sum(sum(dp3))-sum(sum(dp4)));
            if (thresh < dpx) || (thresh < dpy)
                if (dpx+dpy) > 2.35*thresh
                    cornerPoints = [cornerPoints; points(kp)];
                end
            end
        end
    end
    %cornerPoints = cornerPoints(2:size(cornerPoints));
end

end