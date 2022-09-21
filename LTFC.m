function outliersIndex = LTFC(im1,im2,matchedPoints1,matchedPoints2)


matchedPoints2(:,1) = matchedPoints2(:,1) + size(im2,2);
pentes = (matchedPoints2(:,2)-matchedPoints1(:,2))./(matchedPoints2(:,1)-matchedPoints1(:,1));
dist =  sqrt((matchedPoints2(:,2)-matchedPoints1(:,2)).^2 + (matchedPoints2(:,1)-matchedPoints1(:,1)).^2);
vect = [pentes dist];
[vect2, outliersIndex] = rmoutliers(vect)
figure,subplot(1,2,1), plot3([1:size(vect,1)],vect(:,1),vect(:,2));
subplot(1,2,2), plot3([1:size(vect2,1)],vect2(:,1),vect2(:,2)),linkaxes();

