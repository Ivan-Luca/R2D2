function [tform, inliersIndex] = estimateGeoTransPente(moving, fixed, method,distance, angle)

pentes = (moving(:,2) - fixed(:,2))/(moving(:,1) - fixed(:,1));
norms = [sqrt((fixed(:,1) - moving(:,1)).^2) + sqrt((fixed(:,2) - moving(:,2)).^2)];
mpentes = mean(pentes);
mnorms = mean(norms);

inliersIndex = zeros(length(fixed));
for i= 1:length(fixed)
    if
end

