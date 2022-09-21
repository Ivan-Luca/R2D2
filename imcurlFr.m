function [Ixx, Iyy, Ixy, Iyx , C, D] = imcurlFr( image)

G = rescale(double(im2gray(image)));
Fh = [-1 0 1];
Fv = [-1 0 1]';

Ipad = padarray(G,[2 2],0,"post"); % zero padding
Fpadh = padarray(Fh,[size(G,1)+1 size(G,2)-1],0,"post");
Fpadv = padarray(Fv,[size(G,1)-1 size(G,2)+1],0,"post");% zero padding
size(Ipad)
size(Fpadh)
Ifr = fft2(Ipad);
Ffrh = fft2(Fpadh);
Ffrv = fft2(Fpadv);

Ixfr = Ifr .* Ffrh;
Iyfr = Ifr .* Ffrv;
Ixx = ifft2(Ixfr .* Ffrh);
Iyy = ifft2(Iyfr .* Ffrv);
Ixy = ifft2(Ixfr .* Ffrv);
Iyx = ifft2(Iyfr .* Ffrh);

%Opad = ifft2(Offt);
%O = Opad(2:end-1,2:end-1);
% store the values in our output variables, for clarity

D = abs(Ixx + Iyy);
C = abs(Ixy-Iyx);
C = C(2:end-1,2:end-1);
D = D(2:end-1,2:end-1);


figure, imshow(Ixx),colormap('Turbo'),title('Ixx');
figure, imshow(Iyy),colormap('Turbo'),title('Iyy');
figure, imshow(Ixy),colormap('Turbo'),title('Ixy');
figure, imshow(Iyx),colormap('Turbo'),title('Iyx');
figure, imshow(C),colormap('Turbo'),title('C');
figure, imshow(D),colormap('Turbo'),title('D');
end

