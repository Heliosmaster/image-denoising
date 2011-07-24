clear all
close all

Imax = 1000;
obj = double(imread('cameraman.tif'));
objmax = max(obj(:));
Umax = Imax/objmax;
obj = obj*Umax;
% x = [-128:127];
% [X,Y] = meshgrid(x);
% sigma = 1.3;
% psf = exp(-X.^2/(2*sigma^2) -Y.^2/(2*sigma^2))/(2*pi*sigma^2);
% figure, imshow(psf,[]); colorbar; title(['Normalized Gaussian psf, \sigma = ',num2str(sigma)]);
% fprintf('\n Normalizzazione = %g\n',sum(sum(psf)));
% bg = 0;
% TF = fft2(fftshift(psf));
% g = ifft2(fft2(obj).*TF) + bg;
g=obj;
%figure, imshow(g/Umax,[]); colorbar; title('Blurred obj');
g = g/1e12;
gn = imnoise(g,'poisson')*1e12;
obj = obj/Umax; gn = gn/Umax;
figure, imshow(obj,[]); colorbar; title('obj');
figure, imshow(gn,[]); colorbar; title('gn');
save cameraman_1000-2.mat obj gn 
