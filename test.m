clear all;
tic; %initialization of the timer
image = 'Image13.tif';
%load('Dati/airplane.mat');
%load('Dati/airplane_noisy.mat');

gn = 4095-double(imread(image)); 
gn = gn(:,138:404); %% extraction of a square portion of the image

obj = 4095-double(imread('Avg.tif'));
obj = obj(:,138:404); %% same as before

[r,c] = size(gn);

grad1  = @(y)[y(:,2:c)-y(:,1:c-1), y(:,1)-y(:,c)];
grad2  = @(y)[y(2:r,:)-y(1:r-1,:); y(1,:)-y(r,:)];
div    = @(x1,x2)([-x1(:,1)+x1(:,c),-x1(:,2:c)+x1(:,1:c-1)] + [-x2(1,:)+x2(r,:);-x2(2:r,:)+x2(1:r-1,:)]);
H = @(x)x;
HT = @(x)x;

scale = 2/(r^2);

bg = 0;
delta = 0;
NIT = 2000;
tol = 1e-10; %% very small so NIT is reached

zeroindex = gn <= 0;
nonzeroindex = ~zeroindex;
eta = min(gn(nonzeroindex))*ones(size(gn));
eta(zeroindex) = 0;

verbose =false;
Nalpha=10;
epsilon = 1e-4; 

benchsol={obj};

out = [];

tolbeta = 1e-3;

beta_start = 0.001;
beta_end = 0.5;

zeroindex = gn <= 0;
nonzeroindex = ~zeroindex;
rapp = zeros(size(gn));
den = H(obj)+bg;%ifft2(fft2(y).*TF)+bg;
rapp(nonzeroindex) = gn(nonzeroindex)./ den(nonzeroindex);
KLstima = sum( gn(nonzeroindex).* log(rapp(nonzeroindex)) + den(nonzeroindex) - gn(nonzeroindex) );
KLstima = KLstima + sum(den(zeroindex));
KLstima=KLstima*scale;
fprintf('KLSTIMA=%g\n',KLstima);
%%% first and last values are computed to initialize the bisection loop
[x,y1,y2,y3,TimeCost,InnVec,Primal,phi_xy,alpha_vec,err,KKT,KLs] = ...
    AEM(gn, H, HT,bg, beta_start, delta, eta, grad1, grad2, div, NIT, tol, verbose, Nalpha, epsilon, benchsol);
ratio_start = (KLs(end)*scale)-KLstima;
out = [out; beta_start, err{1}(end), KLs(end), phi_xy(end), KLs(end)*scale];
[x,y1,y2,y3,TimeCost,InnVec,Primal,phi_xy,alpha_vec,err,KKT,KLe] = ...
    AEM(gn, H, HT,bg, beta_end, delta, eta, grad1, grad2, div, NIT, tol, verbose, Nalpha, epsilon, benchsol);
ratio_end = (KLe(end)*scale)-KLstima;
out = [out; beta_end, err{1}(end), KLe(end), phi_xy(end), (KLe(end)*scale)];

ratio = ratio_start;

while(abs(beta_start-beta_end)>=tolbeta)
    beta_half = (beta_start+beta_end)/2;
    [x,y1,y2,y3,TimeCost,InnVec,Primal,phi_xy,alpha_vec,err,KKT,KL] = AEM(gn, H, HT,bg, beta_half, delta, eta, grad1, grad2, div, NIT, tol, verbose, Nalpha, epsilon, benchsol);
    ratio_half = (KL(end)*scale)-KLstima;
    out = [out; beta_half, err{1}(end), KL(end), phi_xy(end), (KL(end)*scale)];
    fprintf('beta_start: %g, beta_half:%g, beta_end: %g\n',beta_start,beta_half,beta_end);
    fprintf('ratio_start: %g, ratio_half: %g, ratio_end: %g\n\n',ratio_start,ratio_half,ratio_end); 
    if (ratio_half * ratio_end < 0)
        beta_start = beta_half;
        ratio_start = ratio_half;
    elseif(ratio_start * ratio_half < 0)
        beta_end = beta_half;
        ratio_end = ratio_half;
    else
        break;
    end
    ratio = ratio_half;
 
end

fprintf('Optimal beta value for : %g, with a tolerance of: %g\n',beta_half,tolbeta);
fprintf('time elapsed: %g seconds \n', toc);

disp('beta,error,KL,phi_xy,ratio');
disp(out);
%save outairplane out x KLstima

%out = sortrows(out,1);