clear all
% PIDSplit+ algorithm 
% periodic boundary conditions

%close all

scelta = 'denoising';
boundary = 'periodic';
err_min = 1;

if strcmp(scelta,'deblurring');
    problem   = '../../Dati/micro.mat';
    load(problem);% carica obj, gn e psf
    delta = 1e+000; beta=  0.164050701950861;
    %delta = 1e-002; beta=   0.105598009481053;
    %delta = 1e-004; beta=   0.128466499832609;
    %delta = 1e-006; beta=   0.158973451282237;
    %delta = 1e-008; beta=   0.164573664207980;
    if err_min
        load solmicro.mat
        yopt = x; clear x;
        normyopt = norm(yopt(:));
    end
    TF  = fft2(fftshift(psf));
    eta = zeros(size(gn));
    zeroindex = gn <= eta;
    nonzeroindex = ~zeroindex;
    K  = @(x)ifft2(fft2(x).*TF);
    KT = @(x)ifft2(fft2(x).*conj(TF));
    boundary = 'periodic';
elseif strcmp(scelta,'denoising');

    noisy = 'circles/circles_noisy.mat'; % carica gn e obj
    sol   = 'circles/circles.mat';
    fprintf('%s \n',noisy);
    load(noisy);% carica obj e gn
    load(sol);
    if err_min
        load circles/solcircles.mat
        yopt = x; clear x
        normyopt = norm(yopt(:));
    end

    TF = ones(size(gn));
    zeroindex = gn <= 0;
    nonzeroindex = ~zeroindex;
    eta = min(gn(nonzeroindex))*ones(size(gn));
    eta(zeroindex) = 0;
    delta = 0;
    beta = 0.25;
    bg = 0;
    K  = @(x)x;
    KT = @(x)x;
else
    error('denoising o deblurring?');
end

n = length(gn);                %Assume a square image
if strcmp(boundary,'periodic');
    grad1  = @(y)[y(:,2:n)-y(:,1:n-1), y(:,1)-y(:,n)];
    grad2  = @(y)[y(2:n,:)-y(1:n-1,:); y(1,:)-y(n,:)];
    div    = @(x1,x2)([-x1(:,1)+x1(:,n),-x1(:,2:n)+x1(:,1:n-1)] + [-x2(1,:)+x2(n,:);-x2(2:n,:)+x2(1:n-1,:)]);
elseif strcmp(boundary,'reflexive')
    grad1  = @(y)[y(:,2:n)-y(:,1:n-1), zeros(n,1)];
    grad2  = @(y)[y(2:n,:)-y(1:n-1,:); zeros(1,n)];
    div    = @(x1,x2)([-x1(:,1),-x1(:,2:n)+x1(:,1:n-1)] + [-x2(1,:);-x2(2:n,:)+x2(1:n-1,:)]);
else
    error('periodic or reflexive boundary conditions?');
end
    
%beta = 0.164050701950861;%0.07;2*TV(obj,0)/norm(obj(:)-gn(:))^2;

NIT = 1000;
tol = 0;
verbose = 1;
if err_min
    benchsol = {obj,yopt};
else
    benchsol = {obj};
end
gamma = 5/beta;
[u,w3,TimeCost,fobj,err] = ...
    PIDSplit_plus(gn, TF, K, KT, bg, beta, gamma, eta, grad1, grad2, div, NIT, tol, verbose, benchsol);
