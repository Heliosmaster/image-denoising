% trial of given values for beta for regularization
% using AEM
clear all;

for k=[11,13,22,48,49,50]
    
    image = ['PRISMA/Image' int2str(k) '.tif'];
    gn = 4095-double(imread(image));
    gn = gn(:,138:404);
    
    obj = 4095-double(imread('PRISMA/Avg.tif'));
    obj = obj(:,138:404);
    
    [r] = size(gn,1);
    grad1  = @(y)[y(:,2:r)-y(:,1:r-1), y(:,1)-y(:,r)];
    grad2  = @(y)[y(2:r,:)-y(1:r-1,:); y(1,:)-y(r,:)];
    div    = @(x1,x2)([-x1(:,1)+x1(:,r),-x1(:,2:r)+x1(:,1:r-1)] + ...
        [-x2(1,:)+x2(r,:);-x2(2:r,:)+x2(1:r-1,:)]);
    H = @(x)x;
    HT = @(x)x;
    
    scale = 2/(r^2);
    
    bg = 0;
    delta = 0;
    NIT = 3000;
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
    
    %tolbeta = 1e-3;
    zeroindex = gn <= 0;
    nonzeroindex = ~zeroindex;
    rapp = zeros(size(gn));
    den = H(obj)+bg;
    rapp(nonzeroindex) = gn(nonzeroindex)./ den(nonzeroindex);
    
    KLstima = sum( gn(nonzeroindex).* log(rapp(nonzeroindex)) ...
        + den(nonzeroindex) - gn(nonzeroindex) );
    KLstima = KLstima + sum(den(zeroindex));
    KLstima=KLstima*scale;
    fprintf('KLSTIMA=%g\n',KLstima);
    
    %%% choice of beta
    %beta_start = 0.01;
    %beta_end = 0.02;
    %beta=beta_start:0.01:beta_end;
    beta=[0.01 0.02 0.05 0.07 0.1 0.2 0.3];
    
    %sol=zeros(r,r,length(beta));
    
    for i=1:length(beta)
        [~,~,~,~,~,~,~,phi_xy,~,...
            err,~,KL] = AEM(gn, H, HT,bg, beta(i), delta,...
            eta, grad1, grad2,div, NIT, tol, verbose, Nalpha,...
            epsilon, benchsol);
        %ratio = (KL(end)*scale)-1;
        %sol(:,:,i)=x;
        TV=(phi_xy(end)-KL(end))/beta(i);
        out = [out; beta(i), err{1}(end), KL(end),TV];
        fprintf('beta: %g, err:%g, KL: %g  TV=%g\n',beta(i),...
            err{1}(end),KL(end),TV);
    end
    
    %format long;
    %disp(out);
    save([ 'PRISMA/array2_beta_im' int2str(k)]);
    
end
exit;
