% bisection method to solve the nonlinear
% equation in beta for the regularization
clear all;

for i=[11,13,22,48,49,50]
    
    image = ['PRISMA/Image' int2str(i) '.tif'];
    gn = 4095-double(imread(image));
    gn = gn(:,138:404);
    
    obj = 4095-double(imread('PRISMA/Avg.tif'));
    obj = obj(:,138:404);
    
    [n,n] = size(gn);
    
    grad1  = @(y)[y(:,2:n)-y(:,1:n-1), y(:,1)-y(:,n)];
    grad2  = @(y)[y(2:n,:)-y(1:n-1,:); y(1,:)-y(n,:)];
    div    = @(x1,x2)([-x1(:,1)+x1(:,n),-x1(:,2:n)+...
        x1(:,1:n-1)] +[-x2(1,:)+x2(n,:);-x2(2:n,:)+...
        x2(1:n-1,:)]);
    H = @(x)x;
    HT = @(x)x;
    
    scale = 2/(n^2);
    
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
    
    tolbeta = 1e-3;
    
    beta_start = 0.001;
    beta_end = 0.5;
    
    zeroindex = gn <= 0;
    nonzeroindex = ~zeroindex;
    rapp = zeros(size(gn));
    den = H(obj)+bg;
    rapp(nonzeroindex) = gn(nonzeroindex)./ den(nonzeroindex);
    
    KLstima = sum( gn(nonzeroindex).* log(rapp(nonzeroindex))...
        + den(nonzeroindex) - gn(nonzeroindex) );
    KLstima = KLstima + sum(den(zeroindex));
    KLstima=KLstima*scale;
    
    fprintf('KLstima=%g\n',KLstima);
    
    tic;
    
    % first and last values are computed to initialize the 
    % bisection loop
    [x,y1,y2,y3,TimeCost,InnVec,Primal,phi_xy,alpha_vec,err,...
        KKT,KLs]=AEM(gn, H, HT,bg, beta_start, delta, eta, ...
        grad1, grad2, div,...
        NIT, tol, verbose, Nalpha, epsilon, benchsol);
    ratio_start = (KLs(end)*scale)-1;
    out = [out; beta_start, err{1}(end), KLs(end), ...
        phi_xy(end), KLs(end)*scale];
    [x,y1,y2,y3,TimeCost,InnVec,Primal,phi_xy,alpha_vec,err,...
        KKT,KLe]=AEM(gn, H, HT,bg, beta_end, delta, eta,...
        grad1, grad2, div, ...
        NIT, tol, verbose, Nalpha, epsilon, benchsol);
    ratio_end = (KLe(end)*scale)-1;
    out = [out; beta_end, err{1}(end), KLe(end), ...
        phi_xy(end),KLe(end)*scale];
    
    ratio = ratio_start;
    
    while(abs(beta_start-beta_end)>=tolbeta)
        beta_half = (beta_start+beta_end)/2;
        [x,y1,y2,y3,TimeCost,InnVec,Primal,phi_xy,alpha_vec,...
            err,KKT,KL]=AEM(gn, H, HT,bg, beta_half, delta,...
            eta, grad1, grad2, div,...
            NIT, tol, verbose, Nalpha, epsilon, benchsol);
        ratio_half = (KL(end)*scale)-1;
        out = [out; beta_half, err{1}(end), KL(end),...
            phi_xy(end),ratio_half];
        fprintf('beta_start:%g, beta_half:%g, beta_end:%g\n',...
            beta_start,beta_half,beta_end);
        fprintf('ratio_start:%g,ratio_half:%g,ratio_end:%g\n\n',...
            ratio_start,ratio_half,ratio_end);
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
    
    fprintf('For Image %g \n',i);
    fprintf('Optimal beta value:%g,with a tolerance of:%g\n',...
        beta_half,tolbeta);
    fprintf('time elapsed: %g seconds \n', toc);
    
    save([ 'PRISMA/beta_im' int2str(i)]);
    
end
exit;
