function iter_AEM(dir,im,beta,NIT)
% iter_AEM(dir,im,beta,NIT)
%
% AEM for the given image im in the given directory dir
% with the given beta, with NIT iterations

load([dir '/' im]);
load([dir '/' im '_noisy']);
load([dir '/' 'sol' im]);

[n,n] = size(gn);

grad1  = @(y)[y(:,2:n)-y(:,1:n-1), y(:,1)-y(:,n)];
grad2  = @(y)[y(2:n,:)-y(1:n-1,:); y(1,:)-y(n,:)];
div    = @(x1,x2)([-x1(:,1)+x1(:,n),-x1(:,2:n)+x1(:,1:n-1)] + ...
    [-x2(1,:)+x2(n,:);-x2(2:n,:)+x2(1:n-1,:)]);
H = @(x)x;
HT = @(x)x;

bg = 0;
delta = 0;
tol = 0;

zeroindex =gn <= 0;
nonzeroindex = ~zeroindex;
eta = min(gn(nonzeroindex))*ones(size(gn));
eta(zeroindex) = 0;

verbose=false;
Nalpha=10;
epsilon =1e-4;

benchsol={obj,x};

[u,y1,y2,y3,TimeCost,InnVec,Primal,phi_xy,alpha_vec,err,KKT,KL]=...
    AEM(gn, H, HT,bg, beta, delta, eta, grad1, grad2, div, NIT,...
    tol, verbose, Nalpha, epsilon, benchsol);

save([dir '/' 'data_'  im '_it_' int2str(NIT)]);