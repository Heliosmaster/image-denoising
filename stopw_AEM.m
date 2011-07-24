function stopw_AEM(dir,im,beta,IT)
% stopw_AEM(dir,im,beta,it)
%
% repeat the restoration of the given image im
% in the given directory dir, using the given beta
% IT times, in order to evaluate precisely the
% time required by AEM for such image.


load([dir '/' im]);
load([dir '/' im '_noisy']);
load([dir '/' 'sol' im]);
load([dir '/' 'data_' im ],'TimeCost');

fprintf('Working on %s\n',im);

T1 = TimeCost;

[n,n] = size(gn);

grad1  = @(y)[y(:,2:n)-y(:,1:n-1), y(:,1)-y(:,n)];
grad2  = @(y)[y(2:n,:)-y(1:n-1,:); y(1,:)-y(n,:)];
div    = @(x1,x2)([-x1(:,1)+x1(:,n),-x1(:,2:n)+x1(:,1:n-1)] +...
    [-x2(1,:)+x2(n,:);-x2(2:n,:)+x2(1:n-1,:)]);
H = @(x)x;
HT = @(x)x;

bg = 0;
delta = 0;
NIT = 3000;
tol = 0;

zeroindex =gn <= 0;
nonzeroindex = ~zeroindex;
eta = min(gn(nonzeroindex))*ones(size(gn));
eta(zeroindex) = 0;

verbose=false;
Nalpha=10;
epsilon =1e-4;

benchsol={obj,x};

%TF = ones(size(obj));
for s=1:IT
    [x,y1,y2,y3,TimeCost,InnVec,Primal,phi_xy,alpha_vec,err,...
        KKT,KL] = AEM(gn, H, HT,bg, beta, delta, eta, grad1,...
        grad2, div, NIT, tol, verbose, Nalpha, epsilon,...
        benchsol);
    
    fprintf('%g) Current: %g, Fastest: %g\n',s,...
        TimeCost(end),T1(end));
    if TimeCost(end) < T1(end)
        save([dir '/' 'data_'  im ]);
        T1 = TimeCost;
    end
    clear x y1 y2 y3 TimeCost InnVec Primal phi_xy ...
        alpha vec err KKT KL;
end
