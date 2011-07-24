function app(dir,im,beta)
% app(dir,im,beta)
% computes the AEM and the PID50,PID5,PID1,PID05
% on the given image in the given directory with the given beta
% and saves the results in an appropriate .mat file

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

%scale = 2/(r^2);

bg = 0;
delta = 0;
NIT = 3000;
tol = 1e-6;

zeroindex =gn <= 0;
nonzeroindex = ~zeroindex;
eta = min(gn(nonzeroindex))*ones(size(gn));
eta(zeroindex) = 0;

verbose=true;
Nalpha=10;
epsilon = 1e-4; 

benchsol={obj,x};

TF = ones(size(obj));

tic;

[x,y1,y2,y3,TimeCost,InnVec,Primal,phi_xy,alpha_vec,err,KKT,...
    KL] = AEM(gn, H, HT,bg, beta, delta, eta, grad1, grad2, ...
    div, NIT, tol,verbose, Nalpha, epsilon, benchsol);
save([dir '/' 'data_' im]);
clear x y1 y2 y3 TimeCost InnVec Primal phi_xy alpha vec err KKT KL; 

gamma = 50/beta;
[u_50,w3_50,TimeCost_50,fobj_50,err_50,w2_y_50,b2_x_50,...
    b2_y_50,kkterr_50]  = PIDSplit_plus(gn, TF, H, HT, bg,...
    beta, gamma, ...
    eta, grad1, grad2, div, NIT, tol, verbose, benchsol);
save([dir '/' 'data_' im '_PID_50']);
clear u_50 w3_50 TimeCost_50 fobj_50 err_50 w2_y_50 b2_x_50...
    b2_y_50 kkterr_50; 

gamma = 5/beta;
[u_5,w3_5,TimeCost_5,fobj_5,err_5,w2_y_5,b2_x_5,b2_y_5, ...
    kkterr_5]=PIDSplit_plus(gn, TF, H, HT, bg, beta, gamma, ...
    eta, grad1,grad2, div, NIT, tol, verbose, benchsol);
save([dir '/' 'data_' im '_PID_5']);
clear u_5 w3_5 TimeCost_5 fobj_5 err_5 w2_y_5 b2_x_5 ...
    b2_y_5 kkterr_5; 

gamma = 1/beta;
[u_1,w3_1,TimeCost_1,fobj_1,err_1,w2_y_1,b2_x_1,b2_y_1,...
    kkterr_1]=PIDSplit_plus(gn, TF, H, HT, bg, beta, gamma, ...
    eta, grad1,grad2, div, NIT, tol, verbose, benchsol);
save([dir '/' 'data_' im '_PID_1']);
clear u_1 w3_1 TimeCost_1 fobj_1 err_1 w2_y_1 b2_x_1 b2_y_1 ...
    kkterr_1; 

gamma = 0.5/beta;
[u_05,w3_05,TimeCost_05,fobj_05,err_05,w2_y_05,b2_x_05, ...
    b2_y_05,kkterr_05]  = PIDSplit_plus(gn, TF, H, HT, bg, ...
    beta, gamma,eta, grad1, grad2, div, NIT, tol,  ...
    verbose, benchsol);
save([dir '/' 'data_' im '_PID_05']);

fprintf('Time elapsed: %g seconds',toc);
%exit;
