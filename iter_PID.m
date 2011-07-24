function iter_PID(dir,im,beta,id,NIT)
% iter_PID(dir,im,beta,id,NIT)
%
% PIDSplit+ for the given image im in the given directory dir
% with the given beta, with NIT iterations
% accepts '50' '5' '1' '05' as id.

load([dir '/' im]);
load([dir '/' im '_noisy']);
load([dir '/' 'sol' im]);
fprintf('Working on %s with PID %s\n',im,id);

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

benchsol={obj,x};
TF = ones(size(obj));

if strcmp(id,'50')
    gamma = 50/beta;
    [u_50,w3_50,TimeCost_50,fobj_50,err_50,w2_y_50,b2_x_50,...
        b2_y_50,kkterr_50]=PIDSplit_plus(gn,TF,H,HT,bg,beta,...
        gamma,eta,grad1,grad2,div,NIT,tol,verbose,benchsol);
    save([dir '/' 'data_'  im '_PID_' id '_it_' int2str(NIT)]);
    clear u_50 w3_50 TimeCost_50 fobj_50 err_50 w2_y_50 b2_x_50 ...
        b2_y_50 kkterr_50;
    
elseif strcmp(id,'5')
    gamma = 5/beta;
    [u_5,w3_5,TimeCost_5,fobj_5,err_5,w2_y_5,b2_x_5,b2_y_5,...
        kkterr_5]  = PIDSplit_plus(gn, TF, H, HT, bg, beta, ...
        gamma, eta, grad1, grad2, div, NIT, tol, verbose, benchsol);
    save([dir '/' 'data_'  im '_PID_' id '_it_' int2str(NIT)]);
    clear u_5 w3_5 TimeCost_5 fobj_5 err_5 w2_y_5 b2_x_5 ...
        b2_y_5 kkterr_5;
    
elseif strcmp(id,'1')
    gamma = 1/beta;
    [u_1,w3_1,TimeCost_1,fobj_1,err_1,w2_y_1,b2_x_1,b2_y_1,...
        kkterr_1]  = PIDSplit_plus(gn, TF, H, HT, bg, beta, ...
        gamma, eta, grad1, grad2, div, NIT, tol, verbose, benchsol);
    save([dir '/' 'data_'  im '_PID_' id '_it_' int2str(NIT)]);
    clear u_1 w3_1 TimeCost_1 fobj_1 err_1 w2_y_1 b2_x_1 ...
        b2_y_1 kkterr_1;
elseif strcmp(id,'05')
    gamma = 0.5/beta;
    [u_05,w3_05,TimeCost_05,fobj_05,err_05,w2_y_05,b2_x_05,...
        b2_y_05,kkterr_05]=PIDSplit_plus(gn, TF, H, HT, bg, beta,...
        gamma,eta,grad1,grad2,div,NIT,tol,verbose, benchsol);
    save([dir '/' 'data_'  im '_PID_' id '_it_' int2str(NIT)]);
    clear u_05 w3_05 TimeCost_05 fobj_05 err_05 w2_y_05 b2_x_05 ...
        b2_y_05 kkterr_05;
end
