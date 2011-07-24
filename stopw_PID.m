function stopw_PID(dir,im,beta,id,it)
% stopw_PID(dir,im,beta,id,it)
%
% repeat the restoration of the given image im
% in the given directory dir, using the given beta
% IT times, in order to evaluate precisely the
% time required by PID for such image.
%
% accepts id as '50','5','1','05'

load([dir '/' im]);
load([dir '/' im '_noisy']);
load([dir '/' 'sol' im]);
a = ['TimeCost_' id];
load([dir '/' 'data_' im '_PID_' id ],a);

fprintf('Working on %s with PID %s\n',im,id);

if strcmp(id,'50')
    T1 = TimeCost_50;
elseif strcmp(id,'5')
    T1 = TimeCost_5;
elseif strcmp(id,'1')
    T1 = TimeCost_1;
elseif strcmp(id,'05')
    T1 = TimeCost_05;
end

[n,n] = size(gn);

grad1  = @(y)[y(:,2:n)-y(:,1:n-1), y(:,1)-y(:,n)];
grad2  = @(y)[y(2:n,:)-y(1:n-1,:); y(1,:)-y(n,:)];
div    = @(x1,x2)([-x1(:,1)+x1(:,n),-x1(:,2:n)+x1(:,1:n-1)] + ...
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

benchsol={obj,x};
TF = ones(size(obj));

if strcmp(id,'50')
    gamma = 50/beta;
    for s=1:it
        [u_50,w3_50,TimeCost_50,fobj_50,err_50,w2_y_50,...
            b2_x_50,b2_y_50,kkterr_50]  = PIDSplit_plus(gn,...
            TF, H, HT, bg, beta, gamma, eta, grad1, grad2,...
            div, NIT, tol, verbose, benchsol);
        fprintf('%g) Current: %g, Fastest: %g\n',s,...
            TimeCost_50(end),T1(end));
        if TimeCost_50(end) < T1(end)
            save([dir '/' 'data_'  im '_PID_' id ]);
            T1 = TimeCost_50;
        end
        clear u_50 w3_50 TimeCost_50 fobj_50 err_50 w2_y_50 ...
            b2_x_50 b2_y_50 kkterr_50;
    end
elseif strcmp(id,'5')
    gamma = 5/beta;
    for s=1:it
        [u_5,w3_5,TimeCost_5,fobj_5,err_5,w2_y_5,b2_x_5,...
            b2_y_5,kkterr_5]  = PIDSplit_plus(gn, TF, H,...
            HT, bg, beta, gamma, eta, grad1, grad2, div,...
            NIT, tol, verbose, benchsol);
        fprintf('%g) Current: %g, Fastest: %g\n',s,...
            TimeCost_5(end),T1(end));
        if TimeCost_5(end) < T1(end)
            save([dir '/' 'data_'  im '_PID_' id]);
            T1 = TimeCost_5;
        end
        clear u_5 w3_5 TimeCost_5 fobj_5 err_5 w2_y_5 ...
            b2_x_5 b2_y_5 kkterr_5;
    end
elseif strcmp(id,'1')
    gamma = 1/beta;
    for s=1:it
        [u_1,w3_1,TimeCost_1,fobj_1,err_1,w2_y_1,b2_x_1,...
            b2_y_1,kkterr_1]  = PIDSplit_plus(gn, TF, H, ...
            HT, bg, beta, gamma, eta, grad1, grad2, div, ...
            NIT, tol, verbose, benchsol);
        fprintf('%g) Current: %g, Fastest: %g\n',s,...
            TimeCost_1(end),T1(end));
        if TimeCost_1(end) < T1(end)
            save([dir '/' 'data_'  im '_PID_' id]);
            T1 = TimeCost_1;
        end
        clear u_1 w3_1 TimeCost_1 fobj_1 err_1 w2_y_1 b2_x_1 ...
            b2_y_1 kkterr_1;
    end
elseif strcmp(id,'05')
    gamma = 0.5/beta;
    for s=1:it
        [u_05,w3_05,TimeCost_05,fobj_05,err_05,w2_y_05,...
            b2_x_05,b2_y_05,kkterr_05]  = PIDSplit_plus(gn, ...
            TF, H, HT, bg, beta, gamma, eta, grad1, grad2, ...
            div, NIT, tol, verbose, benchsol);
        fprintf('%g) Current: %g, Fastest: %g\n',s,...
            TimeCost_05(end),T1(end));
        if TimeCost_05(end) < T1(end)
            save([dir '/' 'data_'  im '_PID_' id]);
            T1 = TimeCost_05;
        end
        clear u_05 w3_05 TimeCost_05 fobj_05 err_05 w2_y_05 ...
            b2_x_05 b2_y_05 kkterr_05;
    end
end
