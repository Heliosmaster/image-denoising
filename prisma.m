function prisma(id,i,arrbeta)
%prisma(id,i,arrbeta)
% restoration of the image i from project PRISMA
% with the given array of beta arrbeta
% using both AEM and PIDSplit+
%
% id accepts 'aem','50','5','1','05'.

i = int2str(i);

image = ['PRISMA/Image' i '.tif'];
gn = 4095-double(imread(image));
gn = gn(1:266,138:403);

obj = 4095-double(imread('PRISMA/Avg.tif'));
obj = obj(1:266,138:403);

n = size(gn,1);

grad1  = @(y)[y(:,2:n)-y(:,1:n-1), y(:,1)-y(:,n)];
grad2  = @(y)[y(2:n,:)-y(1:n-1,:); y(1,:)-y(n,:)];
div    = @(x1,x2)([-x1(:,1)+x1(:,n),-x1(:,2:n)+x1(:,1:n-1)] +...
    [-x2(1,:)+x2(n,:);-x2(2:n,:)+x2(1:n-1,:)]);
H = @(x)x;
HT = @(x)x;

bg = 0;
delta = 0;
NIT = 3000;
tol = 1e-6;

zeroindex = gn <= 0;
nonzeroindex = ~zeroindex;
eta = min(gn(nonzeroindex))*ones(size(gn));
eta(zeroindex) = 0;

verbose =false;
Nalpha=10;
epsilon = 1e-4;

benchsol={obj};
TF = ones(size(obj));

for beta = arrbeta;
    
    if strcmp(id,'aem')
        
        fprintf('AEM: Working on Image%s with beta %g\n', i,beta);
        [u,y1,y2,y3,TimeCost,InnVec,Primal,phi_xy,alpha_vec,err,...
            KKT,KL]= AEM(gn, H, HT,bg, beta, delta, eta, grad1,...
            grad2, div, NIT, tol, verbose, Nalpha, epsilon,...
            benchsol);
        fprintf('Number of iteration performed:%g\n',...
            length(err{1})-1);
        save(['PRISMA/' 'data_Image' i '_'...
            num2str(beta) '.mat' ]);
        
    elseif strcmp(id,'50')
        
        gamma = 50/beta;
        fprintf('PID%s: Working on Image%s with beta %g\n',...
            id,i,beta);
        [u_50,w3_50,TimeCost_50,fobj_50,err_50,w2_y_50,...
            b2_x_50,b2_y_50,kkterr_50]  = PIDSplit_plus(gn, TF,...
            H,HT,bg, beta, gamma, eta, grad1, grad2, div, NIT,...
            tol,verbose, benchsol);
        fprintf('Number of iteration performed:%g\n',...
            length(err_50{1})-1);
        save(['PRISMA/' 'data_Image' i '_PID_' id '_' ...
            num2str(beta) '.mat']);
        
    elseif strcmp(id,'5')
        gamma = 5/beta;
        fprintf('PID%s: Working on Image%s with beta %g\n',...
            id,i,beta);
        [u_5,w3_5,TimeCost_5,fobj_5,err_5,w2_y_5,b2_x_5,...
            b2_y_5,kkterr_5]  = PIDSplit_plus(gn, TF, H, HT,...
            bg, beta, gamma, eta, grad1, grad2, div, NIT, ...
            tol, verbose, benchsol);
        fprintf('Number of iteration performed:%g\n',...
            length(err_5{1})-1);
        save(['PRISMA/' 'data_Image' i '_PID_' id '_' ...
            num2str(beta) '.mat']);
        
    elseif strcmp(id,'1')
        gamma = 1/beta;
        fprintf('PID%s: Working on Image%s with beta %g\n',...
            id,i,beta);
        [u_1,w3_1,TimeCost_1,fobj_1,err_1,w2_y_1,b2_x_1,...
            b2_y_1,kkterr_1]  = PIDSplit_plus(gn, TF, H, HT,...
            bg, beta, gamma, eta, grad1, grad2, div, NIT, ...
            tol, verbose, benchsol);
        fprintf('Number of iteration performed:%g\n',...
            length(err_1{1})-1);
        save(['PRISMA/' 'data_Image' i '_PID_' id '_' ...
            num2str(beta) '.mat']);
        
    elseif strcmp(id,'05')
        gamma = 0.5/beta;
        fprintf('PID%s: Working on Image%s with beta %g\n',...
            id,i,beta);
        [u_05,w3_05,TimeCost_05,fobj_05,err_05,w2_y_05,...
            b2_x_05,b2_y_05,kkterr_05]  = PIDSplit_plus(gn,...
            TF, H, HT, bg, beta, gamma, eta, grad1, grad2, ...
            div, NIT, tol, verbose, benchsol);
        fprintf('Number of iteration performed:%g\n',...
            length(err_05{1})-1);
        save(['PRISMA/' 'data_Image' i '_PID_' id '_' ...
            num2str(beta) '.mat']);
        
    end
end
%exit;