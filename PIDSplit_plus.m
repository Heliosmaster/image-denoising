function [u,w3,TimeCost,fobj,err,varargout] = ...
    PIDSplit_plus(f, TF, K, KT, bg, beta, gamma, eta, grad1,...
    grad2, div, NIT, tol, verbose, obj)

%Input parameters
%
%f              data
%TF             Fourier transform of the PSF
%K, KT          function handles to the matrix-vector products 
%               by K and K^T
%bg             background constant term
%beta           regularization parameter
%eta            lower bound
%grad1, grad2   function handles to the gradient with respect to
%               the x and y directions
%div            function handle to the negative divergence operator
%NIT            maximum number of iterations
%tol            tolerance on the distance between two
%               successive iterates
%verbose        flag (=1 print information at each iteration;
%                    =0 no printing)
%obj            cell array containing reference images

%Output parameters
%
%u,w3           computed reconstructions (u may have negative
%               components)
%TimeCost       time per iteration (seconds)
%fobj           objective function value per iteration
%err            cell array containing relative errors per 
%               iteration with respect
%               to each of the reference images in obj

% Code by Silvia Bonettini
% Edited by Davide Taviani


nobj = length(obj);
for i = 1:nobj
    normobj(i) = norm(obj{i}(:));
    err{i} = zeros(NIT+1,1); %arrays to store errors per iteration
end

fobj = zeros(NIT+1,1);
TimeCost = zeros(NIT+1,1);
kkterr=zeros(NIT+1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n=length(f);                %Assume a square image
zeroindex = f <= 0;
nonzeroindex = ~zeroindex;

rapp = zeros(size(f));

%%%Eigenvalues of the Laplacian operator
P = zeros(n);
P(n/2,n/2+1) = 1; P(n/2+1,n/2+1) = -1;
PTF = fft2(fftshift(P));
PTF2 = abs(PTF).^2;
Q   = P';
QTF = fft2(fftshift(Q));
QTF2 = abs(QTF).^2;
TF2  = abs(TF).^2;


%%% Initial points %%
b1 = zeros(n);
b2_x = zeros(n);
b2_y = zeros(n);
b3 = zeros(n);
w1 = K(f);
w2_x = grad1(f);%[f(:,2:n)-f(:,1:n-1), f(:,1)-f(:,n)];
w2_y = grad2(f);%[f(2:n,:)-f(1:n-1,:); f(1,:)-f(n,:)];
w3  = f;
u = f;

Ku  = K(w3);
den = Ku + bg;
rapp = zeros(n);
rapp(nonzeroindex) = f(nonzeroindex)./den(nonzeroindex);
KL = sum( f(nonzeroindex).* log(rapp(nonzeroindex)) + ...
    den(nonzeroindex) - f(nonzeroindex) );
KL = KL + sum(den(zeroindex));
gu_norm = sqrt(w2_x.^2+w2_y.^2);
fobj(1)   = KL + beta*sum(sum(gu_norm));

for i=1:nobj
    err{i}(1) = norm(obj{i}(:)-w3(:))/normobj(i);
end
if verbose
    fprintf('\nInitial: fobj=%8.3e', fobj(1));
    for i=1:nobj
        fprintf(' err{%g} %g', i, err{i}(1));
    end
end

TimeCost(1)=0;
t0 = cputime;                %Start CPU clock

for itr=1:NIT
    % Save old values
    uold = u;
    w1old = w1;
    w2_xold = w2_x;
    w2_yold = w2_y;
    w3old = w3;
    
    b1old = b1;
    
    b2_xold = b2_x;
    b2_yold = b2_y;
    
    b3old = b3;
    
    %Compute u^{(k+1)}
    t1 = KT( w1 - b1 );
    x1 = w2_x - b2_x; x2 = w2_y - b2_y;
    t2 = div(x1,x2);%([-x1(:,1)+x1(:,n),-x1(:,2:n)+ ...
                    %x1(:,1:n-1)]+[-x2(1,:)+x2(n,:);...
                    %-x2(2:n,:)+x2(1:n-1,:)]);
    t3 = w3 - b3;
    t = t1 + t2 + t3;
    
    u = ifft2(fft2(t1+t2+t3)./(1+TF2+PTF2+QTF2));
    
    
    %Update w1
    Ku = K(u) + bg;
    w1 = 0.5*( b1 + Ku - gamma + sqrt( ( b1 + Ku - ...
        gamma ).^2 + 4*gamma*f ) );
    
    %Update w2
    Du_x  = grad1(u);%[u(:,2:n)-u(:,1:n-1), u(:,1)-u(:,n)];
    Du_y  = grad2(u);%[u(2:n,:)-u(1:n-1,:); u(1,:)-u(n,:)];
    w2_x = b2_x + Du_x;
    w2_y = b2_y + Du_y;
    % shrink
    norm_w2 = sqrt(w2_x.^2 + w2_y.^2);
    ij  = norm_w2 >= gamma*beta;
    nij = ~ij;
    w2_x(ij) = w2_x(ij) - gamma*beta*w2_x(ij)./norm_w2(ij);
    w2_x(nij)  = 0;
    w2_y(ij) = w2_y(ij) - gamma*beta*w2_y(ij)./norm_w2(ij);
    w2_y(nij)  = 0;
    
    %update w3
    w3 = max(b3 + u, eta);
    
    %update b
    b1   = b1 + Ku - w1;
    b2_x = b2_x + Du_x - w2_x;
    b2_y = b2_y + Du_y - w2_y;
    b3   = b3 + u - w3;
    
    Kw3 = K(w3);
    den = Kw3 + bg;
    rapp(nonzeroindex) = f(nonzeroindex)./den(nonzeroindex);
    KL = sum( f(nonzeroindex).* log(rapp(nonzeroindex)) + ...
        den(nonzeroindex) - f(nonzeroindex) );
    KL = KL + sum(den(zeroindex));
    Dw3_x  = grad1(w3);%[w3(:,2:n)-w3(:,1:n-1), w3(:,1)-w3(:,n)];
    Dw3_y  = grad2(w3);%[w3(2:n,:)-w3(1:n-1,:); w3(1,:)-w3(n,:)];
    
    gu_norm = sqrt(Dw3_x.^2+Dw3_y.^2);
    fobj(itr+1)   = KL + beta*sum(sum(gu_norm));
    TimeCost(itr+1) = cputime - t0;
    
    for i=1:nobj
        err{i}(itr + 1) = norm(obj{i}(:)-w3(:))/normobj(i);
    end
    
    if verbose
        fprintf('\n%4d): f(x)=%g ', itr, ...
            fobj(itr+1) );
        for i=1:nobj
            fprintf(' err{%g} %g', i, err{i}(itr + 1));
        end
    end
    w = [u(:);w3(:);w2_x(:);w2_y(:);w1(:);b3(:);b2_x(:);...
        b2_y(:);b1(:)];
    wold = [uold(:);w3old(:);w2_xold(:);w2_yold(:);w1old(:);...
        b3old(:);b2_xold(:);b2_yold(:);b1old(:)];
    kkterr(itr+1) = norm(w(:)-wold(:))/norm(w(:));
    %kkterr(itr+1)=norm(w3(:)-w3old(:))/norm(w3(:));
    if kkterr(itr+1)< tol
        break
    end
end
% end of the main loop
fobj(itr +2:end) = [];
for i = 1:nobj
    err{i}(itr + 2:end) = [];
end
TimeCost(itr+2:end) = [];
kkterr(itr+2:end)=[];
if nargout > 0
    varargout{1} = w2_x;
    varargout{2} = w2_y;
    varargout{3} = b2_x;
    varargout{4} = b2_y;
    varargout{5} =kkterr;
end