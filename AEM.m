function [x,y1,y2,y3,TimeCost,InnVec,Primal,phi_xy,alpha_vec,err,...
    varargout] = ...
    AEM(z, H, HT,bg, beta, delta, eta, grad1, grad2, div, NIT, ...
    tol, verbose, Nalpha, epsilon, obj)

%Nalpha = -1 per non tenere memoria%
%epsilon = 1e-4;%5e-1; %provare
% obj is a cell array of images to compute the relative difference

nobj = length(obj);
for i = 1:nobj
    normobj(i) = norm(obj{i}(:));
    err{i}     = zeros(NIT+1,1); %arrays to store errors per iteration
end

Primal = zeros(NIT+1,1);
phi_xy = zeros(NIT+1,1);
alpha_vec = zeros(NIT+1,1);
InnVec = zeros(NIT+1,1);
TimeCost = zeros(NIT+1,1);
kkt_vec  = zeros(NIT+1,1);
KL_vec   = zeros(NIT+1,1);

gamma   = 0.99;
theta   = 0.99;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


totred = 0;
n=length(z);                %Assume a square image
zeroindex = z <= 0;
nonzeroindex = ~zeroindex;
rapp = zeros(size(z));

%%% Initial point %%
x  = max(eta,z);
y1 = zeros(size(z));
y2 = y1;
y3 = zeros(size(z));

%%% Projection of the initial x
ynorm= max(1, sqrt(y1.^2+y2.^2+y3.^2));
y1 = y1./ynorm;
y2 = y2./ynorm;
y3 = y3./ynorm;
% periodic boundary conditions
dx1 = grad1(x);
dx2 = grad2(x);

gu_norm = sqrt(dx1.^2+dx2.^2+delta^2); %total variarion
den = H(x) + bg;%ifft2(fft2(y).*TF)+bg;
rapp(nonzeroindex) = z(nonzeroindex)./ den(nonzeroindex);
KL = sum( z(nonzeroindex).* log(rapp(nonzeroindex)) + ...
    den(nonzeroindex) - z(nonzeroindex) );
KL = KL + sum(den(zeroindex));
g_KL = 1 - HT(rapp);%1 - ifft2(conj(TF).*fft2(rapp));

Ay = div(y1,y2);
g = g_KL + beta*Ay;

Primal(1) = beta*sum(sum(gu_norm)) + KL;
phi_xy(1) = beta*sum(sum(Ay.*x))+ beta*delta* sum(sum(y3)) +KL;
KL_vec(1) = KL;
for i=1:nobj
    err{i}(1) = norm(obj{i}(:)-x(:))/normobj(i);
end

if verbose
    fprintf('\nInitial: Primal=%8.3e', Primal(1));
    for i=1:nobj
        fprintf(' err{%g} %g', i, err{i}(1));
    end
    
end
alpha = 1;
alphaOK(1) = alpha;
alpha_vec(1) = alpha;
TimeCost(1)=0;
t0 = cputime;                %Start CPU clock

for itr=1:NIT
    
    y1old = y1; y2old = y2; y3old = y3;
    
    y1 = y1 + alpha*beta*dx1;
    y2 = y2 + alpha*beta*dx2;
    y3 = y3 + alpha*beta*delta;
    % Projection on the set X={||x||<=1}
    ynorm= max(1, sqrt(y1.^2+y2.^2+y3.^2));
    y1 = y1./ynorm;
    y2 = y2./ynorm;
    y3 = y3./ynorm;
    Ay = div(y1,y2);
    
    xold = x;
    g_KLold = g_KL;
    dx1old = dx1; dx2old = dx2;
    
    x = x - alpha*( g_KL + beta*Ay );
    x = max(x,eta);
    
    dx1 = grad1(x);
    dx2 = grad2(x);
    
    den  = H(x) + bg;
    rapp(nonzeroindex) = z(nonzeroindex)./ den(nonzeroindex);
    g_KL = 1-HT(rapp);
    
    Deltax = x - xold; normDeltax = norm(Deltax(:));
    Ak = norm(g_KLold(:) - g_KL(:)) / normDeltax;
    Bk = beta*sqrt( norm(dx1(:)-dx1old(:))^2 + norm(dx2(:)-...
        dx2old(:))^2 )/normDeltax;
    cond = 1-2*alpha*Ak-2*alpha^2*Bk^2;
    %fprintf('\n L1=%g L3=%g cond %g\n',L1,L3,cond);
    alphabar = gamma*(sqrt(Ak^2+2*Bk^2*(1-epsilon))-Ak)/(2*Bk^2);
    ired = 0;
    while cond < epsilon
        alpha = min(alphabar,theta*alpha);
        y1 = y1old + alpha*beta*dx1old;
        y2 = y2old + alpha*beta*dx2old;
        y3 = y3old + alpha*beta*delta;
        % Projection on the set X={||x||<=1}
        ynorm= max(1, sqrt(y1.^2+y2.^2+y3.^2));
        y1 = y1./ynorm;
        y2 = y2./ynorm;
        y3 = y3./ynorm;
        Ay = div(y1,y2);
        
        x = xold - alpha*( g_KLold + beta*Ay );
        x = max(x,eta);
        
        dx1 = grad1(x);
        dx2 = grad2(x);
        
        den  = H(x)+bg;
        rapp(nonzeroindex) = z(nonzeroindex)./ den(nonzeroindex);
        g_KL = 1-HT(rapp);
        
        Deltax = x - xold; normDeltax = norm(Deltax(:));
        Ak = norm(g_KLold(:) - g_KL(:)) / normDeltax;
        Bk = beta*sqrt( norm(dx1(:)-dx1old(:))^2 + norm(dx2(:)-...
            dx2old(:))^2 )/normDeltax;
        cond = 1-2*alpha*Ak-2*alpha^2*Bk^2;
        alphabar = gamma*(sqrt(Ak^2+2*Bk^2*(1-epsilon))-...
            Ak)/(2*Bk^2);
        ired = ired + 1;
        %            return
    end
    alpha_vec(itr + 1) = alpha;
    totred = totred + ired;
    InnVec(itr + 1) = ired;
    y1 = y1old + alpha*beta*dx1;
    y2 = y2old + alpha*beta*dx2;
    y3 = y3old + alpha*beta*delta;
    % Projection on the set X={||x||<=1}
    ynorm = max(1, sqrt(y1.^2+y2.^2+y3.^2));
    y1 = y1./ynorm;
    y2 = y2./ynorm;
    y3 = y3./ynorm;
    Ay = div(y1,y2);
    KL = sum( z(nonzeroindex).* log(rapp(nonzeroindex)) + ...
        den(nonzeroindex) - z(nonzeroindex) );
    KL = KL + sum(den(zeroindex));
    
    gu_norm = sqrt(dx1.^2+dx2.^2+delta^2);
    alphaOK(itr) = alpha;
    alpha = mean([alphaOK(max(1,itr-Nalpha):itr) alphabar]);
    
    
    Primal(itr+1) = beta*sum(sum(gu_norm)) + KL;
    phi_xy(itr+1) = beta*sum(sum(Ay.*x))+ beta*delta* ...
        sum(sum(y3)) + KL ;
    KL_vec(itr+1) = KL;
    TimeCost(itr+1)=cputime-t0;
    
    for i=1:nobj
        err{i}(itr + 1) = norm(obj{i}(:)-x(:))/normobj(i);
    end
    
    if verbose
        fprintf('\n%4d): f(x)=%g Phi(x,y)=%g alpha %g', itr, ...
            Primal(itr+1), phi_xy(itr+1), alpha );
        for i=1:nobj
            fprintf(' err{%g} %g', i, err{i}(itr + 1));
        end
    end
    normDeltay =sqrt(sum((y1(:)-y1old(:)).^2) + sum((y2(:)-...
        y2old(:)).^2) + sum((y3(:)-y3old(:)).^2));
    normy = sqrt(sum(y1(:).^2) + sum(y2(:).^2) + sum(y3(:).^2));
    kkt_vec(itr+1) = sqrt(normDeltax^2+normDeltay^2)/...
        sqrt(norm(x(:))^2+normy^2);
    if  kkt_vec(itr+1) < tol
        break;
    end
end
%end of the main loop
Primal(itr+2:end) = [];
for i = 1:nobj
    err{i}(itr + 2:end) = [];
end
phi_xy(itr+2:end) = [];
alpha_vec(itr+2:end) = [];
InnVec(itr+1:end) = [];
TimeCost(itr+2:end) = [];
if nargout >= 11
    varargout{1} = kkt_vec(1:itr+1);
    if nargout == 12
        varargout{2} = KL_vec(1:itr +1);
    end
end
if verbose
    fprintf('\n');
end
