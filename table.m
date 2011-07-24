clear all;
close all;

dir = 'circles';
im = 'circles_d5';

load([dir '/' 'data_' im]);
load([dir '/' 'data_' im '_PID_50']);
load([dir '/' 'data_' im '_PID_5']);
load([dir '/' 'data_' im '_PID_1']);
load([dir '/' 'data_' im '_PID_05']);

disp(im);

for tol=[1e-2,1e-4,1e-5,5e-6]
    disp(tol);
    i = find(err{2}<tol,1,'first');
    if isempty(i)
        i = 0;
        j = 0;
        k = 0;
    else
        j = TimeCost(i);
        k = err{1}(i);
    end
    fprintf('AEM: iter %g, time %g, err %g \n',i,j,k);
    
    i = find(err_50{2}<tol,1,'first');
    if isempty(i)
        i = 0;
        j = 0;
        k = 0;
    else
        j = TimeCost_50(i);
        k = err_50{1}(i);
    end
    fprintf('PID50: iter %g, time %g,err %g\n',i,j,k);
    
     i = find(err_5{2}<tol,1,'first');
    if isempty(i)
        i = 0;
        j = 0;
        k = 0;
    else
        j = TimeCost_5(i);
        k = err_5{1}(i);
    end
    fprintf('PID5: iter %g, time %g, err %g\n',i,j,k);
    
     i = find(err_1{2}<tol,1,'first');
    if isempty(i)
        i = 0;
        j = 0;
        k = 0;
    else
        j = TimeCost_1(i);
        k = err_1{1}(i);
    end
    fprintf('PID1: iter %g, time %g, err %g\n',i,j,k);
    
     i = find(err_05{2}<tol,1,'first');
    if isempty(i)
        i = 0;
        j = 0;
        k = 0;
    else
        j = TimeCost_05(i);
        k = err_05{1}(i);
    end
    fprintf('PID05: iter %g, time %g, err %g\n\n',i,j,k);
    
end

fprintf('After %g iterations\n',length(err{1})-1);
fprintf('AEM: time %g, e %g, err %g\n',TimeCost(end),err{2}(end),err{1}(end));
fprintf('PID50: time %g, e %g, err %g\n',TimeCost_50(end),err_50{2}(end),err_50{1}(end));
fprintf('PID5: time %g, e %g, err %g\n',TimeCost_5(end),err_5{2}(end),err_5{1}(end));
fprintf('PID1: time %g, e %g, err %g\n',TimeCost_1(end),err_1{2}(end),err_1{1}(end));
fprintf('PID05: time %g, e %g, err %g\n',TimeCost_05(end),err_05{2}(end),err_05{1}(end));
    