function load_prisma(id,i,arrbeta)
% load_prisma(id,i,arrbeta)
% function used to print values to form the relative table

i = int2str(i);

for beta = arrbeta;
    
    fprintf('beta: %g\n',beta);
    
    if strcmp(id,'aem')
        load(['PRISMA/' 'data_Image' i '_' num2str(beta) '.mat' ]);
        fprintf('AEM: time %g, iter %g, err: %g, f:%g\n',...
            TimeCost(end),length(err{1})-1,err{1}(end),Primal(end));
    elseif strcmp(id,'50')
        load(['PRISMA/' 'data_Image' i '_PID_' id '_' ...
            num2str(beta) '.mat']);
        fprintf('PID%s: time %g, iter %g, err: %g, f:%g\n',id,...
            TimeCost_50(end),length(err_50{1})-1,err_50{1}(end),...
            fobj_50(end));
    elseif strcmp(id,'5')
        load(['PRISMA/' 'data_Image' i '_PID_' id '_' ...
            num2str(beta) '.mat']);
        fprintf('PID%s: time %g, iter %g, err: %g, f:%g\n',id,...
            TimeCost_5(end),length(err_5{1})-1,err_5{1}(end),...
            fobj_5(end));
        
    elseif strcmp(id,'1')
        
        load(['PRISMA/' 'data_Image' i '_PID_' id '_' ...
            num2str(beta) '.mat']);
        fprintf('PID%s: time %g, iter %g, err: %g, f:%g\n',id,...
            TimeCost_1(end),length(err_1{1})-1,err_1{1}(end),...
            fobj_1(end));
        
    elseif strcmp(id,'05')
        load(['PRISMA/' 'data_Image' i '_PID_' id '_' ...
            num2str(beta) '.mat']);
        fprintf('PID%s: time %g, iter %g, err: %g, f:%g\n',id,...
            TimeCost_05(end),length(err_05{1})-1,err_05{1}(end),...
            fobj_05(end));
        
    end
end
%exit;