function vis_iter(dir,im,id,it)

it = int2str(it);

if strcmp(id,'aem')
    
    load ([dir '/data_' im '_it_' it]);
   
    fprintf('AEM: it %s, err: %g\n',it,err{1}(end));

elseif strcmp(id,'50')
        load ([dir '/data_' im '_PID_' id '_it_' it]); 
fprintf('PID%s: it %s, err: %g\n',id,it,err_50{1}(end));
elseif strcmp(id,'5')
    load ([dir '/data_' im '_PID_' id '_it_' it]); 
    fprintf('PID%s: it %s, err: %g\n',id,it,err_5{1}(end));

elseif strcmp(id,'1')
    load ([dir '/data_' im '_PID_' id '_it_' it]); 
fprintf('PID%s: it %s, err: %g\n',id,it,err_1{1}(end));

elseif strcmp(id,'05')
    load ([dir '/data_' im '_PID_' id '_it_' it]); 
fprintf('PID%s: it %s, err: %g\n',id,it,err_05{1}(end));
end