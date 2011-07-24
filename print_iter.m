function print_iter(dir,im,id,it)

it = int2str(it);


if strcmp(id,'aem')
    
load ([dir '/data_' im '_it_' it]);
figure;
imshow(u,[],'Border','tight');
a = [im '_it_' it];
print('-depsc','-r300',a);
  
elseif strcmp(id,'50')
    
   load ([dir '/data_' im '_PID_' id '_it_' it]); 
    figure;
imshow(u_50,[],'Border','tight');
a = [im '_PID_' id '_it_' it];
print('-depsc','-r300',a);

elseif strcmp(id,'5')
    load ([dir '/data_' im '_PID_' id '_it_' it]); 

figure;
imshow(u_5,[],'Border','tight');
a = [im '_PID_' id '_it_' it];
print('-depsc','-r300',a);

elseif strcmp(id,'1')
    load ([dir '/data_' im '_PID_' id '_it_' it]); 
figure;
imshow(u_1,[],'Border','tight');
a = [im '_PID_' id '_it_' it];
print('-depsc','-r300',a);

elseif strcmp(id,'05')
    load ([dir '/data_' im '_PID_' id '_it_' it]); 
figure;
imshow(u_05,[],'Border','tight');
a = [im '_PID_' id '_it_' it];
print('-depsc','-r300',a);

end