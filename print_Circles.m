function print_Circles(im,id,it)

it = int2str(it);


if strcmp(id,'aem')
    
    load (['circles/data_' im '_it_' it]);
    
m1 = min(u(128,:));
m2 = max(u(128,:));

figure;
plot(u(128,:));
axis([-2 260 0.5*m1 1.05*m2]);
a = [im '_it_' it];
print('-depsc','-r300',a);

elseif strcmp(id,'50')
    
   load (['circles/data_' im '_PID_' id '_it_' it]); 
    
m1 = min(u_50(128,:));
m2 = max(u_50(128,:));

figure;
plot(u_50(128,:));
axis([-2 260 0.5*m1 1.05*m2]);
a = [im '_PID_' id '_it_' it];
print('-depsc','-r300',a);


elseif strcmp(id,'5')
    load (['circles/data_' im '_PID_' id '_it_' it]); 

m1 = min(u_5(128,:));
m2 = max(u_5(128,:));

figure;
plot(u_5(128,:));
axis([-2 260 0.5*m1 1.05*m2]);
a = [im '_PID_' id '_it_' it];
print('-depsc','-r300',a);


elseif strcmp(id,'1')
    load (['circles/data_' im '_PID_' id '_it_' it]); 
m1 = min(u_1(128,:));
m2 = max(u_1(128,:));

figure;
plot(u_1(128,:));
axis([-2 260 0.5*m1 1.05*m2]);
a = [im '_PID_' id '_it_' it];
print('-depsc','-r300',a);


elseif strcmp(id,'05')
    load (['circles/data_' im '_PID_' id '_it_' it]); 

m1 = min(u_05(128,:));
m2 = max(u_05(128,:));

figure(1);
plot(u_05(128,:));
axis([-2 260 0.5*m1 1.05*m2]);
a = [im '_PID_' id '_it_' it];
print('-depsc','-r300',a);


end