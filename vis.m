function vis(dir,im)

load([dir '/' 'data_' im]);
%clear x y1 y2 y3 TimeCost InnVec Primal phi_xy alpha vec err KKT KL; 


load([dir '/' 'data_' im '_PID_50']);

%clear u_50 w3_50 TimeCost_50 fobj_50 err_50 w2_y_50 b2_x_50 b2_y_50 kkterr_50; 

load([dir '/' 'data_' im '_PID_5']);
%clear u_5 w3_5 TimeCost_5 fobj_5 err_5 w2_y_5 b2_x_5 b2_y_5 kkterr_5; 


load([dir '/' 'data_' im '_PID_1']);
load([dir '/' 'data_' im '_PID_05']);

n = length(err{1});

%figure(1);
%title(['Restoration error for ' im]);
%hold on;
%plot(1:n,err{1});
%plot(1:n,err_50{1},'g');
%plot(1:n,err_5{1},'r');
%plot(1:n,err_1{1},'c');
%plot(1:n,err_05{1},'k');
%legend('AEM','PID50','PID5','PID1','PID05');

figure(2);
m2 = max([max(err{2}),max(err_50{2}),max(err_5{2}),max(err_1{2}),max(err_05{2})]);

t2 = max([max(TimeCost),max(TimeCost_50),max(TimeCost_5),max(TimeCost_1),max(TimeCost_05)]);

semilogy(TimeCost,err{2},'--');
hold on;
semilogy(TimeCost_50,err_50{2},'g');
semilogy(TimeCost_5,err_5{2},'r-.');
semilogy(TimeCost_1,err_1{2},'k');
semilogy(TimeCost_05,err_05{2},'k:');
axis([0 t2 0 m2]);
%axis tight
legend('AEM','PID50','PID5','PID1','PID05');
% title(['Error from ideal solution for ' im]);

%exit;

% l=5;
% figure(3);
% %title(['Error for first ' int2str(l) ' seconds on ' im]);
% hold on;
% xlim = ([0 5]);
% i = find(TimeCost>l,1,'first');
% i_50 = find(TimeCost_50>l,1,'first');
% i_5 = find(TimeCost_5>l,1,'first');
% i_1 = find(TimeCost_1>l,1,'first');
% i_05 = find(TimeCost_05>l,1,'first');
% 
% m3 = max([max(err{1}(1:i)),max(err_50{1}(1:i_50)),max(err_5{1}(1:i_5)),max(err_1{1}(1:i_1)),max(err_05{1}(1:i_05))]);
% m4 = min([min(err{1}(1:i)),min(err_50{1}(1:i_50)),min(err_5{1}(1:i_5)),min(err_1{1}(1:i_1)),min(err_05{1}(1:i_05))]);
% 
% plot(TimeCost(1:i),err{1}(1:i),'--');
% plot(TimeCost_50(1:i_50),err_50{1}(1:i_50),'g');
% plot(TimeCost_5(1:i_5),err_5{1}(1:i_5),'r-.');
% plot(TimeCost_1(1:i_1),err_1{1}(1:i_1),'k');
% plot(TimeCost_05(1:i_05),err_05{1}(1:i_05),'k:');
% axis([0 l 0.9*m4 m3]);
% %axis tight
% legend('AEM','PID50','PID5','PID1','PID05');

% a = [im '_err_' int2str(l) '_sec'];
b = [im '_err2'];

% print(3,'-depsc','-r300',a);
print(2,'-depsc','-r300',b);

%fprintf('for %s,\n',im);
%fprintf('AEM: iter %g, err %g\n',i,err{1}(i));
%fprintf('PID50: iter %g, err %g\n',i_50,err_50{1}(i_50));
%fprintf('PID5: iter %g, err %g\n',i_5,err_5{1}(i_5));
%fprintf('PID1: iter %g, err %g\n',i_1,err_1{1}(i_1));
%fprintf('PID05: iter %g, err %g\n',i_05,err_05{1}(i_05));