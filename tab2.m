clear all;
close all;

dir = 'circles';
im = 'circles';


for thr=[1e-2,1e-4,1e-5]
    fprintf('\n');
    disp(thr);
    for l=[3,4,5,6,7,2,13,14,15,16,17]
        a = [dir '/' 'data_' int2str(l) '_' im];
        load(a);
        i = find(err{2}<thr,1,'first');
        j = TimeCost(i);
        fprintf('Im%g %g\n',l,j);
    end
   
end
    