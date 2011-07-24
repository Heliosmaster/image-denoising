for i=[11,13,22,48,49,50]
    i = int2str(i);
    a = ['PRISMA/data_Image' i '_PID_5_0.02.mat'];
    load(a);
    
    figure(1);
    imshow(u_5,[],'Border','tight');
    h = ['Image' i ' _PID_5_0.02.eps'];
    print(1,'-depsc','-r300',h);
    close all;
   
   b  = ['PRISMA/data_Image' i '_0.02.mat'];
   load(b);
    
    figure(2);
    imshow(u,[],'Border','tight');
    k = ['Image' i ' _0.02.eps'];
        print(2,'-depsc','-r300',k);
    close all;
end
