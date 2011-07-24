for i=[11,13,22,48,49,50];
    im = ['Image' int2str(i) '.tif'];
    image = 4095-double(imread(im));
    image = image(:,138:404);
    figure(i);
    imshow(image,[],'Border','tight');
    print(i,'-depsc','-r300',int2str(i))
 %   pause;
end

close all;

im ='Avg.tif';
image = 4095-double(imread(im));
image = image(:,138:404);
figure(1);
imshow(image,[],'Border','tight');
print(1,'-depsc','-r300','avg');
