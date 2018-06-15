clear all

Pixel_Size=0.45;    %micron

Image=double(imread('5.jpg'))';
%%

imagesc(Image);
colormap(gray);
axis equal
axis off
xlim([0 size(Image,2)]);
ylim([0 size(Image,1)]);

Center=[332 8];

Radius=486;

X=repmat([1:648]',[1 488]);
Y=repmat([1:488],[648 1]);

Mask=ones(648,488);
Mask(((X-Center(1)).^2+(Y-Center(2)).^2)>Radius^2)=0;

Masked_Image=Image.*Mask;
imagesc(Masked_Image);
colormap(gray);
axis equal
axis off
xlim([0 size(Image,2)]);
ylim([0 size(Image,1)]);

Diameter=Radius*2*Pixel_Size