clear all


Image=imread('D:\160613_SFR test for C4\Previous\15_ROI.bmp');
%%
N=23;

Image_Blur=conv2(Image,ones(N,N)/N^2);

imagesc(Image_Blur);
%%
IntX=sum(Image_Blur,2);
IntY=sum(Image_Blur,1);

plot(IntX);