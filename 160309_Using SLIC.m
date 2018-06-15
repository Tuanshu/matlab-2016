clear all

cd('D:\MATLAB\SLIC\');
frame_width=7501;
frame_height=3757;

folder_path_without_index='G:\Everday Experiements\Wide_20160308_110428_R60_4x_15-5609VD\_stiched_image\delta_X_-14_delta_Y_10\';

file_path=[folder_path_without_index 'stiched_image_offset8.8 micron_X6_Y8'];
fin=fopen(file_path);
A=fread(fin,[frame_width,frame_height],'float32','b');
%%
Start_X=1500;
Start_Y=1000;


A_partial=A(Start_X+(1:648),Start_Y+(1:488));

imagesc(A_partial);
%%
im=repmat(A_partial,[1 1 3]);
k=25;
m=5;

[l, Am, Sp, d]=slic(im, k, m);

%%

imagesc(l);