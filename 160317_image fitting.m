clear all

ImA=sum(imread('D:\160317_OM\15-33455a\Blank\Image46.TIF'),3);
%ImA=max(ImA(:))-ImA;
ImA=ImA/max(ImA(:));
%ImB=double(imread('G:\Everday Experiements\WideScan_20160316_163444_15-33455a_C3_S70_R65_blaken_1x\_stiched_image\stiched_image_offset40.0 micron_X12_Y16.png'));
ImB=sum(imread('D:\160317_OM\15-33455a\Blank\Image46.TIF'),3);
ImB=ImB/max(ImB(:));
%% Test- also locate A
Pre_angle=0;
X_start=150;   %1: 1050 2: 1100    3: 2235     4: 1050
Y_start=50;   %1: 875 2: 1375     3: 2950     4: 2940
X_size=700*1*0.6;
Y_size=X_size/3*4;

ImA=ImA(X_start:(X_start+X_size-1),Y_start:(Y_start+Y_size-1));
ImA=imrotate(ImA,Pre_angle,'crop');


%%       % Patial
Pre_angle=0;
X_start=100;   %1: 1050 2: 1100    3: 2235     4: 1050
Y_start=100;   %1: 875 2: 1375     3: 2950     4: 2940
X_size=700*1*0.6;
Y_size=X_size/3*4;

ImB_Patial=ImB(X_start:(X_start+X_size-1),Y_start:(Y_start+Y_size-1));
ImB_Prerotate=imrotate(ImB_Patial,Pre_angle,'crop');


%imwrite(stiched_image_rotated_patial,[sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_offset%.1f micron_X%d_Y%d_rotated',(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index),'.png']);

%
subplot(2,1,1)
imagesc(ImA);
colormap(gray);
axis equal
xlim([0 size(ImA,2)]);
ylim([0 size(ImA,1)]);

subplot(2,1,2)
imagesc(ImB_Prerotate);
colormap(gray);
axis equal
xlim([0 size(ImB_Prerotate,2)]);
ylim([0 size(ImB_Prerotate,1)]);

%% Blur
% Blur_Window_Size=100;
% filter=fspecial('gaussian',Blur_Window_Size,round(Blur_Window_Size/2));
% filter=filter/sum(sum(filter));
% 
% ImA_Blur=conv2(ImA,filter,'same');
% 
% ImB_Blur=conv2(ImB_Patial,filter,'same');
% %%
% subplot(2,1,1)
% imagesc(ImA_Blur);
% colormap(gray);
% axis equal
% xlim([0 size(ImA_Blur,2)]);
% ylim([0 size(ImA_Blur,1)]);
% 
% subplot(2,1,2)
% imagesc(ImB_Blur);
% colormap(gray);
% axis equal
% xlim([0 size(ImB_Blur,2)]);
% ylim([0 size(ImB_Blur,1)]);

%% Resolution reduction
Factor_B=10;
Factor_A=round(Factor_B*size(ImA,2)/size(ImB_Prerotate,2));

ImA_Reduced=downsample(downsample(ImA,Factor_A)',Factor_A)';
ImB_Reduced=downsample(downsample(ImB_Prerotate,Factor_B)',Factor_B)';
ImA_Reduced=ImA_Reduced/max(ImA_Reduced(:));
ImB_Reduced=ImB_Reduced/max(ImB_Reduced(:));

subplot(2,1,1)
imagesc(ImA_Reduced);
colormap(gray);
axis equal
xlim([0 size(ImA_Reduced,2)]);
ylim([0 size(ImA_Reduced,1)]);

subplot(2,1,2)
imagesc(ImB_Reduced);
colormap(gray);
axis equal
xlim([0 size(ImB_Reduced,2)]);
ylim([0 size(ImB_Reduced,1)]);

%% BW
Th_A=0.73;
Th_B=0.73;
TH_Size=100;

ImA_BW=bwareaopen(im2bw(ImA_Reduced,Th_A),TH_Size);
ImB_BW=bwareaopen(im2bw(ImB_Reduced,Th_B),TH_Size);

subplot(2,1,1)
imagesc(ImA_BW);
colormap(gray);
axis equal
xlim([0 size(ImA_BW,2)]);
ylim([0 size(ImA_BW,1)]);

subplot(2,1,2)
imagesc(ImB_BW);
colormap(gray);
axis equal
xlim([0 size(ImB_BW,2)]);
ylim([0 size(ImB_BW,1)]);
%%
% angle=7;
% ImB_Rotated=imrotate(ImB_Reduced,angle);
% subplot(2,1,1)
% imagesc(ImA_Reduced);
% colormap(gray);
% axis equal
% xlim([0 size(ImA_Reduced,2)]);
% ylim([0 size(ImA_Reduced,1)]);
% 
% subplot(2,1,2)
% imagesc(ImB_Rotated);
% colormap(gray);
% axis equal
% xlim([0 size(ImB_Rotated,2)]);
% ylim([0 size(ImB_Rotated,1)]);


%%
% ImA_edge=edge(ImA_Reduced,'canny');
% ImB_edge=edge(ImB_Reduced,'canny');
% %%
% subplot(2,1,1)
% imagesc(ImA_edge);
% colormap(gray);
% axis equal
% xlim([0 size(ImA_edge,2)]);
% ylim([0 size(ImA_edge,1)]);
% 
% subplot(2,1,2)
% imagesc(ImB_edge);
% colormap(gray);
% axis equal
% xlim([0 size(ImB_edge,2)]);
% ylim([0 size(ImB_edge,1)]);
%%

ImA_Dist=bwdist(ImA_BW);
ImB_Dist=bwdist(ImB_BW);
%%
subplot(2,1,1)
imagesc(ImA_Dist);
colormap(gray);
axis equal
xlim([0 size(ImA_Dist,2)]);
ylim([0 size(ImA_Dist,1)]);

subplot(2,1,2)
imagesc(ImB_Dist);
colormap(gray);
axis equal
xlim([0 size(ImB_Dist,2)]);
ylim([0 size(ImB_Dist,1)]);

%%  先旋轉, offset就依賴數學的力量
Range_Angle=11;  
Delta_Angle=0.1;

Angle_Array=-1*Range_Angle:Delta_Angle:Range_Angle;
Window_Size_X=min(size(ImA_Dist,1),size(ImB_Dist,1));
Window_Size_Y=min(size(ImA_Dist,2),size(ImB_Dist,2));

%Window_Size_X=round(min(size(ImA_Dist,1),size(ImB_Dist,1))-size(ImB_Dist,2)*sin(Range_Angle/180*pi));
%Window_Size_Y=round(min(size(ImA_Dist,2),size(ImB_Dist,2))-size(ImB_Dist,1)*sin(Range_Angle/180*pi));

%%
% 
% subplot(2,1,1)
% imagesc(Temp_A);
% colormap(gray);
% axis equal
% xlim([0 size(Temp_A,2)]);
% ylim([0 size(Temp_A,1)]);
% 
% subplot(2,1,2)
% imagesc(Temp_B);
% colormap(gray);
% axis equal
% xlim([0 size(Temp_B,2)]);
% ylim([0 size(Temp_B,1)]);

%%
Temp_A=ImA_Reduced(round((size(ImA_Dist,1)-Window_Size_X)/2)+(1:Window_Size_X),round((size(ImA_Dist,2)-Window_Size_Y)/2)+(1:Window_Size_Y));
Temp_B=ImB_Reduced(round((size(ImB_Dist,1)-Window_Size_X)/2)+(1:Window_Size_X),round((size(ImB_Dist,2)-Window_Size_Y)/2)+(1:Window_Size_Y));

FFT_A=fft2(Temp_A);
FFT_B=fft2(Temp_B);            
V_Temp=ifftshift(ifftshift(ifft2(conj(FFT_A).*FFT_B),1),2);

imagesc(V_Temp);
%% Grid generation
[value index]=max(V_Temp(:));



X_grid=repmat((1:size(V_Temp,1))',1,size(V_Temp,2));
Y_grid=repmat((1:size(V_Temp,2)),size(V_Temp,1),1);
X_grid_array=X_grid(:);
Y_grid_array=Y_grid(:);

%V_X=sum(sum(V_Temp.*X_grid))/sum(sum(V_Temp))-size(V_Temp,1)/2;
%V_Y=sum(sum(V_Temp.*Y_grid))/sum(sum(V_Temp))-size(V_Temp,2)/2;

%%
V=0;
Temp_A=ImA_Reduced(round((size(ImA_Dist,1)-Window_Size_X)/2)+(1:Window_Size_X),round((size(ImA_Dist,2)-Window_Size_Y)/2)+(1:Window_Size_Y));


Zero_Padding_N=1;

Temp_A_ZP=zeros((2*Zero_Padding_N+1)*size(Temp_A,1),(2*Zero_Padding_N+1)*size(Temp_A,2));
Temp_B_ZP=zeros((2*Zero_Padding_N+1)*size(Temp_B,1),(2*Zero_Padding_N+1)*size(Temp_B,2));
Temp_A_ZP((Zero_Padding_N*size(Temp_A,1)+1):((Zero_Padding_N+1)*size(Temp_A,1)),(Zero_Padding_N*size(Temp_A,2)+1):((Zero_Padding_N+1)*size(Temp_A,2)))=Temp_A;

Value_Temp=0;
for p=1:length(Angle_Array)
    ImB_Rotate=imrotate(ImB_Reduced,Angle_Array(p));
    Temp_B=ImB_Rotate(round((size(ImB_Rotate,1)-Window_Size_X)/2)+(1:Window_Size_X),round((size(ImB_Rotate,2)-Window_Size_Y)/2)+(1:Window_Size_Y));
    Temp_B_ZP((Zero_Padding_N*size(Temp_B,1)+1):((Zero_Padding_N+1)*size(Temp_B,1)),(Zero_Padding_N*size(Temp_B,2)+1):((Zero_Padding_N+1)*size(Temp_B,2)))=Temp_B;
    FFT_A=fft2(Temp_A_ZP);
    FFT_B=fft2(Temp_B_ZP);            
    V_Temp=ifftshift(ifftshift(ifft2(conj(FFT_A).*FFT_B),1),2);
    V=V_Temp((Zero_Padding_N*size(Temp_B,1)+1):((Zero_Padding_N+1)*size(Temp_B,1)),(Zero_Padding_N*size(Temp_B,2)+1):((Zero_Padding_N+1)*size(Temp_B,2)));
    [value index]=max(V(:));
    if value>Value_Temp
        Angle_record=Angle_Array(p);
        V_X=X_grid_array(index)-size(V,1)/2;
        V_Y=Y_grid_array(index)-size(V,2)/2;
        Value_Temp=value;
    end

    
end
V_X
V_Y
Angle_record


%%
% X_Delete=25;
% Y_Delete=25;
% 
% 
% angle=7;
% 
% ImB_Rotated=imrotate(ImB_Dist,angle);
% 
% ImB_Center=ImB_Rotated((X_Delete+1):(size(ImB_Rotated,1)-X_Delete),(Y_Delete+1):(size(ImB_Rotated,2)-Y_Delete));
% 
% X_Delete_A=round((size(ImA_Dist,1)-size(ImB_Rotated_Center,2))/2);
% Y_Delete_A=round((size(ImA_Dist,2)-size(ImB_Rotated_Center,2))/2);
% 
% 
% ImA_Center=ImA_Dist(X_Delete_A+1:(X_Delete_A+size(ImB_Rotated_Center,1)),Y_Delete_A+1:(Y_Delete_A+size(ImB_Rotated_Center,2)));
% %%
% subplot(2,1,1)
% imagesc(ImA_Center);
% colormap(gray);
% axis equal
% xlim([0 size(ImA_Center,2)]);
% ylim([0 size(ImA_Center,1)]);
% 
% subplot(2,1,2)
% imagesc(ImB_Center);
% colormap(gray);
% axis equal
% xlim([0 size(ImB_Center,2)]);
% ylim([0 size(ImB_Center,1)]);
% 
% 
% %%
% 
% 
% 
%         Fitting=ifftshift(ifft(conj(FFT_1_X).*FFT_2_X,[],1),1);
% 
% 
% angle=7;
% stiched_image_rotated=imrotate(stiched_image_write,angle);
