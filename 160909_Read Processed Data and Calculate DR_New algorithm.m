clear all
fclose('all')
%%
Output_File_Path='I:\PCO\Processed Data\';
last_folder_name='in vivo_21_after adj 7';
%Axial_Decimation_Factor=[1 2 4 8 16 32 64];
ave_factor=[16];
for p=1:length(ave_factor)
    Data_Save_Folder='I:\PCO\Processed Data\';

    Processed_Data_Path=[Data_Save_Folder last_folder_name sprintf('_Ave_Factor_%d.raw',ave_factor(p))];


    Row=1000;
    Colomn=8;


    fin = fopen(Processed_Data_Path);
    Image_Temp=fread(fin,[Row,Inf],'double');
    fclose(fin);
    %
    Colomn_Total=size(Image_Temp,2);

    Frame=Colomn_Total/Colomn;
    Image_Stack=zeros(Row,Colomn,Frame);
    for r=1:Frame
        Image_Stack(:,:,r)=Image_Temp(:,(1+(r-1)*Colomn):(r*Colomn));
    end
    %

%     X_Show=size(Image_Stack,1)/2+10;
%     Y_Show=size(Image_Stack,2)/2+10;
% 
%     Array_Show=squeeze(Image_Stack(X_Show,Y_Show,:));
%     plot(Array_Mean);

    Array_Mean=squeeze(mean(mean(Image_Stack,1),2));

    Signal_Start_Index=30;
    Signal_End_Index=50;
    
    Signal_Max=max(Image_Stack(:,:,Signal_Start_Index:Signal_End_Index),[],3);
    
    Mean_Signal_Max(p)=mean(Signal_Max(:));
    %
    Noise_Start_Index=600;%100;
    Noise_End_Index=800;%400;

    Glass_Interface_Position=214;

    Noise_STD(:,:)=std(Image_Stack(:,:,Noise_Start_Index:Noise_End_Index),0,3);

    Mean(:,:)=mean(Image_Stack(:,:,Noise_Start_Index:Noise_End_Index),3);
    
    Noise_STD_norm=Noise_STD./Mean;
    
    
    imagesc(Noise_STD);
    imagesc(Noise_STD_norm);

    Mean_Noise_STD(p)=mean(Noise_STD(:))
    %Mean_Noise_STD(p)=Noise_STD(X_Show,Y_Show)
    Mean_Noise_STD_norm(p)=mean(Noise_STD_norm(:))
    Mean_Mean(p)=mean(Mean(:))

end

%%
QQQ=1;
C_max=20;
C_min=6;

Crosssection_Show(:,:)=Image_Stack(:,QQQ,:);
Crosssection_Show_Norm=(Crosssection_Show'-C_min)/(C_max-C_min);
imagesc(Crosssection_Show_Norm);
colormap(gray);
caxis([0 1]);
axis equal
xlim([1 size(Crosssection_Show_Norm,1)]);
ylim([1 size(Crosssection_Show_Norm,2)]);
axis off

%% Identify the glass interface position
Min_Distance_between_Glass_and_Tissue=8;
[Max_Map_value Max_Map_index]=max(Image_Stack,[],3);

se = strel('disk',10);
Max_Map_index_Opened=imopen(Max_Map_index,se);

imagesc(Max_Map_index_Opened);

Starting_index_Map=Max_Map_index_Opened+Min_Distance_between_Glass_and_Tissue;

%% To generate data with isotropic resolution
Axial_ave_Factor=2;
Temp=0;
Axial_Length_Original=size(Image_Stack,3);
Axial_Length_Used=floor(Axial_Length_Original/Axial_ave_Factor)*Axial_ave_Factor;
Reduced_Length=Axial_Length_Used/Axial_ave_Factor;
for p=1:Axial_ave_Factor
   Temp=Temp+Image_Stack(:,:,(Axial_ave_Factor-(p-1)):Axial_ave_Factor:(Axial_ave_Factor*Reduced_Length)-(p-1));
end
Reduced_Stack=Temp/Axial_ave_Factor;
Reduced_Image=squeeze(mean(Reduced_Stack,2))';
Reduced_Image=Reduced_Image(1:size(Reduced_Image,1),:);


%% Blurring the 3D data
Blur_Size=4;
H=fspecial('gaussian',Blur_Size*4,Blur_Size);
%H=H/max(H(:));
clear Image_Stack_Blur
for p=1:size(Reduced_Stack,2)
    Image_Stack_Blur(:,p,:) = imfilter(Reduced_Stack(:,p,:) ,H);
end

%% Opening the 3D data
Open_Size=5;
se = strel('disk',Open_Size);
Image_Stack_Opened=imopen(Image_Stack_Blur,se);

%% Closeing the 3D data
Close_Size=10;
se = strel('disk',Close_Size);
Image_Stack_Closed=imclose(Image_Stack_Blur,se);

%% Find the BND value
BND_Frame_Index_Range=[1 15];

BND_Stack=Image_Stack_Closed(:,:,BND_Frame_Index_Range(1):BND_Frame_Index_Range(2));
BND=mean(BND_Stack(:));


Image_Stack_Closed_BND=Image_Stack_Closed-BND;
%%

QQQ=1;
C_max=3;
C_min=0;

Crosssection_Show(:,:)=Image_Stack_Blur(:,QQQ,:);
Crosssection_Show_Norm=(Crosssection_Show'-C_min)/(C_max-C_min);
Crosssection_Show_Norm(Crosssection_Show_Norm>1)=1;
Crosssection_Show_Norm(Crosssection_Show_Norm<0)=0;

Crosssection_Show_Open(:,:)=Image_Stack_Opened(:,QQQ,:);
Crosssection_Show_Norm_Open=(Crosssection_Show_Open'-C_min)/(C_max-C_min);
Crosssection_Show_Norm_Open(Crosssection_Show_Norm_Open>1)=1;
Crosssection_Show_Norm_Open(Crosssection_Show_Norm_Open<0)=0;

Crosssection_Show_Closed(:,:)=Image_Stack_Closed(:,QQQ,:);
Crosssection_Show_Norm_Closed=(Crosssection_Show_Closed'-C_min)/(C_max-C_min);
Crosssection_Show_Norm_Closed(Crosssection_Show_Norm_Closed>1)=1;
Crosssection_Show_Norm_Closed(Crosssection_Show_Norm_Closed<0)=0;


Crosssection_Show_Closed_BND(:,:)=Image_Stack_Closed_BND(:,QQQ,:);
Crosssection_Show_Norm_Closed_BND=(Crosssection_Show_Closed_BND'-C_min)/(C_max-C_min);
Crosssection_Show_Norm_Closed_BND(Crosssection_Show_Norm_Closed_BND>1)=1;
Crosssection_Show_Norm_Closed_BND(Crosssection_Show_Norm_Closed_BND<0)=0;




subplot(3,1,1)
imagesc(Crosssection_Show_Norm);
colormap(gray);
caxis([0 1]);
axis equal
xlim([0 size(Crosssection_Show_Norm,2)]);
ylim([0 size(Crosssection_Show_Norm,1)]);
axis off

subplot(3,1,2)
imagesc(Crosssection_Show_Norm_Closed);
colormap(gray);
caxis([0 1]);
axis equal
xlim([0 size(Crosssection_Show_Norm_Closed,2)]);
ylim([0 size(Crosssection_Show_Norm_Closed,1)]);
axis off

subplot(3,1,3)
imagesc(Crosssection_Show_Norm_Closed_BND);
colormap(gray);
caxis([0 1]);
axis equal
xlim([0 size(Crosssection_Show_Norm_Closed_BND,2)]);
ylim([0 size(Crosssection_Show_Norm_Closed_BND,1)]);
axis off

%Image_Stack_Glass_Interface_Removed=Image_Stack(:,:,)

%  plot(ave_factor,Mean_Noise_STD_norm);
%  xlim([0 70]);
% 
%  xlabel('Averaging Factor');
%  ylabel('Normalized Noise (ADU)');
% 
%  
%  
%  plot(ave_factor,Mean_Mean.*ave_factor);
%  xlim([0 70]);
%  xlabel('Averaging Factor');
%  ylabel('Mean Value of Backgroud after 4-point Calculation (ADU)');
%  
%  
%  plot(ave_factor,Mean_Noise_STD.*ave_factor);
%  xlim([0 70]);
%  ylim([0 80]);
%   xlabel('Averaging Factor');
% 
%  ylabel('Noise of Backgroud after 4-point Calculation (ADU)');
% 
%  
%     plot(Array_Mean);
%     xlim([0 1000]);
%     ylim([0 200]);
% 
%  
% plot(ave_factor,Mean_Noise_STD./Mean_Mean);
% ylim([0.3 0.5]);
% %xlim([0 70]);
% xlabel('Averaging Factor');
% ylabel('Normalized Noise (Noise/BND)');

 
 Output=[ave_factor' Mean_Mean' Mean_Noise_STD' Mean_Noise_STD_norm' Mean_Signal_Max'];
 
dlmwrite([Output_File_Path last_folder_name '.txt'],Output,'delimiter','\t','newline','pc','precision', '%.6f');
