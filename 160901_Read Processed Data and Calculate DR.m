clear all
fclose('all')
%%
Output_File_Path='C:\Users\Owner\Desktop\160901_Text files\';
last_folder_name='160830_PCO_4160fps_4x_200microsec_200x200_forearm';
%Axial_Decimation_Factor=[1 2 4 8 16 32 64];
ave_factor=[1 2:2:64];
for p=1:length(ave_factor)
    Data_Save_Folder='F:\P1.2 Test\Processed Data\';

    Processed_Data_Path=[Data_Save_Folder last_folder_name sprintf('_Ave_Factor_%d.bin',ave_factor(p))];


    Row=200;
    Colomn=200;


    fin = fopen(Processed_Data_Path);
    Image_Temp=fread(fin,[Row,inf],'double');
    fclose(fin);
    %
    Colomn_Total=size(Image_Temp,2);

    Frame=Colomn_Total/Colomn;
    Image_Stack=zeros(Row,Colomn,Frame);
    for r=1:Frame
        Image_Stack(:,:,r)=Image_Temp(:,(1+(r-1)*Colomn):(r*Colomn));
    end
    %

    X_Show=size(Image_Stack,1)/2+10;
    Y_Show=size(Image_Stack,2)/2+10;

    Array_Show=squeeze(Image_Stack(X_Show,Y_Show,:));
    Array_Mean=squeeze(mean(mean(Image_Stack,1),2));
    plot(Array_Mean);

    
    Signal_Start_Index=30;
    Signal_End_Index=50;
    
    Signal_Max=max(Image_Stack(:,:,Signal_Start_Index:Signal_End_Index),[],3);
    
    Mean_Signal_Max(p)=mean(Signal_Max(:));
    %
    Noise_Start_Index=600;%100;
    Noise_End_Index=900;%400;

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

 plot(ave_factor,Mean_Noise_STD_norm);
 xlim([0 70]);

 xlabel('Averaging Factor');
 ylabel('Normalized Noise (ADU)');

 
 
 plot(ave_factor,Mean_Mean.*ave_factor);
 xlim([0 70]);
 xlabel('Averaging Factor');
 ylabel('Mean Value of Backgroud after 4-point Calculation (ADU)');
 
 
 plot(ave_factor,Mean_Noise_STD.*ave_factor);
 xlim([0 70]);
 ylim([0 80]);
  xlabel('Averaging Factor');

 ylabel('Noise of Backgroud after 4-point Calculation (ADU)');

 
    plot(Array_Mean);
    xlim([0 1000]);
    ylim([0 200]);

 
    plot(ave_factor,Mean_Signal_Max);
 xlim([0 70]);
  xlabel('Averaging Factor');

 
 Output=[ave_factor' Mean_Mean' Mean_Noise_STD' Mean_Noise_STD_norm' Mean_Signal_Max'];
 
dlmwrite([Output_File_Path last_folder_name '.txt'],Output,'delimiter','\t','newline','pc','precision', '%.6f');
