clear all
%%
root_folder_path='I:\161228_e2v EM1 new setup\tissue';
%root_folder_path='I:\161223_e2v EM1\CeYAG\Different Exp Time';
Current_Array=[22 24 25 26 27 28];%[3:3:24];
it=210;
Horizontal_Size=1024;
Size=25;

% Image Generation
SR_Lateral=0.888;
SR_Axial=0.279;
Axial_ave_Factor=1;

%

H_Offset_View=200;
H_Range_View=600;
Horizontal_ROI_View=[H_Offset_View H_Offset_View+H_Range_View];

%
H_Offset=200;
H_Range=600;
Horizontal_ROI=[H_Offset H_Offset+H_Range];

If_norm=1; %0: no normalization, 1: with normalization

% Feature for Filtering mode
Filtering_Ratio=0.02;       %意思是, 若Array長度為A, 則將外側(低頻)之A*Filtering_Ratio個pixels設為0
% Feature for 4-point mode
Averaging_Factor=16;
Product_of_Axial_Decimation_Factor_and_Ave_Factor=16;
N=4;          %for N-point
Unbias_Estimator=0.9213;

Mean_Array_Map=zeros(Horizontal_ROI(2)-Horizontal_ROI(1)+1,length(Current_Array));
STD_Array_Map=zeros(Horizontal_ROI(2)-Horizontal_ROI(1)+1,length(Current_Array));

Mean_Array_for_norm=zeros(1,length(Current_Array));

for p=1:length(Current_Array)
    
     %% file reading
    file_name=sprintf('%g_it%g',Current_Array(p),it);
    folder_path=[root_folder_path '\'];
    file_path=[folder_path file_name];

    cd(folder_path);

    Data_Image_Raw=double(imread(file_path,'tiff'));
    Data_Image_Raw_Decimated=Data_Image_Raw(1:floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/Averaging_Factor):end);
       
    if If_norm == 0
        Mean_Array_Map(Horizontal_ROI(1):Horizontal_ROI(2),p)=mean(Data_Image_Raw(:,Horizontal_ROI(1):Horizontal_ROI(2)),1);
        STD_Array_Map(Horizontal_ROI(1):Horizontal_ROI(2),p)=std(Data_Image_Raw(:,Horizontal_ROI(1):Horizontal_ROI(2)),1,1);
        Data_Image=Data_Image_Raw;
    elseif If_norm == 1
        Mean_Array_for_norm=mean(Data_Image_Raw(:,Horizontal_ROI(1):Horizontal_ROI(2)),2); % 改成可調ROI
        All_Mean=mean(Mean_Array_for_norm);
        Mean_Image_for_norm=repmat(Mean_Array_for_norm,[1 size(Data_Image_Raw,2)]);
        Data_Image=Data_Image_Raw./Mean_Image_for_norm*All_Mean;
    end
    %% Ave before N-point and N-point
    Averaged_Length=floor(size(Data_Image,1)/Averaging_Factor);
    Temp=0;
    for q=1:Averaging_Factor
       Temp=Temp+Data_Image((Averaging_Factor-(q-1)):Averaging_Factor:(Averaging_Factor*Averaged_Length)-(q-1),:);
    end
    Data_Image_Ave=Temp/Averaging_Factor;
    Data_Image_Npoint=zeros(floor(size(Data_Image_Ave,1)/N),size(Data_Image_Ave,2));
    for q=1:size(Data_Image_Npoint,1)
        Temp=Data_Image_Ave(((q-1)*N+1):(q*N),:);
        Data_Image_Npoint(q,:)=((N*sum(Temp.^2,1)-sum(Temp,1).^2).^0.5)*(2^0.5)/N;
    end
    
    %% Axial Ave after N-point
    Temp=0;
    Axial_Length_Original=size(Data_Image_Npoint,1);
    Axial_Length_Used=floor(Axial_Length_Original/Axial_ave_Factor)*Axial_ave_Factor;
    Reduced_Length=Axial_Length_Used/Axial_ave_Factor;
    for q=1:Axial_ave_Factor
    Temp=Temp+Data_Image_Npoint((Axial_ave_Factor-(q-1)):Axial_ave_Factor:(Axial_ave_Factor*Reduced_Length)-(q-1),:);
    end
    Data_Image_Reduced=Temp/Axial_ave_Factor;

    %% Image generation
    C_max=32;
    C_min=0.8;

    Data_Image_Npoint_normalized=(Data_Image_Reduced(:,Horizontal_ROI_View(1):Horizontal_ROI_View(2))-C_min)/(C_max-C_min);
    Data_Image_Npoint_normalized(Data_Image_Npoint_normalized<0)=0;
    Data_Image_Npoint_normalized(Data_Image_Npoint_normalized>1)=1;



    %%  Data Saving (Image only)

    Data_Save_Folder='I:\Processed Data\';
    Data_Name='161229_e2v crosssection';
    Processed_Data_Path=[Data_Save_Folder Data_Name '_' file_name];

    if If_norm ==0
        imwrite(Data_Image_Npoint_normalized,[Processed_Data_Path sprintf('_Axial Ave %g',Axial_ave_Factor) sprintf('_Ave %g',Averaging_Factor) '_Bscan.png'],'png');
    elseif If_norm ==1
        imwrite(Data_Image_Npoint_normalized,[Processed_Data_Path sprintf('_Axial Ave %g',Axial_ave_Factor) sprintf('_Ave %g',Averaging_Factor) '_Bscan_norm.png'],'png');
    end
    
    disp(p);
end

subplot(1,1,1)
imagesc(Data_Image_Npoint_normalized);
%axis off
%axis equal
caxis([0 1]);
X_Label_Array=[0  100 200 300 400 500];
X_Label_Index_Array=X_Label_Array/SR_Lateral;

Y_Label_Array=[0 50 100 150 200 250 300 350];
Y_Label_Index_Array=round(Y_Label_Array/SR_Axial/Axial_ave_Factor);
set(gca,'xtick',X_Label_Index_Array,'xticklabel',X_Label_Array,'ytick',Y_Label_Index_Array,'yticklabel',Y_Label_Array);
xlabel('Lateral Position (micron)');
ylabel('Axial Position (micron)');
colormap(gray);
