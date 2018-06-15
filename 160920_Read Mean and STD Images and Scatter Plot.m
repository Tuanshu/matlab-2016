clear all


Model='Mikrotron-New-2';

Power_Array=[1:2:11];

if strcmp(Model,'Mikrotron')
    Second_last_folder_name='160920_Diff Power';
    Bin_Adj_Factor=16;
    Row=1024;
    Colomn=8;
    Power_Array=[0 1:1:13];

elseif strcmp(Model,'Mikrotron-Old')
    Second_last_folder_name='160919_Diff Power';
    Bin_Adj_Factor=16;
    Row=1024;
    Colomn=8;
    Power_Array=[1:2:11];

elseif strcmp(Model,'Mikrotron-New')
    Second_last_folder_name='160921_Diff Power';
    Bin_Adj_Factor=16;
    Row=1024;
    Colomn=8;
    Power_Array=[0 1:1:13];

elseif strcmp(Model,'Mikrotron-New-2')
    Second_last_folder_name='160921_Diff Power_2';
    Bin_Adj_Factor=16;
    Row=1024;
    Colomn=8;
    Power_Array=[0 1:1:14];

    
elseif strcmp(Model,'PCO')
    Second_last_folder_name='160826_HS4_200x200_defocus_with C mount_4160fps_10microsec';
    Bin_Adj_Factor=1;
    Row=200;
    Colomn=200;
    Power_Array=[1:2:11];


elseif strcmp(Model,'Imperx')
    Second_last_folder_name='On System_Different Power';
    Bin_Adj_Factor=1;
    Row=648;
    Colomn=3;
    Power_Array=[1:2:13];

end

Data_Save_Folder=['F:\Processed Data\' Second_last_folder_name '\'];

If_ROI=1;
If_Central_ROI=1;

Width_ROI=[504];
Height_ROI=[5];
Central_ROI_Half_Width=20;
Central_ROI_Half_Height=1;
X_Offset=1;
Y_Offset=1;
Size=25;



Color_Matrix=   [1 1 0;
                 1 0 1;
                 0 1 1;
                 1 0 0;
                 0 1 0;
                 0 0 1;];

Mean_All=[];
STD_All=[];
for r=1:length(Power_Array)

    last_folder_name=sprintf('%dmW',Power_Array(r));

    fid = fopen([Data_Save_Folder 'Diff_Power_' last_folder_name '_Mean.bin']);
    Mean=fread(fid, [Row,Colomn], 'double');
    fclose(fid);

    fid = fopen([Data_Save_Folder 'Diff_Power_' last_folder_name '_STD.bin']);
    STD=fread(fid, [Row,Colomn], 'double');
    fclose(fid);

    if r==1
        X_Central=ceil(size(Mean,1)/2);
        Y_Central=ceil(size(Mean,2)/2);
        
        if If_Central_ROI ==1
            Width_ROI=[floor(X_Central+X_Offset-Central_ROI_Half_Width):(floor(X_Central+X_Offset+Central_ROI_Half_Width)-1)];
            Height_ROI=[floor(Y_Central+Y_Offset-Central_ROI_Half_Height):(floor(Y_Central+Y_Offset+Central_ROI_Half_Height)-1)];

        end
    end
    
    if If_ROI == 1
        Mean_ROI=Mean(Width_ROI,Height_ROI);
        STD_ROI=STD(Width_ROI,Height_ROI);
    else
        Mean_ROI=Mean;
        STD_ROI=STD;
    end

    
    Mean_All(:,:,r)=Mean_ROI;
    STD_All(:,:,r)=STD_ROI;
    scatter((Mean_ROI(:)*Bin_Adj_Factor),((STD_ROI(:)*Bin_Adj_Factor).^2),Size,Color_Matrix(rem(r,size(Color_Matrix,1))+1,:),'filled');
    
    
    
    hold on
end


xlim([0 4096])
ylim([0 400])
xlabel('Signal Mean (DN)');
ylabel('Signal Variance (DN^2)');

hold off
%% Noise Fitting
Fit_Range=[1 1100];
Mean_All_Array=Mean_All(:);

STD_All_Array=STD_All(:);

[Mean_All_Array_Sort index]=sort(Mean_All_Array);
Mean_All_Array_Sort_For_Fit=Mean_All_Array_Sort((Mean_All_Array_Sort-mean(Fit_Range/Bin_Adj_Factor))<abs(mean(Fit_Range/Bin_Adj_Factor)-Fit_Range(1)/Bin_Adj_Factor));
index=index((Mean_All_Array_Sort-mean(Fit_Range/Bin_Adj_Factor))<abs(mean(Fit_Range/Bin_Adj_Factor)-Fit_Range(1)/Bin_Adj_Factor));

STD_All_Array_Sort_For_Fit=STD_All_Array(index);
g = fittype( @(a,b,c,x) (a*x.^2+b*x+c)); %a: excess, b: shot, c: dark
FitResult = fit(Mean_All_Array_Sort_For_Fit*Bin_Adj_Factor, (STD_All_Array_Sort_For_Fit*Bin_Adj_Factor).^2, g, 'StartPoint', [0.015, 0.0569, -0]);
FitResult.a
FitResult.b
FitResult.c

STD_All_Array_fit=(FitResult.a*(Mean_All_Array_Sort*Bin_Adj_Factor).^2+FitResult.b*Mean_All_Array_Sort*Bin_Adj_Factor+FitResult.c).^0.5;

hold on
plot(Mean_All_Array_Sort*Bin_Adj_Factor,(STD_All_Array_fit).^2,'-');


title(sprintf('%s, 2nd Coef=%g',Model,FitResult.b)) 

hold off

%% Pixel-by-pixel Fitting
a_Map=zeros(size(Mean_All,1),size(Mean_All,2));
b_Map=zeros(size(Mean_All,1),size(Mean_All,2));
c_Map=zeros(size(Mean_All,1),size(Mean_All,2));
MAX_STD_VAR_Map=zeros(size(Mean_All,1),size(Mean_All,2));
p_record=[];
q_record=[];
V_Th=425;
for p=1:size(Mean_ROI,1)
    for q=1:size(Mean_ROI,2)
        Mean_Array=squeeze(Mean_All(p,q,:));
        STD_Array=squeeze(STD_All(p,q,:));
        if max(STD_Array*Bin_Adj_Factor)^2>V_Th
            p_record=[p_record p];
            q_record=[q_record q];
        end
            
        FitResult = fit(Mean_Array, STD_Array.^2, g, 'StartPoint', [0.015, 0.0569, -0]);
        VAR_Array_fit=FitResult.a*Mean_Array.^2+FitResult.b*Mean_Array+FitResult.c;
        a_Map(p,q)=FitResult.a;
        b_Map(p,q)=FitResult.b;
        c_Map(p,q)=FitResult.c;
        MAX_STD_VAR_Map(p,q)=max(STD_Array.^2-VAR_Array_fit);
        h=plot(Mean_Array*Bin_Adj_Factor,(STD_Array*Bin_Adj_Factor).^2);
        set(h,'Color',Color_Matrix(rem((p-1)*size(Mean_ROI,2)+q,size(Color_Matrix,1))+1,:))
        hold on
    end
    disp(p)
end
hold off

xlim([0 4096])
ylim([0 400])
xlabel('Signal Mean (DN)');
ylabel('Signal Variance (DN^2)');
title(sprintf('%s, Pixel by Pixel',Model)) 

% 
% subplot(3,1,1)
% imagesc(a_Map);
% colorbar 
% subplot(3,1,2)
% imagesc(b_Map);
% colorbar 
% 
% subplot(3,1,3)
% imagesc(c_Map);
% colorbar 
%%
imwrite(MAX_STD_VAR_Map,'C:\Users\Owner\Desktop\160920_Images\MAX_STD_VAR_Map.png')