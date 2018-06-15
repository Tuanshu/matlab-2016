clear all
%%
Second_last_folder_name='160922_Perfect Dark';

root_folder_path=['F:\Imperx\' Second_last_folder_name '\'];
Data_Save_Folder=['F:\Processed Data\' Second_last_folder_name '\'];

mkdir(Data_Save_Folder);

Power_Array=[1:2:13];
If_smooth = 0;


Std_Mean_Array_1=ones([length(Power_Array) 1]);
Mean_Array_1=ones([length(Power_Array) 1]);
Std_Mean_Array_2=ones([length(Power_Array) 1]);
Mean_Array_2=ones([length(Power_Array) 1]);


for r=1:length(Power_Array)

    last_folder_name=sprintf('%dmW',Power_Array(r));
    folder_path=[root_folder_path last_folder_name '\'];

    Processed_Data_Path=[Data_Save_Folder last_folder_name '.bin'];

    cd(folder_path);

    % Data format related
    Row=648;
    Colomn=3;
    Byte_Skip=0;
    % Processing related
    ave_factor=1;
    N=4;
    micron_per_frame=0.2/ave_factor/N;
    Offset_1=0;
    Offset_2=0; 
    Max_Frame_Number=2000;
    %%
    file_list=dir(folder_path);
    file_list=file_list(4:end); %多skip一個rec檔
    Frame=min(length(file_list),Max_Frame_Number);
    %%
    Ave_Frame_Length=floor(Frame/ave_factor);
    Ave_Image_Stack=zeros(Row,Colomn,Ave_Frame_Length);
    Ave_Temp=zeros(Row,Colomn,ave_factor);
    X=[1:Frame];
    for p=1:Ave_Frame_Length
        for q=1:ave_factor

            file_path=[folder_path file_list((p-1)*ave_factor+q+Offset_1).name];

            fin=fopen(file_path);

            fseek(fin, Byte_Skip, 'bof');

            Ave_Temp(:,:,q)=fread(fin,[Row,Colomn],'uint16'); %*Frame   不知為何, 看起來就像是要除16
            
            fclose(fin);
        end
        Ave_Image_Stack(:,:,p)=mean(Ave_Temp,3);
        disp(p);
        X_ave(p)=mean(X(((p-1)*ave_factor+1+Offset_1):((p)*ave_factor+Offset_1)));
    end

    %% Show Image

    NN=1;

    Image_Show=Ave_Temp(:,:,1);

    imagesc(Image_Show);


    %% Z-plot along certian pixel
    X_Show=Row/2;
    Y_Show=Colomn;

    Array_Show=squeeze(Ave_Image_Stack(X_Show,Y_Show,:));

    plot(Array_Show);

    %%
    Window_Size=101;
    %yy = smooth(y,span)
    %W = smooth3(V,'filter',size,sd) 
    Array_Show_Smooth=smooth(Array_Show,Window_Size);
    X=1:length(Array_Show_Smooth);
    plot(X,Array_Show,X,Array_Show_Smooth);

    Array_Show_Div=Array_Show-Array_Show_Smooth;
    plot(X,Array_Show_Div);
    Array_Show_Div_Avaliable=Array_Show_Div((Window_Size+1):(length(Array_Show_Div)-Window_Size));
    Std=std(Array_Show_Div_Avaliable);

%     %%
%     Row_Center_1=500;
%     Colomn_Center_1=4;
%     Half_Size_1=4;
% 
% 
%     Row_Center_2=100;
%     Colomn_Center_2=100;
%     Half_Size_2=10;
% 
% 
%     Row_ROI_1=Row_Center_1+[-1*(Half_Size_1-1):Half_Size_1];
%     Colomn_ROI_1=Colomn_Center_1+[-1*(Half_Size_1-1):Half_Size_1];
% 
%     Row_ROI_2=Row_Center_2+[-1*(Half_Size_2-1):Half_Size_2];
%     Colomn_ROI_2=Colomn_Center_2+[-1*(Half_Size_2-1):Half_Size_2];

%     Ave_Image_Stack_ROI_1=Ave_Image_Stack(Row_ROI_1,Colomn_ROI_1,:);
    
    Ave_Image_Stack_ROI_1=Ave_Image_Stack;
%     Ave_Image_Stack_ROI_2=Ave_Image_Stack(Row_ROI_2,Colomn_ROI_2,:);
    %Ave_Image_Stack_ROI_1=Ave_Image_Stack;
    %Ave_Image_Stack_ROI_2=Ave_Image_Stack(Row_ROI_2,Colomn_ROI_2,:);

    % Show Image

    NN=285;
    Image_Show_1=Ave_Image_Stack_ROI_1(:,:,NN);
    %Image_Show_2=Ave_Image_Stack_ROI_2(:,:,NN);


    %
    if If_smooth == 1
        Ave_Image_Stack_Smooth_1=smooth3(Ave_Image_Stack_ROI_1,'box',[1 1 Window_Size]);
    else
        Ave_Image_Stack_Smooth_1=0;
    end
    
    Ave_Image_Stack_Div_1=Ave_Image_Stack_ROI_1-Ave_Image_Stack_Smooth_1;
    Ave_Image_Stack_Div_Avaliable_1=Ave_Image_Stack_Div_1(:,:,(Window_Size+1):(length(Array_Show_Div)-Window_Size));
    Ave_Image_Stack_Avaliable_1=Ave_Image_Stack_ROI_1(:,:,(Window_Size+1):(length(Array_Show_Div)-Window_Size));
    Ave_Image_Stack_Mean_1=mean(Ave_Image_Stack_Avaliable_1,3);
    Ave_Image_Stack_Std_1=std(Ave_Image_Stack_Div_Avaliable_1,0,3);

    fid = fopen([Data_Save_Folder 'Diff_Power_' last_folder_name '_Mean.bin'], 'w+');
    fwrite(fid, Ave_Image_Stack_Mean_1, 'double');
    fclose(fid);
    
    fid = fopen([Data_Save_Folder 'Diff_Power_' last_folder_name '_STD.bin'], 'w+');
    fwrite(fid, Ave_Image_Stack_Std_1, 'double');
    fclose(fid);
    

%     Ave_Image_Stack_Smooth_2=smooth3(Ave_Image_Stack_ROI_2,'box',[1 1 Window_Size]);
%     Ave_Image_Stack_Div_2=Ave_Image_Stack_ROI_2-Ave_Image_Stack_Smooth_2;
%     Ave_Image_Stack_Div_Avaliable_2=Ave_Image_Stack_Div_2(:,:,(Window_Size+1):(length(Array_Show_Div)-Window_Size));
%     Ave_Image_Stack_Std_2=std(Ave_Image_Stack_Div_Avaliable_2,0,3);





    subplot(2,1,1)
    imagesc(Image_Show_1);
    axis equal
    xlim([1 size(Image_Show_1,1)]);
    ylim([1 size(Image_Show_1,2)]);
    %caxis([min(min(Image_Show_1(:)),min(Image_Show_2(:))) max(max(Image_Show_1(:)),max(Image_Show_2(:)))*1.2]);
    axis off

%     subplot(2,2,2)
%     imagesc(Image_Show_2);
%     axis equal
%     xlim([1 size(Image_Show_2,1)]);
%     ylim([1 size(Image_Show_2,2)]);
%     caxis([min(min(Image_Show_1(:)),min(Image_Show_2(:))) max(max(Image_Show_1(:)),max(Image_Show_2(:)))*1.2]);
%     axis off


    subplot(2,1,2)
    imagesc(Ave_Image_Stack_Std_1);
    axis equal
    xlim([1 size(Image_Show_1,1)]);
    ylim([1 size(Image_Show_1,2)]);
    axis off

%     subplot(2,2,4)
%     imagesc(Ave_Image_Stack_Std_2);
%     axis equal
%     xlim([1 size(Image_Show_2,1)]);
%     ylim([1 size(Image_Show_2,2)]);
%     axis off


    %%

    Std_Mean_1=mean(Ave_Image_Stack_Std_1(:));
    Max_1=max(Ave_Image_Stack_Smooth_1(:));
    Mean_1=mean(Ave_Image_Stack_ROI_1(:));
% 
% 
%     Std_Mean_2=mean(Ave_Image_Stack_Std_2(:));
%     Max_2=max(Ave_Image_Stack_Smooth_2(:));
%     Mean_2=mean(Ave_Image_Stack_Smooth_2(:));
%     
    Std_Mean_Array_1(r)=Std_Mean_1;
    Mean_Array_1(r)=Mean_1;
%     
%     Std_Mean_Array_2(r)=Std_Mean_2;
%     Mean_Array_2(r)=Mean_2;
    
end    
    
disp(sprintf('STD1=%g',Std_Mean_1));
% disp(sprintf('MEAN1=%g',Mean_1));
% disp(sprintf('STD2=%g',Std_Mean_2));
% disp(sprintf('MEAN2=%g',Mean_2));
%%
FWC_Array=20000:100:40000;
Max_ADU=4096;

%QE_Array=0.50;
QE=0.5;
MSE_Array=zeros([length(FWC_Array) 1]);
for p=1:length(FWC_Array)
    FWC=FWC_Array(p);
    Std_Mean_Array_1_p=Std_Mean_Array_1/Max_ADU*FWC/QE;
    Mean_Array_1_p=Mean_Array_1/Max_ADU*FWC/QE;
     Std_Mean_Array_2_p=Std_Mean_Array_2/Max_ADU*FWC/QE;
%     Mean_Array_2_p=Mean_Array_2/Max_ADU*FWC/QE;

    Shot_1_p=Mean_Array_1_p.^0.5;
%     Shot_2_p=Mean_Array_2_p.^0.5;

    Std_Mean_Array_wo_Shot_1_p=(Std_Mean_Array_1_p.^2-Shot_1_p.^2).^0.5;
%     Std_Mean_Array_wo_Shot_2_p=(Std_Mean_Array_2_p.^2-Shot_2_p.^2).^0.5;

%     plot(Power_Array,Std_Mean_Array_1);
% 
%     plot(Power_Array,Mean_Array_1);
%     plot(Power_Array,Std_Mean_Array_1_p,Power_Array,Std_Mean_Array_wo_Shot_1_p);

    Matrix=[Power_Array; ones([1 length(Power_Array)])]';
    beta=Matrix\Std_Mean_Array_wo_Shot_1_p;
    Std_Mean_Array_wo_Shot_1_p_fit=beta(1)*Power_Array+beta(2);

%     plot(Power_Array,Std_Mean_Array_wo_Shot_1_p,Power_Array,Std_Mean_Array_wo_Shot_1_p_fit);

    MSE_Array(p)=sum(((Std_Mean_Array_wo_Shot_1_p'-Std_Mean_Array_wo_Shot_1_p_fit)./Std_Mean_Array_wo_Shot_1_p_fit).^2)^0.5;
end

[minvalue minindex]=min(MSE_Array);

FWC_Best=FWC_Array(minindex);
%%

FWC=FWC_Best;
Std_Mean_Array_1_p=Std_Mean_Array_1/Max_ADU*FWC/QE;
Mean_Array_1_p=Mean_Array_1/Max_ADU*FWC/QE;
% Std_Mean_Array_2_p=Std_Mean_Array_2/Max_ADU*FWC/QE;
% Mean_Array_2_p=Mean_Array_2/Max_ADU*FWC/QE;

Shot_1_p=Mean_Array_1_p.^0.5;
% Shot_2_p=Mean_Array_2_p.^0.5;

Std_Mean_Array_wo_Shot_1_p=(Std_Mean_Array_1_p.^2-Shot_1_p.^2).^0.5;
% Std_Mean_Array_wo_Shot_2_p=(Std_Mean_Array_2_p.^2-Shot_2_p.^2).^0.5;

%     plot(Power_Array,Std_Mean_Array_1);
% 
%     plot(Power_Array,Mean_Array_1);
%     plot(Power_Array,Std_Mean_Array_1_p,Power_Array,Std_Mean_Array_wo_Shot_1_p);

Matrix=[Power_Array; ones([1 length(Power_Array)])]';
beta=Matrix\Std_Mean_Array_wo_Shot_1_p;
Std_Mean_Array_wo_Shot_1_p_fit=beta(1)*Power_Array+beta(2);

plot(Power_Array,Std_Mean_Array_wo_Shot_1_p,Power_Array,Std_Mean_Array_wo_Shot_1_p_fit);
plot(Mean_Array_1,Std_Mean_Array_wo_Shot_1_p,Mean_Array_1,Std_Mean_Array_wo_Shot_1_p_fit);

%% Noise Fitting
Bin_Adj_Factor=1;

plot(Mean_Array_1,Std_Mean_Array_1);

g = fittype( @(a,b,c,x) (a*x.^2+b*x+c)); %a: excess, b: shot, c: dark
FitResult = fit(Mean_Array_1, Std_Mean_Array_1.^2, g, 'StartPoint', [0.015, 0.0569, -0]);
FitResult.a
FitResult.b
FitResult.c

% 從b推估FWC (assume QE=0.5)
FWC_estimated=QE*Max_ADU/FitResult.b;


Std_Mean_Array_without_shot_1=(Std_Mean_Array_1.^2-FitResult.b*Mean_Array_1).^0.5;
%Std_Mean_Array_fit=(FitResult.a*Mean_Array_1.^2+FitResult.b*Mean_Array_1+FitResult.c).^0.5;
Std_Mean_Array_fit=(FitResult.a*Mean_Array_1.^2+Mean_Array_1*(1/Max_ADU*FWC_estimated/QE)^-1+FitResult.c).^0.5;

subplot(1,1,1)
plot(Power_Array,Std_Mean_Array_1*Bin_Adj_Factor,Power_Array,Std_Mean_Array_fit*Bin_Adj_Factor,'-o');
xlabel('Reference Power (mW)');
ylabel('Noise (ADU)');
legend('Measured','Fitted','Location','Best');


plot(Mean_Array_1*Bin_Adj_Factor,Std_Mean_Array_1*Bin_Adj_Factor,Mean_Array_1*Bin_Adj_Factor,Std_Mean_Array_fit*Bin_Adj_Factor,'-o');
xlabel('Signal (ADU)');
ylabel('Noise (ADU)');
xlim([0 4096]);
ylim([0 20]);
legend('Measured',sprintf('Fitted, 2nd Coef=%g',FitResult.b),'Location','NorthWest');

Output=[Mean_Array_1 Std_Mean_Array_1 Std_Mean_Array_fit];

dlmwrite([Data_Save_Folder 'Output_Different Power.txt'],Output,'delimiter','\t','newline','pc','precision', '%.6f');

plot(Power_Array,Mean_Array_1,'-o');
xlabel('Input Power (Reference) (mW)');
ylabel('Signal (ADU)');



plot(Mean_Array_1*Bin_Adj_Factor,(Std_Mean_Array_1*Bin_Adj_Factor).^2,Mean_Array_1*Bin_Adj_Factor,(Std_Mean_Array_fit*Bin_Adj_Factor).^2,'-o');
xlabel('Signal (ADU)');
ylabel('Noise (ADU^2)');
xlim([0 4096]);
ylim([0 400]);
legend('Measured',sprintf('Fitted, 2nd Coef=%g',FitResult.b),'Location','NorthWest');
