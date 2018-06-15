clear all
%%

Frame_Rate=10000;             %Hz
Time_Spacing=1/Frame_Rate;  %sec
Frequency_Max=Frame_Rate;

Save_File_Name='161012_TiSa test_Mikrotron';
root_folder_path='I:\Everday Experiements\161021_Power Stability Issue\Mikrotron\';
last_folder_name='Ce_Images_4';
folder_path=[root_folder_path last_folder_name '\'];

If_EffMap=-1;
If_Norm=0;
Data_Save_Folder='I:\Processed Data\';
MIN=100;

cd(folder_path);

% Data format related
Row=160;
Colomn=4;
Byte_Skip=0;
% Processing related
ave_factor=[1];
Product_of_Axial_Decimation_Factor_and_Ave_Factor=1;

for QQQ=1:length(ave_factor)
    N=1;
    micron_per_frame=0.2/ave_factor(QQQ)/N;
    Offset_1=0;
    Offset_2=0; 

    Max_Number_of_Frame=3000000;
    if If_EffMap ==1
        Processed_Data_Path=[Data_Save_Folder Save_File_Name last_folder_name sprintf('_EffMap_Ave_Factor_%d.raw',ave_factor(QQQ))];
    elseif If_EffMap ==0
        Processed_Data_Path=[Data_Save_Folder Save_File_Name last_folder_name sprintf('_Ave_Factor_%d.raw',ave_factor(QQQ))];
    elseif If_EffMap ==-1
        Processed_Data_Path=[Data_Save_Folder Save_File_Name last_folder_name sprintf('Test_Ave_Factor_%d.raw',ave_factor(QQQ))];
    end
   %%
    file_list=dir(folder_path);
    file_list=file_list(4:end);
  %% fiel name sorting
    name = {file_list.name};
    str  = sprintf('%s ', name{:});
    num  = sscanf(str, '%u.bmp ');
    [dummy, index] = sort(num);
    file_list = file_list(index);
%%
    file_list=downsample(file_list,floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/ave_factor(QQQ)));
    Frame=length(file_list);
    %%
    Ave_Temp=zeros(Row,Colomn,ave_factor(QQQ));
    After_Npoint_Frame_Length=floor(Frame/N/ave_factor(QQQ));
    After_Npoint_Image_Stack=zeros(Row,Colomn,min(Max_Number_of_Frame,After_Npoint_Frame_Length));
    Npoint_Temp=zeros(Row,Colomn,N);

    X=[1:Frame];
    Frame_Ave=zeros([min(Max_Number_of_Frame,After_Npoint_Frame_Length)*N 1]);
    for p=1:min(Max_Number_of_Frame,After_Npoint_Frame_Length)
        for q=1:N
            for r=1:ave_factor(QQQ)

                file_path=[folder_path file_list((p-1)*N*ave_factor(QQQ)+(q-1)*ave_factor(QQQ)+r).name];

                Ave_Temp(:,:,r)=double(imread(file_path,'TIFF'))'*16;
            end
            Npoint_Temp(:,:,q)=mean(Ave_Temp,3);
            Frame_Ave((p-1)*N+q)=mean(mean(Npoint_Temp(:,:,q)));

            if If_Norm ==1
                Npoint_Temp(:,:,q)=(Npoint_Temp(:,:,q)-MIN)/(Frame_Ave((p-1)*N+q)-MIN);
            end
        end
        
        if If_EffMap ==1
            After_Npoint_Image_Stack(:,:,p)=((N*sum(Npoint_Temp.^2,3)-sum(Npoint_Temp,3).^2).^0.5)*(2^0.5)/N./mean(Npoint_Temp,3);
        elseif If_EffMap ==0
            After_Npoint_Image_Stack(:,:,p)=((N*sum(Npoint_Temp.^2,3)-sum(Npoint_Temp,3).^2).^0.5)*(2^0.5)/N;
        elseif If_EffMap ==-1
            After_Npoint_Image_Stack(:,:,p)=mean(Npoint_Temp,3);
        end
        
        disp(p);
    end

    if If_Norm ==1
        After_Npoint_Image_Stack=After_Npoint_Image_Stack*mean(Frame_Ave);
    end
    %
    fid = fopen(Processed_Data_Path, 'w+');
    fwrite(fid, After_Npoint_Image_Stack, 'double');
    fclose(fid);
    disp(QQQ);
    
    %%
    Axial_ave_Factor=1;
    Maximum_Axial_Frame=200000;
    Reduced_Stack=0;
    Axial_Length_Original=size(After_Npoint_Image_Stack,3);
    Axial_Length_Used=floor(Axial_Length_Original/Axial_ave_Factor)*Axial_ave_Factor;
    Reduced_Length=Axial_Length_Used/Axial_ave_Factor;
    for p=1:Axial_ave_Factor
       Reduced_Stack=Reduced_Stack+After_Npoint_Image_Stack(:,:,(Axial_ave_Factor-(p-1)):Axial_ave_Factor:(Axial_ave_Factor*Reduced_Length)-(p-1));
    end
    Reduced_Stack=Reduced_Stack/Axial_ave_Factor;
    
    clear After_Npoint_Image_Stack

    
    %%
    
    Start_Y_Index=1;
    Y_Thickness=4;
    
    if If_EffMap==-1
        Reduced_Image=squeeze(Reduced_Stack(:,1,:))';
    else
        Reduced_Image=squeeze(mean(Reduced_Stack(:,(Start_Y_Index:(Start_Y_Index+Y_Thickness-1)),:),2))';
    end

    
    Reduced_Image=Reduced_Image(1:min(size(Reduced_Image,1),Maximum_Axial_Frame),:);
%% To generate noise spectrum
  %%
Reduced_Stack(end,:,:)=[];
X_Plot=50;
Y_Plot=4;
Z_Profile=squeeze(Reduced_Stack(X_Plot,Y_Plot,:));
plot(Z_Profile);
STD_Map=std(Reduced_Stack,0,3);
imagesc(STD_Map);
Reduced_Line=squeeze(mean(mean(Reduced_Stack,1),2));
STD=std(Reduced_Line);
Noise_Spectrum=fft(Reduced_Line);
Frequency_Spacing=Frequency_Max/length(Noise_Spectrum);
Frequency_Array=Frequency_Spacing*[1:length(Noise_Spectrum)];
Noise_Spectrum_Norm=abs(Noise_Spectrum)/5E5;
figure('Position', [100, 100, 800, 500]);
plot(Frequency_Array,Noise_Spectrum_Norm,'LineWidth',1.5);
xlim([0 1000]);
ylim([0 0.1]);
xlabel('Frequency (Hz)','fontsize',15);
ylabel('Normalized Spectral Power (abs.)','fontsize',15);
    
    % IE map
    
    %[Max_Map Max_Index_Map]=max(Reduced_Stack,[],3);
    %Mean_Map=mean(Reduced_Stack,3);
    
    %Mean_Map=mean(Reduced_Stack(:,:,(Max_Index_Map-500):(Max_Index_Map+500)),3);

    %IE_Map=(Max_Map-Mean_Map)./Mean_Map;
    %imagesc(IE_Map);
    
    %
    
    ROI_Depth_Glass=[1 5];
    ROI_Width_Glass=[200 250];
    
    ROI_Depth_Sig=[111 115];
    ROI_Width_Sig=[200 250];
  
    ROI_Depth_BND=[21 25];
    ROI_Width_BND=[200 250];
    
    if If_EffMap ==1
        C_max=0.01;
        C_min=0.001;
        Reduced_Image_normalized=(Reduced_Image-C_min)/(C_max-C_min);

        Reduced_Image_normalized(Reduced_Image_normalized<0)=0;
        Reduced_Image_normalized(Reduced_Image_normalized>1)=1;

    elseif If_EffMap ==0
        C_max=32;
        C_min=3.2;
        Reduced_Image_normalized=(Reduced_Image-C_min)/(C_max-C_min);

        Reduced_Image_normalized(Reduced_Image_normalized<0)=0;
        Reduced_Image_normalized(Reduced_Image_normalized>1)=1;
    elseif If_EffMap ==-1
        C_max=4096;
        C_min=0;
        Reduced_Image_normalized=Reduced_Image;

    end
    
%     subplot(1,1,1)
%     imagesc(Reduced_Image_normalized);
%     %caxis([0 1]);
%     axis equal
%     axis off
%     set(gca,'xtick',[0  20  40  60  80  100],'xticklabel',[0  20  40  60  80  100]*0.33,'ytick',[0  100  200  300  400  500 600],'yticklabel',[0  100  200  300  400  500 600]*0.2645);
%     xlabel('Lateral Position (micron)');
%     ylabel('Axial Position (micron)');
%     xlim([1 size(Reduced_Image_normalized,2)])
%     ylim([1 size(Reduced_Image_normalized,1)])
%     colormap(gray);
    %%
    Reduced_Line_normalized=Reduced_Image_normalized(:,size(Reduced_Image_normalized,2)/2);
    dlmwrite([Data_Save_Folder 'PSF_' last_folder_name '.txt'],Reduced_Line_normalized,'delimiter','\t','newline','pc','precision', '%.6f');

    %%plot(Reduced_Line_normalized);
    if If_EffMap ==1
        imwrite(Reduced_Image_normalized,[Data_Save_Folder Save_File_Name last_folder_name '_Bscan_EffMap.png'],'png');
    elseif If_EffMap ==0
        imwrite(Reduced_Image_normalized,[Data_Save_Folder Save_File_Name last_folder_name '_Bscan.png'],'png');
    elseif If_EffMap ==-1
        imwrite(Reduced_Image_normalized,[Data_Save_Folder Save_File_Name last_folder_name '_Bscan_Test.png'],'png');
     end
%     Eff_Coef_Glass(QQQ)=max(mean(Reduced_Image(ROI_Depth_Glass(1):ROI_Depth_Glass(2),ROI_Width_Glass(1):ROI_Width_Glass(2))))
%     Eff_Coef_Sig(QQQ)=max(mean(Reduced_Image(ROI_Depth_Sig(1):ROI_Depth_Sig(2),ROI_Width_Sig(1):ROI_Width_Sig(2))))
%     Eff_Coef_BND(QQQ)=max(mean(Reduced_Image(ROI_Depth_BND(1):ROI_Depth_BND(2),ROI_Width_BND(1):ROI_Width_BND(2))))
%     
%     Eff_Coef_Glass_BND_Sub(QQQ)=max(mean(((Reduced_Image(ROI_Depth_Glass(1):ROI_Depth_Glass(2),ROI_Width_Glass(1):ROI_Width_Glass(2))).^2-Eff_Coef_BND(QQQ)^2).^0.5));
%     Eff_Coef_Sig_BND_Sub(QQQ)=max(mean(((Reduced_Image(ROI_Depth_Sig(1):ROI_Depth_Sig(2),ROI_Width_Sig(1):ROI_Width_Sig(2))).^2-Eff_Coef_BND(QQQ)).^0.5));

%% For Noise Analysis
    ROI_Width=[1 8];
    ROI_Height=[1 8];
    Noise_ROI=Reduced_Image(ROI_Height(1):ROI_Height(2),ROI_Width(1):ROI_Width(2));
    Noise_EffMap=mean(Noise_ROI(:))
    %imagesc(Reduced_Image);

    %colormap(gray);
%% IE calc   
    
 end
% subplot(3,1,1)
% plot(ave_factor,(Eff_Coef_Glass.^2-Eff_Coef_BND.^2).^0.5);
% plot(ave_factor,Eff_Coef_Glass_BND_Sub);
% ylim([0 0.02]);
% 
% subplot(3,1,2)
% plot(ave_factor,(Eff_Coef_Sig.^2-Eff_Coef_BND.^2).^0.5);
% ylim([0 0.01]);
% 
% subplot(3,1,3)
% plot(ave_factor,Eff_Coef_Sig-Eff_Coef_BND);
% %plot(ave_factor,Eff_Coef_BND);
% ylim([0 0.01]);

fclose('all');

