clear all
%%
root_folder_path='F:\160905_Imperx Noise\';
last_folder_name='For After N-point Noise Calc_648x3';
folder_path=[root_folder_path last_folder_name '\'];

Data_Save_Folder='F:\Processed Data\';


cd(folder_path);

% Data format related
Row=648;
Colomn=3;
Byte_Skip=0;
% Processing related
ave_factor=[1 2 4 8 16 32 64];
Product_of_Axial_Decimation_Factor_and_Ave_Factor=64;

for QQQ=1:length(ave_factor)
    N=4;
    micron_per_frame=0.2/ave_factor(QQQ)/N;
    Offset_1=0;
    Offset_2=0; 

    Max_Number_of_Frame=2000;

    Processed_Data_Path=[Data_Save_Folder last_folder_name sprintf('_Ave_Factor_%d.raw',ave_factor(QQQ))];

    %%
    file_list=dir(folder_path);
    file_list=downsample(file_list(4:end),floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/ave_factor(QQQ)));
    Frame=length(file_list);
    %%
    Ave_Temp=zeros(Row,Colomn,ave_factor(QQQ));
    After_Npoint_Frame_Length=floor(Frame/N/ave_factor(QQQ));
    After_Npoint_Image_Stack=zeros(Row,Colomn,min(Max_Number_of_Frame,After_Npoint_Frame_Length));
    Npoint_Temp=zeros(Row,Colomn,N);

    X=[1:Frame];
    for p=1:min(Max_Number_of_Frame,After_Npoint_Frame_Length)
        for q=1:N
            for r=1:ave_factor(QQQ)

                file_path=[folder_path file_list((p-1)*N*ave_factor(QQQ)+(q-1)*ave_factor(QQQ)+r).name];

                fin=fopen(file_path);

                fseek(fin, Byte_Skip, 'bof');

                Ave_Temp(:,:,r)=fread(fin,[Row,Colomn],'uint16'); %*Frame   不知為何, 看起來就像是要除16
                fclose(fin);
            end
            Npoint_Temp(:,:,q)=mean(Ave_Temp,3);
        end
        After_Npoint_Image_Stack(:,:,p)=((N*sum(Npoint_Temp.^2,3)-sum(Npoint_Temp,3).^2).^0.5)*(2^0.5)/N;
        disp(p);
    end

%     %% Show Image
% 
%     NN=1;
% 
%     Image_Show=Ave_Temp(:,:,NN);
% 
%     imagesc(Image_Show);
%     colorbar
% 
%     %% Show Image
% 
%     NN=1;
% 
%     Image_Show=Npoint_Temp(:,:,NN);
% 
%     imagesc(Image_Show);
%     colorbar
%     %% Show Image
% 
%     NN=40;
% 
%     Image_Show=After_Npoint_Image_Stack(:,:,NN);
% 
%     imagesc(Image_Show);
%     caxis([0 40]);
%     colorbar
%     %% STD
% 
%     STD_Map=std(After_Npoint_Image_Stack,0,3);
%     Mean_Map=mean(After_Npoint_Image_Stack,3);
% 
%     Normalized_STD_Map=STD_Map./Mean_Map;
%     imagesc(STD_Map);
% 
%     imagesc(Normalized_STD_Map);
%     caxis([0 0.02]);
%     xlabel('Normalized STD');
%     %% Z-plot along certian pixel
%     X_Show=100;%Row/2+23;
%     Y_Show=100;%Colomn/2+108;
% 
%     Array_Show=squeeze(After_Npoint_Image_Stack(X_Show,Y_Show,:));
% 
%     plot(Array_Show);
% 
%     %% Z-plot along certian pixel
%     X_Show=size(After_Npoint_Image_Stack,1)/2;
%     Y_Show=size(After_Npoint_Image_Stack,1)/2;
% 
%     Array_Show=squeeze(After_Npoint_Image_Stack(X_Show,Y_Show,:));
% 
%     plot(Array_Show);
% 
% 
%     %%
%     NN=74;
%     Max=16;
%     Min=6;
% 
% 
%     Image_Show(:)=After_Npoint_Image_Stack(:,:,NN);
% 
%     Image_Show_norm=(Image_Show-Min)/(Max-Min);
% 
%     Image_Show_norm(Image_Show_norm>1)=1;
%     Image_Show_norm(Image_Show_norm<0)=0;
% 
%     imagesc(Image_Show_norm);
%     colormap(gray);
    %%
    fid = fopen(Processed_Data_Path, 'w+');
    fwrite(fid, After_Npoint_Image_Stack, 'double');
    fclose(fid);
    disp(QQQ);
end
%%