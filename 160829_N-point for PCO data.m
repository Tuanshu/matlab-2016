clear all
%%
root_folder_path='F:\P1.2 Test\';
last_folder_name='160829_PCO_4160fps_4x_200microsec_forearm_no scan_with glycerol_for BND';
folder_path=[root_folder_path last_folder_name '\'];

Data_Save_Folder='F:\P1.2 Test\Processed Data\';

Processed_Data_Path=[Data_Save_Folder last_folder_name '.bin'];

cd(folder_path);

% Data format related
Row=1000;
Colomn=1000;
Byte_Skip=1024;
% Processing related
ave_factor=64;
N=4;
micron_per_frame=0.2/ave_factor/N;
Offset_1=0;
Offset_2=0; 

Start_Frame=0;


Max_Number_of_Frame=1000;


%%
file_list=dir(folder_path);
file_list=file_list(4:end);
Frame=length(file_list);
%%
Ave_Frame_Length=floor(Frame/ave_factor);
Ave_Image_Stack=zeros(Row,Colomn,min(Ave_Frame_Length,Max_Number_of_Frame));
Ave_Temp=zeros(Row,Colomn,ave_factor);
X=[1:Frame];
for p=1:min(Ave_Frame_Length,Max_Number_of_Frame)
    for q=1:ave_factor

        file_path=[folder_path file_list((p-1+Start_Frame)*ave_factor+q+Offset_1).name];

        fin=fopen(file_path);

        fseek(fin, Byte_Skip, 'bof');

        Ave_Temp(:,:,q)=fread(fin,[Row,inf],'uint16')/16; %*Frame   不知為何, 看起來就像是要除16
        fclose(fin);
    end
    Ave_Image_Stack(:,:,p)=mean(Ave_Temp,3);
    disp(p);
    X_ave(p)=mean(X(((p-1)*ave_factor+1+Offset_1):((p)*ave_factor+Offset_1)));
end

%% Show Image

NN=200;

Image_Show=Ave_Image_Stack(:,:,NN);

imagesc(Image_Show);
%% STD

STD_Map=std(Ave_Image_Stack,0,3);
Mean_Map=mean(Ave_Image_Stack,3);

Normalized_STD_Map=STD_Map./Mean_Map;
imagesc(STD_Map);

imagesc(Normalized_STD_Map);
caxis([0 0.02]);
xlabel('Normalized STD');
%% Z-plot along certian pixel
X_Show=470;%Row/2+23;
Y_Show=717;%Colomn/2+108;

Array_Show=squeeze(Ave_Image_Stack(X_Show,Y_Show,:));

plot(Array_Show);


%% N-point

Ave_Image_Stack_Npoint=zeros([Row Colomn floor(min(Ave_Frame_Length,Max_Number_of_Frame)/N)-1]);
X_Npoint=zeros([floor(min(Ave_Frame_Length,Max_Number_of_Frame)/N)-1 1]);

for p=1:size(Ave_Image_Stack_Npoint,3)
    subarray=(Ave_Image_Stack(:,:,(N*(p-1)+1+Offset_2):(N*p+Offset_2)));
    Ave_Image_Stack_Npoint(:,:,p)=((N*sum(subarray.^2,3)-sum(subarray,3).^2).^0.5)*(2^0.5)/N;
    X_Npoint(p)=mean(X_ave((N*(p-1)+1+Offset_2):(N*p+Offset_2)));
    disp(p);
end

%% Z-plot along certian pixel
X_Show=size(Ave_Image_Stack_Npoint,1)/2;
Y_Show=size(Ave_Image_Stack_Npoint,1)/2;

Array_Show=squeeze(Ave_Image_Stack_Npoint(X_Show,Y_Show,:));

plot(Array_Show);


%%
NN=115;
Max=16;
Min=6;


Image_Show(:)=Ave_Image_Stack_Npoint(:,:,NN);

Image_Show_norm=(Image_Show-Min)/(Max-Min);

Image_Show_norm(Image_Show_norm>1)=1;
Image_Show_norm(Image_Show_norm<0)=0;

imagesc(Image_Show_norm);
colormap(gray);
%%
fid = fopen(Processed_Data_Path, 'w+');
fwrite(fid, Ave_Image_Stack_Npoint, 'double');
fclose(fid);

%%
Data_Npoint=squeeze(Ave_Image_Stack_Npoint(size(Ave_Image_Stack_Npoint,1)/2,size(Ave_Image_Stack_Npoint,1)/2,:));
Data_Npoint_norm=Data_Npoint/max(Data_Npoint);
Data_Npoint_norm_interp=interp1(X_Npoint,Data_Npoint_norm,X);

FWHM_Np=X(find(Data_Npoint_norm_interp>0.5,1,'last')-find(Data_Npoint_norm_interp>0.5,1,'first'))*micron_per_frame;
