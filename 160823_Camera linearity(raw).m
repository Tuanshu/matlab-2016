clear all

folder_path='C:\Users\Owner\Desktop\160823_Exposure time test\';

Exposure_time_Index=[0:30];
dET=100;
Exposure_time=Exposure_time_Index*dET;
X=100;
Y=100;

Row=648;
Colomn=488;

for p=1:length(Exposure_time_Index)
    File_name=sprintf('%d.raw',Exposure_time_Index(p));
    fin=fopen([folder_path File_name]);
    Image=fread(fin,[Row,Colomn],'uint16'); %*Frame   不知為何, 看起來就像是要除16

    %Image=imread([folder_path File_name]);
    %Image_mean_Array(p)=mean(mean(Image(250:350,225:325)));
    Image_mean_Array(p)=mean(mean(Image(275:325,250:300)));
    Image_Intensity_Array(p)=double(Image(X,Y));
end

%%



plot(Exposure_time,Image_mean_Array);

Image_mean_Array_Fit=interp1(Exposure_time,Image_mean_Array,[Exposure_time(1) Exposure_time(end)]);


plot(Exposure_time,Image_mean_Array,[Exposure_time(1) Exposure_time(end)],Image_mean_Array_Fit);
%%
imagesc(Image);


%%

dADU=0;


plot(Exposure_time,Image_mean_Array);

%plot(Image_mean_Array,Image_mean_Array./(Exposure_time+dET));

%ylim([0.22 0.3]);



%%
