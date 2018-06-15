clear all

folder_path='C:\Users\Owner\Desktop\160823_Exposure time\160803_ADU vs exposure time\';

Exposure_time=[2 4 8 16 32 64 128 256 512];

X=100;
Y=100;

for p=1:length(Exposure_time)
    File_name=sprintf('it%d.bmp',Exposure_time(p));
    Image=imread([folder_path File_name]);
    Image_mean_Array(p)=mean(Image(:));
    Image_Intensity_Array(p)=double(Image(X,Y));
end

plot(Exposure_time,Image_mean_Array);

Image_mean_Array_Fit=interp1(Exposure_time,Image_mean_Array,[Exposure_time(1) Exposure_time(end)]);


plot(Exposure_time,Image_mean_Array,[Exposure_time(1) Exposure_time(end)],Image_mean_Array_Fit);
%%
dADU=-2;


plot(Exposure_time,(Image_mean_Array+dADU)./Exposure_time);

%plot(Image_mean_Array,Image_mean_Array./(Exposure_time+dET));

xlim([50 550]);
%ylim([0.22 0.3]);



%%
plot(Exposure_time,Image_Intensity_Array./Exposure_time);
