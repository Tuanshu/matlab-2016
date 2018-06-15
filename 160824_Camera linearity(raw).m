clear all

%folder_path='C:\Users\Owner\Desktop\160824_Exp Test_0.125A\';
%folder_path='C:\Users\Owner\Desktop\160824_Exp Test_0.125A_defocus\';
folder_path='C:\Users\Owner\Desktop\160824_Exp Test_0.125A_defocus_LUT off\';

Exposure_time_Index=[1:20];
dET=50;
Exposure_time=(Exposure_time_Index-1)*dET;
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
    %Image_mean_Array(p)=mean(mean(Image(275:325,250:300)));
    Image_mean_Array(p)=mean(mean(Image(275:325,50:100)));
    Image_Intensity_Array(p)=double(Image(X,Y));
end

%%

Matrix=[Image_mean_Array' ones([length(Image_mean_Array) 1])];

beta=Exposure_time'\Matrix;

Image_mean_Array_fit=beta(1)*(Exposure_time)+beta(2);


%Image_mean_Array_Fit=interp1(Exposure_time,Image_mean_Array,[Exposure_time(1) Exposure_time(end)]);

plot(Exposure_time,Image_mean_Array,Exposure_time,Image_mean_Array_fit);
%plot(Exposure_time,Image_mean_Array,[Exposure_time(1) Exposure_time(end)],Image_mean_Array_Fit);
%%
imagesc(Image);


%%

dADU=0;


plot(Exposure_time,Image_mean_Array);
plot(Exposure_time,Image_mean_Array,Exposure_time,Image_mean_Array_fit);

%plot(Image_mean_Array,Image_mean_Array./(Exposure_time+dET));

%ylim([0.22 0.3]);
xlabel('Exposure_time (microsec)');
ylabel('Signal (ADU)');
legend('Measured','Fitted','Location','NorthWest');


%%
