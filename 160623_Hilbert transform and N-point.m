clear all

cd('D:\160608_PSF_Oil and vibration\');

Set=2;     %1 for FW4, 2 for FW5

FW4=dlmread('FW4.txt'); 
FW4=FW4(:,2);
FW5=dlmread('FW5.txt'); 
FW5=FW5(:,2);

X=1:2000;

plot(FW4);
plot(FW5);

if Set == 1
    [value index]=max(FW4);
    Data=FW4(index-1000:index+999);
elseif Set == 2
    [value index]=max(FW5);
    
    Data=FW5(index-1000:index+999);
end

plot(Data);
xlabel('Frame #');
ylabel('LSB');
%% Hilbert
SP=20;
LP=40;

Data_spectral=fft(Data);

plot(real(Data_spectral));
xlabel('Frequency (bin)');
ylabel('Spectral Signal')
xlim([0 60]);
ylim([-9000 9000]);
Data_spectral_filtered=Data_spectral;

Data_spectral_filtered(1:SP)=0;
Data_spectral_filtered(LP:end)=0;

plot(real(Data_spectral_filtered));
xlabel('Frequency (bin)');
ylabel('Spectral Signal')
xlim([0 60]);
ylim([-9000 9000]);

Data_filtered=ifft(Data_spectral_filtered)*2;
Data_filtered_car=real(Data_filtered);
Data_filtered_env=abs(Data_filtered);


plot(X,Data_filtered_car,X,Data_filtered_env);

%% N-point
ave_factor=19;
micron_per_frame=0.2/ave_factor/4;
N=4;
Offset_1=1;
Offset_2=1;

Data_ave=zeros([round(length(Data)/ave_factor)-1 1]);
X_ave=zeros([round(length(Data)/ave_factor)-1 1]);


for p=1:length(Data_ave)
    Data_ave(p)=mean(Data((ave_factor*(p-1)+1+Offset_1):(ave_factor*p+Offset_1)));
    X_ave(p)=mean(X((ave_factor*(p-1)+1+Offset_1):(ave_factor*p+Offset_1)));

end

plot(X,Data,X_ave,Data_ave)

Data_Npoint=zeros([round(length(Data_ave)/N)-1 1]);
X_Npoint=zeros([round(length(Data_ave)/N)-1 1]);



for p=1:length(Data_Npoint)
    subarray=(Data_ave((N*(p-1)+1+Offset_2):(N*p+Offset_2)));
    
    Data_Npoint(p)=((N*sum(subarray.^2)-sum(subarray)^2).^0.5)*(2^0.5)/N;
    X_Npoint(p)=mean(X_ave((N*(p-1)+1+Offset_2):(N*p+Offset_2)));

end

plot(X,Data-mean(Data_ave),X_ave,Data_ave-mean(Data_ave),X_Npoint,Data_Npoint);

plot(X,Data_filtered_car,X,Data_filtered_env,X_Npoint,Data_Npoint);


plot((X-1000)*micron_per_frame,Data_filtered_car/max(Data_filtered_car),(X-1000)*micron_per_frame,Data_filtered_env/max(Data_filtered_car),(X_Npoint-1000)*micron_per_frame,Data_Npoint/max(Data_filtered_car));
plot((X-1000)*micron_per_frame,Data_filtered_car/max(Data_filtered_car),(X-1000)*micron_per_frame,Data_filtered_env/max(Data_filtered_car));
xlim([-2.5 2.5]);
ylim([-1 1]);
xlabel('Position (\mum)');
ylabel('Signal (a.u.)');


Data_filtered_env_norm=Data_filtered_env/max(Data_filtered_env);

Data_Npoint_norm=Data_Npoint/max(Data_Npoint);

Data_Npoint_norm_interp=interp1(X_Npoint,Data_Npoint_norm,X);

FWHM_Hilbert=X(find(Data_filtered_env_norm>0.5,1,'last')-find(Data_filtered_env_norm>0.5,1,'first'))*micron_per_frame;
FWHM_Np=X(find(Data_Npoint_norm_interp>0.5,1,'last')-find(Data_Npoint_norm_interp>0.5,1,'first'))*micron_per_frame;

Ratio=FWHM_Np/FWHM_Hilbert;

plot((X-1000)*micron_per_frame,Data_filtered_env_norm,(X-1000)*micron_per_frame,Data_Npoint_norm_interp);
xlim([-2.5 2.5]);
ylim([0 1]);
xlabel('Position (\mum)');
ylabel('Signal (norm)');