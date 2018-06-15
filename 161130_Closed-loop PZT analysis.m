clear all

Sampling_Rate=200000;   %Hz

Data_Number=5;

cd('D:\161129_PZT test_PID\');
Data=dlmread(sprintf('D:\\161129_PZT test_PID\\161129_1st PID test_Setting %d.txt',Data_Number),'\t');

TH=7;   %擷取從第一個超過TH的位置開始, 在最後一個超過TH的位置結束

Start_Index=find(Data>TH,1,'first');
End_Index=12*1E5;%find(Data>TH,1,'last');
Data_Range=[Start_Index:End_Index];


Signal=Data(Data_Range);
plot(Signal)
clear Data

%%
Time=(([1:length(Signal)]-1)/Sampling_Rate)';

plot(Time,Signal);
xlabel('Time (second)');
ylabel('Signal');

Signal_Sub=Signal-mean(Signal);
%% Use Running Ave to sub DC (Fail)
% DC_Window_Size=1000;
% 
% 
% Temp_Signal_DC=0;
% Reduced_Length=length(Signal)/DC_Window_Size;
% for p=1:DC_Window_Size
%     Temp_Signal_DC=Temp_Signal_DC+Signal((DC_Window_Size-(p-1)):DC_Window_Size:(DC_Window_Size*Reduced_Length)-(p-1));
% end
% Signal_DC=Temp_Signal_DC/DC_Window_Size;
% 
% 
%     
% Signal_Sub=Signal-Signal_DC;
% plot(Time,Signal_Sub);
% 
% clear Signal_DC

%%
FFT=fft(Signal_Sub);
plot(real(FFT))
%%
SP=1;
LP=5000;
Half_Wavelength=0.633/2;

FFT_Filtered=FFT;

FFT_Filtered(1:SP)=0;

FFT_Filtered(LP:end)=0;

Signal_Reconstructed=ifft(FFT_Filtered)*2;
Phase=unwrap(angle(Signal_Reconstructed));
Phase=Phase-min(Phase);
Phase_Raw=(angle(Signal_Reconstructed));
plot(Time,Signal_Sub,Time,Signal_Reconstructed);
%% Position Waveform
plot(Time,Signal_Sub/max(Signal_Sub),Time,Phase_Raw/max(Phase_Raw));

Position=Phase*Half_Wavelength/pi/2;

plot(Time,Position);

%% Velocity Profile
Velocity=[0;diff(Position)./diff(Time)];

plot(Time,Velocity);
plot(Position,Velocity);

%% Reduced Position Sampling Resolution
Position_Sampling_Resolution=0.633/2/4;

Ave_d_Position=(max(Position)-min(Position))/length(Position);
Window_Size=round(Position_Sampling_Resolution/Ave_d_Position);

Temp_Position=0;
Temp_Time=0;
Temp_Signal=0;
Reduced_Length=length(Position)/Window_Size;
for p=1:Window_Size
    Temp_Time=Temp_Time+Time((Window_Size-(p-1)):Window_Size:(Window_Size*Reduced_Length)-(p-1));
    Temp_Position=Temp_Position+Position((Window_Size-(p-1)):Window_Size:(Window_Size*Reduced_Length)-(p-1));
    Temp_Signal=Temp_Signal+Signal((Window_Size-(p-1)):Window_Size:(Window_Size*Reduced_Length)-(p-1));
end
Position_Reduced=Temp_Position/Window_Size;
Time_Reduced=Temp_Time/Window_Size;
Signal_Reduced=Temp_Signal/Window_Size;

plot(Time,Position,Time_Reduced,Position_Reduced);
xlabel('Time (second)','fontsize',15)
ylabel('Position (micron)','fontsize',15)
%% Velocity Profile based on reduced data
Position_Range=100+[0 100];

Position_Linear=interp1([Time_Reduced(1) Time_Reduced(end)],[Position_Reduced(1) Position_Reduced(end)],Time_Reduced);
Position_Reduced_Sub=Position_Reduced-Position_Linear;

Velocity_Reduced=[0;diff(Position_Reduced)./diff(Time_Reduced)];

plot(Time_Reduced,Velocity_Reduced);

figure('Position', [100, 100, 800, 500]);
title(sprintf('Data Number=%d',Data_Number),'fontsize',15,'color','b');

subplot(2,1,1)
plot(Time_Reduced,Position_Reduced_Sub);
title(sprintf('Data Number=%d',Data_Number),'fontsize',15,'color','b');
xlabel('Time (second)','fontsize',15)
ylabel('Residual Position (micron)','fontsize',15)
xlim([1 2]);
%ylim([1 2]);
set(gca,'fontsize',15)

subplot(2,1,2)
plot(Position_Reduced,Velocity_Reduced);
xlabel('Position (micron)','fontsize',15)
ylabel('Velocity (micron/sec)','fontsize',15)
xlim([Position_Range(1) Position_Range(2)]);
ylim([0 200]);
set(gca,'fontsize',15)

%figure('Position', [100, 100, 800, 500]);
% subplot(3,1,3)
% plot(Time_Reduced,Velocity_Reduced);
% %title(sprintf('Data Number=%d',Data_Number),'fontsize',15,'color','b');
% xlabel('Time (sec)','fontsize',15)
% ylabel('Velocity (micron/sec)','fontsize',15)
% xlim([1 2]);
% ylim([0 200]);
% set(gca,'fontsize',15)

%% N-point test
N=4;
Signal_Reduced_Npoint=zeros(floor(length(Signal_Reduced)/N),1);
Time_Reduced_Npoint=zeros(floor(length(Time_Reduced)/N),1);
for p=1:length(Signal_Reduced_Npoint)
    Temp=Signal_Reduced((1+(p-1)*N):p*N);
    Temp_Time_2=Time_Reduced((1+(p-1)*N):p*N);
    Signal_Reduced_Npoint(p)=((N*sum(Temp.^2)-sum(Temp).^2).^0.5)*(2^0.5)/N;
    Time_Reduced_Npoint(p)=mean(Temp_Time_2);
end

% plot(Time,Signal_Sub,Time_Reduced_Npoint,Signal_Reduced_Npoint);
