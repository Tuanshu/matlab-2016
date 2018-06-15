clear all

Sampling_Rate=200000;   %Hz

Data_Number=8;


Data=dlmread(sprintf('D:\\161124_Closed-loop PZT analysis\\161125_PZT test by LabVIEW_%d.txt',Data_Number),'\t');

if Data_Number==6
    Data_Range=[4.2E5:9.9E5];
elseif Data_Number==8
    Data_Range=[3.2E5:14.9E5];
elseif Data_Number==9
    Data_Range=[4.4E5:5.6E5];
else    
    Data_Range=[1:length(Data)];
end


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
Half_Wavelength=0.315;

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
Position_Sampling_Resolution=0.78/2/4;

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

%% Velocity Profile based on reduced data

Velocity_Reduced=[0;diff(Position_Reduced)./diff(Time_Reduced)];

plot(Time_Reduced,Velocity_Reduced);

figure('Position', [100, 100, 800, 500]);
plot(Position_Reduced,Velocity_Reduced);
%title(sprintf('Data Number=%d',Data_Number),'fontsize',15,'color','b');
xlabel('Position (micron)','fontsize',15)
ylabel('Velocity (micron/sec)','fontsize',15)
xlim([100 200]);
set(gca,'fontsize',15)

figure('Position', [100, 100, 800, 500]);
plot(Time_Reduced,Velocity_Reduced);
%title(sprintf('Data Number=%d',Data_Number),'fontsize',15,'color','b');
xlabel('Time (sec)','fontsize',15)
ylabel('Velocity (micron/sec)','fontsize',15)
xlim([1 2]);
set(gca,'fontsize',15)

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

plot(Time,Signal_Sub,Time_Reduced_Npoint,Signal_Reduced_Npoint);
