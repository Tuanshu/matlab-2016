clear all

%%% To find the phase from min voltage to max voltage ONLY
%%%突然發現我該紀錄的好像應該是給的電壓 而不是回傳的電壓

Path='D:\AMO\160108_PZT (P4)\';

Parameter='Frequency';

Order=0;
Frequency=[0.004];
if Frequency==0.004
    N_Waveform=1000000/4;%400000;%1000000;  %Full period
elseif Frequency==0.01
    N_Waveform=400000;%400000;%1000000;  %Full period
end
%Frequency=[0.005];
Sampling_rate=500; %Hz
smooth_window=10000;
Range=10;        %V
min_voltage=(10-Range)/2;
max_voltage=min_voltage+Range;
dvoltage_min=0.1;
dvoltage_max=0.1;
Wavelength_Laser=0.6328;     %micron

Lower_Band=20;
Upper_Band=600;%fix(length(Spectrum)/2);

dt=1/Sampling_rate;
cd(Path);
clf
%gcf: handle of current fiqure
%gca: handle of current axes
%%%rect = [left, bottom, width, height]
set(gcf,'Position',get(0,'ScreenSize')+[200 100 -1000 -600])             %0可能是螢幕, gca是目前圖形的handle (如果沒有圖的話會自己開新的一個handle)
set(gca,'FontSize',16);

Time_length_reduction_factor=1;

required_number_of_loop=100;

Data=dlmread(sprintf('130103_f%g_V%d_O%d.txt',Frequency,Range,Order));                   %Asuume: line 1: time, line 2: interference signal, line 3: voltage 121219_f0.001_V2to8.txt怪怪

Voltage_all=Data(:,1);
Signal_all=Data(:,2);

clear Data

Scaling_Factor=1.0027;

Voltage_all=Voltage_all*Scaling_Factor;           %%normalization and scaling!

Estimated_Number_of_Pixel_for_One_Trip=1/Frequency/2*Sampling_rate;
Voltage_UniformBase=min_voltage:(max_voltage-min_voltage)/(Estimated_Number_of_Pixel_for_One_Trip-1):max_voltage;
Voltage_UniformBase=0:(max(Voltage_all))/(Estimated_Number_of_Pixel_for_One_Trip-1):max(Voltage_all);

p=1;
Voltage_previous_max_index=1;
Voltage_min_index=1;
Voltage_next_max_index=1;

while Voltage_previous_max_index<size(Signal_all,1)-2*(Voltage_next_max_index-Voltage_min_index)

    Voltage_min_index_temp=Voltage_previous_max_index+find(Voltage_all((Voltage_previous_max_index+1):end)<(min(Voltage_all)+dvoltage_min),1,'first');
    Voltage_next_max_index=Voltage_min_index_temp+find(Voltage_all((Voltage_min_index_temp+1):end)>(max(Voltage_all)-dvoltage_max),1,'first');
    [value Voltage_min_index]=min(Voltage_all((Voltage_min_index_temp+1):Voltage_next_max_index));
    Voltage_min_index=Voltage_min_index_temp+Voltage_min_index;
    Voltage_previous_max_index=Voltage_next_max_index;
    
    Voltage=Voltage_all(Voltage_min_index:Voltage_next_max_index);
    Signal=Signal_all(Voltage_min_index:Voltage_next_max_index);

    Time_index=1:length(Signal);

    Spectrum=ifft(Signal);
    Spectrum(1:Lower_Band)=0;
    Spectrum(Upper_Band:end)=0;

    Signal_New=fft(Spectrum,[],1);
    Phase_original=angle(Signal_New);
    Phase=unwrap(Phase_original);
    Phase_5V=Phase_original(find(Voltage>5,1,'first'));
    if p==1
        Phase_5V_0=Phase_5V;
    end
    Phase_5V_relative=Phase_5V-Phase_5V_0;
    Position_5V_relative=-1*Wavelength_Laser*Phase_5V_relative/(2*pi)/2;                %
    Position=-1*Wavelength_Laser*Phase/(2*pi)/2;
    %Position=Position-Position(find(Voltage>5,1,'first')); 因為反正要做內插,在內插後做這件事比較好
    Velocity=diff(Position);
    Velocity(length(Velocity)+1)=Velocity(length(Velocity));
    Velocity=smooth(Position,smooth_window);
    Position=smooth(Position,smooth_window);
    Voltage=smooth(Voltage,smooth_window);
   
    %plot(Voltage,Position);
    
    Velocity_UniformBase_temp=interp1(Voltage,Velocity,Voltage_UniformBase);
    Position_UniformBase_temp=interp1(Voltage,Position,Voltage_UniformBase);
    Velocity_UniformBase(:,p)=Velocity_UniformBase_temp;
    Position_UniformBase(:,p)=Position_UniformBase_temp-Position_UniformBase_temp(find(Voltage_UniformBase>=5,1,'first'))+Position_5V_relative; %其實這裡有點近似, -的是內插後,+卻是加內插前
    disp(p);
    p=p+1;
    
end
%clear Voltage_all Signal_all
q=1;
while q<size(Position_UniformBase,2)
    if q==1
        phase_temp=Position_UniformBase(1,q);
    elseif abs(Position_UniformBase(1,q)-phase_temp)>2
        Position_UniformBase(:,q)=[];
        q=q+1;
    end
    q=q+1;
end


plot(Voltage_UniformBase,mean(Position_UniformBase,2));
xlabel('Driving Voltage (V)','FontSize',18);
ylabel('PZT Position (micron)','FontSize',18);

plot(Voltage_UniformBase,std(Position_UniformBase,0,2));
xlabel('Driving Voltage (V)','FontSize',18);
ylabel('SD of PZT Position (micron)','FontSize',18);

%plot(Voltage_UniformBase,std(Velocity_UniformBase,0,2));
%xlabel('Driving Voltage (V)','FontSize',18);
%ylabel('SD of PZT Position (micron)','FontSize',18);

Acquired_Number_of_Loop=p;

Position_ave=mean(Position_UniformBase,2);
first_NOT_NAN_index=find(isnan(Position_ave)==0,1,'first');
Position_used_temp=Position_ave(first_NOT_NAN_index:end);
Voltage_used_temp=Voltage_UniformBase(first_NOT_NAN_index:end);
first_NAN_index=find(isnan(Position_used_temp)>0,1,'first');
Position_used=Position_used_temp(1:(first_NAN_index-1));
Voltage_used=Voltage_used_temp(1:(first_NAN_index-1));


Position_max=max(Position_used);
Position_min=min(Position_used);
Position_range=Position_max-Position_min;
    
Predict_Speed=Position_range/(1/Frequency/2);
    
Time_old=Position_used/Predict_Speed;
    
Time_waveform_1=(Time_old(1)*Time_length_reduction_factor+((Time_old(end)*Time_length_reduction_factor-Time_old(1)*Time_length_reduction_factor)/(N_Waveform/2))):(Time_old(end)*Time_length_reduction_factor-Time_old(1)*Time_length_reduction_factor)/(N_Waveform/2):Time_old(end)*Time_length_reduction_factor;
    
Voltage_waveform_1=interp1(Time_old,Voltage_used,Time_waveform_1);
    
Time_waveform_2=Time_waveform_1-Time_waveform_1(1);
Time_waveform=Time_waveform_2';
Time_waveform((length(Time_waveform_2)+1):(2*length(Time_waveform_2)))=Time_waveform_2(end)+Time_waveform_2';
Voltage_waveform=[Voltage_waveform_1 max(Voltage_waveform_1)-min(Voltage_waveform_1)-Voltage_waveform_1]';

 
    
%% To generate time array

%xlabel('Voltage (V)','FontSize',18);
%ylabel('PZT Position (micron)','FontSize',18);
%hleg=legend('Frequency = 0.001 Hz','Frequency = 0.002 Hz','Frequency = 0.005 Hz','Frequency = 0.01 Hz','Frequency = 0.02 Hz','Frequency = 0.05 Hz','Positon','Best');

plot(Time_waveform,Voltage_waveform);
xlabel('Time (sec)','FontSize',18);
ylabel('Voltage waveform (V)','FontSize',18);
%cd('C:\Users\NTU\Desktop');
%dlmwrite('Time_waveform.txt',Time_waveform,'delimiter','\t','newline','pc', 'precision', '%.6f');
%dlmwrite('Voltage_waveform.txt',Voltage_waveform,'delimiter','\t','newline','pc','precision', '%.6f');
Waveform=[Time_waveform Voltage_waveform];

dlmwrite(sprintf('Waveform_f%g_V%d_O%d.txt',Frequency,Range,Order+1),Waveform,'delimiter','\t','newline','pc','precision', '%.6f');

corr(Position,[1:length(Position)]')

%% to find the linearity (represent by the residual of PZT drift)
Average_Position_Array=[min(Position):(max(Position)-min(Position))/(length(Position)-1):max(Position)]';
Residual=Position-Average_Position_Array;
Average_Time_Array=(1/Sampling_rate):(1/Sampling_rate):(1/Sampling_rate)*length(Average_Position_Array);

plot(Average_Time_Array,Residual);
xlabel('Time (sec)','FontSize',18);
ylabel('Relative PZT Displacement (micron)','FontSize',18);

%Velocity=diff(Residual)./diff(Average_Time_Array');
%Velocity(length(Velocity)+1)=Velocity(length(Velocity));
Error_Fiqure=Velocity(round(length(Velocity)/2));

dlmwrite('Average_Time_Array.txt',Average_Time_Array,'delimiter','\t','newline','pc','precision', '%.6f');
dlmwrite('Residual.txt',Residual,'delimiter','\t','newline','pc','precision', '%.6f');

SD=std(Residual);