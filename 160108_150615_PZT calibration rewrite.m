clear all
cd('D:\AMO\150526_PZT');
%%
Sampling_rate=1000;         %points/sec


%%
Wavelength_Laser=0.6328;     %micron

%%

SPF=1000;
LPF=10000;

Start_Pixel=4000;
End_Pixel=2150000;

Start_Voltage=0;
Range=10;

End_Voltage=Start_Voltage+Range;

Smooth_Window=10;

Voltage_cutoff=0; %how many voltage change near the turning points is neglected.

%%
Data=dlmread('140307_3rd_measurement (1micronsec)_SR1000.txt');        
plot(Data(:,2));
Voltage_read=Data(Start_Pixel:End_Pixel,2);
Signal_read=Data(Start_Pixel:End_Pixel,1);
plot(Voltage_read);
%%  filtering the signal
FFT_Signal_read=fft(Signal_read);
plot(real(FFT_Signal_read));
FFT_Signal_read(1:SPF)=0;
FFT_Signal_read(LPF:end)=0;
Signal_read_new=ifft(FFT_Signal_read);

%% to delete elements from array
Voltage=smooth(Voltage_read,Smooth_Window);
Signal=smooth(Signal_read_new,Smooth_Window);
Voltage_for_turningpoint=smooth(Voltage_read,Smooth_Window*40); %*10 for 4 micron per sec, *20 for 2 micron per second, *40 for 1 micron per second

%Signal(Voltage>(End_Voltage-Voltage_cutoff))=[];    %注意順序
%Voltage(Voltage>(End_Voltage-Voltage_cutoff))=[];

%% phase wunwrapping
Phase_waveform_original=unwrap(angle(Signal));

%% 去找turning points
Voltage_invert=max(Voltage_for_turningpoint)-Voltage_for_turningpoint;
[minvalue minindex]=findpeaks(Voltage_invert,'MINPEAKDISTANCE',10000,'MINPEAKHEIGHT',End_Voltage-1);
[maxvalue maxindex]=findpeaks(Voltage_for_turningpoint,'MINPEAKDISTANCE',10000,'MINPEAKHEIGHT',End_Voltage-1);
index_all_turningpoints=sort([minindex; maxindex]);

%% to turn the direction of phase
Phase_waveform=Phase_waveform_original;
for p=1:length(index_all_turningpoints)
    Phase_waveform((index_all_turningpoints(p)+1):end)=2*Phase_waveform((index_all_turningpoints(p)))-(Phase_waveform((index_all_turningpoints(p)+1):end));
end

%% Position

Position_waveform=Wavelength_Laser/2*Phase_waveform/2/pi;   %由於起點抓的位置不同 有些時候這裡要乘上-1
Time_waveform=((1/Sampling_rate):(1/Sampling_rate):(1/Sampling_rate)*length(Position_waveform))';

plot(Position_waveform);
plot(Voltage,Position_waveform);
xlabel('Voltage (V)');
ylabel('Position (micron)');
plot(Time_waveform,Position_waveform);
xlabel('Time (second)');
ylabel('Position (micron)');


Voltage_Patial=Voltage(Time_waveform<131);
Position_waveform_Patial=Position_waveform(Time_waveform<131);
Time_waveform_Patial=Time_waveform(Time_waveform<131);
Voltage_Patial=Voltage_Patial(Time_waveform_Patial>10.1);
Position_waveform_Patial=Position_waveform_Patial(Time_waveform_Patial>10.1);

plot(Voltage_Patial,Position_waveform_Patial);
xlabel('Voltage (V)');
ylabel('Position (micron)');
xlim([-2 10]);




%% Calibration
Order=3;
Velocity=1;                 %micron/sec
Position_Sampling_Resolution=Velocity/Sampling_rate;

%% Forward
Nth_Start_Turning_Point=2;  %1~2 or 2~3
Nth_End_Turning_Point=3;

Position_Patial=Position_waveform(index_all_turningpoints(Nth_Start_Turning_Point):(index_all_turningpoints(Nth_End_Turning_Point)-1));
Voltage_Patial=Voltage(index_all_turningpoints(Nth_Start_Turning_Point):(index_all_turningpoints(Nth_End_Turning_Point)-1));
Time_Patial=Time_waveform(index_all_turningpoints(Nth_Start_Turning_Point):(index_all_turningpoints(Nth_End_Turning_Point)-1));

Position_Patial=Position_Patial(Voltage_Patial<10);
Time_Patial=Time_Patial(Voltage_Patial<10);
Voltage_Patial=Voltage_Patial(Voltage_Patial<10);

plot(Time_Patial,Position_Patial);
[Position_Patial_Min index_Position_Patial_Min]=min(Position_Patial);
[Position_Patial_Max index_Position_Patial_Max]=max(Position_Patial);

Position_Patial_Ideal=interp1([Time_Patial(index_Position_Patial_Min) Time_Patial(index_Position_Patial_Max)],[Position_Patial(index_Position_Patial_Min) Position_Patial(index_Position_Patial_Max)],Time_Patial);

plot(Time_Patial,Position_Patial_Ideal,Time_Patial,Position_Patial);

plot(Voltage_Patial,Position_Patial);
Position_UniformBase=(min(Position_Patial):Position_Sampling_Resolution:max(Position_Patial))';
Time_UniformBase=((1/Sampling_rate):(1/Sampling_rate):(1/Sampling_rate)*length(Position_UniformBase))';
Voltage_UniformBase=interp1(Position_Patial,Voltage_Patial,Position_UniformBase);

%% (NEW) take only the positive voltage

Voltage_min=0;
Voltage_max=10;

index_start=find(Voltage_Patial>Voltage_min,1,'first');
Voltage_Patial_New=Voltage_Patial(index_start:end);
Position_Patial_New=Position_Patial(index_start:end);
Position_Patial_New=Position_Patial_New-min(Position_Patial_New);
plot(Voltage_Patial_New,Position_Patial_New);
%% (16/01/08 NEW!) To Decomposing the V-X relation (Initial loading curve)
V_measured=Voltage_Patial_New;
X_measured=Position_Patial_New;
plot(Voltage_Patial_New,Position_Patial_New);

delta_V=1;

Vth_measured=delta_V:delta_V:max(V_measured);
Vth_measured(length(Vth_measured)+1)=10;

w_measured=zeros(length(Vth_measured),1);

for p=1:length(Vth_measured)
    w_measured(p)=








%% (NEW!) to generate the new V-X relation
Expansion_Ratio=460/330*(64/68.5);

Position_Patial_New_Used=(Position_Patial_New)*Expansion_Ratio;
Voltage_Patial_New_Used=Voltage_Patial_New;
plot(Voltage_Patial_New_Used,Position_Patial_New_Used);


%% (NEW) To generate a waveform of specified velocity

Wished_Velocity=0.813; %micron/sec
Max_position=max(Position_Patial_New_Used);
Max_time_ms=Max_position/Wished_Velocity*1000;

Time_ms=0:Max_time_ms;

Position_Uniform_Based_Used=Time_ms*Wished_Velocity/1000;
Voltage_Uniform_Based_Used=interp1(Position_Patial_New_Used,Voltage_Patial_New_Used,Position_Uniform_Based_Used);

plot(Position_Patial_New_Used,Voltage_Patial_New_Used,Position_Uniform_Based_Used,Voltage_Uniform_Based_Used);

plot(Time_ms,Voltage_Uniform_Based_Used);

%% (NEW) Fitting

Time=Time_ms';
Voltage=Voltage_Uniform_Based_Used';

plot(Time,Voltage);
gs = fittype( @(a, b, c, x) a*x+b*x.^0.5+c);
g4 = fittype( @(a, b, c, d, e, x) a*x.^4+b*x.^3+c*x.^2+d*x+e);
g3 = fittype( @(a, b, c, d, x) a*x.^3+b*x.^2+c*x+d);
g2 = fittype( @(a, b, c, x) a*x.^2+b*x+c);
g1 = fittype( @(a, b, x) a*x+b);

%COEF3=fit(Time,Voltage,g3);
%COEF1=fit(Time,Voltage,g1);
%COEF2=fit(Time,Voltage,g2,'StartPoint',[0 COEF1.a COEF1.b]);
%COEF4=fit(Time,Voltage,g4);%,'StartPoint',[0 0 COEF2.a COEF2.b COEF2.c]);
COEFS=fit(Time,Voltage,gs,'StartPoint',[0 0.0161 -2.6919]);%,'StartPoint',[0 0 COEF2.a COEF2.b COEF2.c]);

fprintf('%s %s %s',COEFS.a,COEF.b, COEF.c);

%Voltage_fit1=COEF1.a*Time+COEF1.b;
%Voltage_fit2=COEF2.a*Time.^2+COEF2.b*Time+COEF2.c;
%Voltage_fit3=COEF3.a*Time.^3+COEF3.b*Time.^2+COEF3.c*Time+COEF3.d;
%Voltage_fit4=COEF4.a*(Time).^4+COEF4.b*(Time).^3+COEF4.c*(Time).^2+COEF4.d*(Time)+COEF4.e;
Voltage_fits=COEFS.a*Time+COEFS.b*Time.^0.5+COEFS.c;

%plot(Time,Voltage,Time,Voltage_fit3);
%plot(Time,Voltage,Time,Voltage_fit2,Time,Voltage_fit4);
plot(Time,Voltage,Time,Voltage_fits);

legend('Original','fitted');

%% Error Calculation (consider only the forward)

Error=max(abs(Position_Patial_Ideal-Position_Patial));

plot(Time_Patial,Position_Patial_Ideal,Time_Patial,Position_Patial);
xlabel('Time (second)');
ylabel('Position (micron)');

%% Backward

Nth_Start_Turning_Point_2=3;    %2~3 or 3~4
Nth_End_Turning_Point_2=4;

Position_Patial_2=Position_waveform(index_all_turningpoints(Nth_Start_Turning_Point_2):(index_all_turningpoints(Nth_End_Turning_Point_2)-1));
Voltage_Patial_2=Voltage(index_all_turningpoints(Nth_Start_Turning_Point_2):(index_all_turningpoints(Nth_End_Turning_Point_2)-1));
Position_Patial_2=Position_Patial_2(Voltage_Patial_2<10);
Voltage_Patial_2=Voltage_Patial_2(Voltage_Patial_2<10);

plot(Voltage_Patial,Position_Patial);
Position_UniformBase_2=(max(Position_Patial_2):(-1*Position_Sampling_Resolution):min(Position_Patial_2))';
Time_UniformBase_2=((1/Sampling_rate):(1/Sampling_rate):(1/Sampling_rate)*length(Position_UniformBase_2))';
Voltage_UniformBase_2=interp1(Position_Patial_2,Voltage_Patial_2,Position_UniformBase_2);



%% to generate the 3 loops + backwalking
Number_of_Cycle=3;
Pre_Voltage=0;
%%
Voltage_half_cycle_forward=Voltage_UniformBase;
Voltage_half_cycle_backward=Voltage_UniformBase_2;
Voltage_one_cycle=[Voltage_half_cycle_forward; Voltage_half_cycle_backward];
Voltage=repmat(Voltage_one_cycle,[Number_of_Cycle 1]);

%% backwalking
Voltage_pre_cycle_backward=Voltage_half_cycle_backward(Voltage_half_cycle_backward<Pre_Voltage);

%%
Voltage_with_pre=[Voltage_pre_cycle_backward; Voltage];
Time_with_pre=(1/Sampling_rate)*(1:length(Voltage_with_pre))';

plot(Time_with_pre,Voltage_with_pre);



Output=[Time_with_pre Voltage_with_pre];

dlmwrite(sprintf('Waveform_SR%dHz_V%dmicronsec_Order%d_withbackwalking.txt',Sampling_rate,Velocity,Order),Output,'delimiter','\t','newline','pc','precision', '%.6f');

