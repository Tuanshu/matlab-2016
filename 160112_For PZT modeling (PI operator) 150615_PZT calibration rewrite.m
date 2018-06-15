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
 index_negative_w=13;

V_measured=Voltage_Patial_New;
X_measured=Position_Patial_New;
plot(V_measured,X_measured);
xlabel('V');
ylabel('X');
delta_V=0.5;

vth_measured=0:delta_V:max(V_measured);

w_measured=zeros(length(vth_measured),1);
index_vth=1;
for p=1:length(vth_measured)
    if p==1
        index_vth_next=find(V_measured>vth_measured(p+1),1,'first');
        w_measured(p)=(X_measured(index_vth_next)-X_measured(index_vth))/delta_V;
    elseif p<length(vth_measured)
        index_vth_next=find(V_measured>vth_measured(p+1),1,'first');
        w_measured(p)=(X_measured(index_vth_next)-X_measured(index_vth))/delta_V-sum(w_measured(1:p-1));
    else
        index_vth_next=length(V_measured);
        w_measured(p)=(X_measured(index_vth_next)-X_measured(index_vth))/delta_V-sum(w_measured(1:p-1));
    end

    index_vth=index_vth_next;
end
%index_negative_w=find(w_measured<0,1,'first');
if isempty(index_negative_w)==0
    vth_measured=vth_measured(1:(index_negative_w-1));
    w_measured=w_measured(1:(index_negative_w-1));
end

%% (16/01/08 NEW!) to graph out this result



V_start=0;
V_max_1=10;
V_max_2=9;
V_max_3=8;
V_min_1=0;
V_min_2=0;
V_min_3=0;
d_V=0.1;
V_forward_1=V_start:d_V:V_max_1;
V_backward_1=(V_max_1-d_V):-d_V:V_min_1;
V_forward_2=(V_min_1+d_V):d_V:V_max_2;
V_backward_2=(V_max_2-d_V):-d_V:V_min_2;
V_forward_3=(V_min_2+d_V):d_V:V_max_3;
V_backward_3=(V_max_3-d_V):-d_V:V_min_3;

V=[V_forward_1 V_backward_1 V_forward_2 V_backward_2 V_forward_3 V_backward_3];
plot(V)

Rate=0.01; %V/sec
Time=Rate:Rate:Rate*length(V);


X=zeros(length(V),length(w_measured));
for s=1:length(w_measured)
    for p=1:length(V)
        if p==1
            X(p,s)=max(w_measured(s)*(V(p)-vth_measured(s)),min(w_measured(s)*(V(p)+vth_measured(s)),0));
        else
            X(p,s)=max(w_measured(s)*(V(p)-vth_measured(s)),min(w_measured(s)*(V(p)+vth_measured(s)),X(p-1,s)));
        end
    end
end
X_total=sum(X,2);

plot(V,X,V,X_total);
        xlabel('V (volt)');
    ylabel('X (micron)');
legend('X_1','X_2','X_3','X_4','X_5','X_6','X_7','X_8','X_9','X_1_0','X_1_1','X_1_2','X_t_o_t_a_l');

%
plot(V,X_total,V_measured,X_measured);
%plot(V,X_total,Voltage,Position_waveform);

        xlabel('V (volt)');
    ylabel('X (micron)');
    legend('V-X (calculated)','V-X (measured)');
    
    
    plot(Voltage,Position_waveform);
    
    
    
    plot(Time,V);
    xlabel('Time (second)');
    ylabel('V (volt)');
    
    plot(Time,X_total);
    xlabel('Time (second)');
    ylabel('Displacement (micron)');
    
    plot(V,X_total);
        xlabel('V (volt)');
    ylabel('X (micron)');
    %% Try to find the memory-free nonlinearity
    Position_Segments=cell(length(index_all_turningpoints));    %2 for x and y coordinate

    for p=1:length(index_all_turningpoints)
        if p==1
            Position_Segments{1}=Position_waveform(1:index_all_turningpoints(1));
        else
            Position_Segments{p}=Position_waveform(index_all_turningpoints(p-1):index_all_turningpoints(p));
        end
    end

    
    %%
    Max_array_Size=max([index_all_turningpoints(1) diff(index_all_turningpoints)']);
    for p=1:length(index_all_turningpoints)
        if length(Position_Segments{QQQ})<Max_array_Size
            Position_Segments{QQQ}(length(Position_Segments{QQQ}):Max_array_Size)=0;
        end
    end
            %%
        QQQ=1;
    plot(Position_Segments{QQQ});
    %%
index_Voltage_Zero=find(Voltage>0,1,'first');
plot(V,X_total,Voltage,Position_waveform-Position_waveform(index_Voltage_Zero));
    xlabel('V');
    ylabel('X');
    
plot(V,X_total);
    xlabel('V');
    ylabel('X');
%
    
    %%
%Output=[Time_with_pre Voltage_with_pre];

%dlmwrite(sprintf('Waveform_SR%dHz_V%dmicronsec_Order%d_withbackwalking.txt',Sampling_rate,Velocity,Order),Output,'delimiter','\t','newline','pc','precision', '%.6f');

