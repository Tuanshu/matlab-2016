clear all
cd('D:\160112_PZT calibration\');
%%
Sampling_rate=100;         %points/sec

If_position_inverted=1;
%%
Wavelength_Laser=0.6328;     %micron

%%

SPF=4000;
LPF=16000;

Start_Pixel=1;
End_Pixel=47.5E5;

Start_Voltage=0;
Range=10;

End_Voltage=Start_Voltage+Range;

Smooth_Window=10;

Voltage_cutoff=0; %how many voltage change near the turning points is neglected.

%%
Data=dlmread('160121_0to10V.txt');        
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
plot(Voltage,Signal);

%% phase wunwrapping
Phase_waveform_original=unwrap(angle(Signal));

%% 去找turning points
Voltage_invert=max(Voltage_for_turningpoint)-Voltage_for_turningpoint;
[minvalue minindex]=findpeaks(Voltage_invert,'MINPEAKDISTANCE',10000,'MINPEAKHEIGHT',End_Voltage-6);
[maxvalue maxindex]=findpeaks(Voltage_for_turningpoint,'MINPEAKDISTANCE',10000,'MINPEAKHEIGHT',End_Voltage-6);
index_all_turningpoints=sort([minindex; maxindex]);

%% to turn the direction of phase
Phase_waveform=Phase_waveform_original;
for p=1:length(index_all_turningpoints)
    Phase_waveform((index_all_turningpoints(p)+1):end)=2*Phase_waveform((index_all_turningpoints(p)))-(Phase_waveform((index_all_turningpoints(p)+1):end));
end
if If_position_inverted==1
    Phase_waveform=-1*Phase_waveform;
end
%%
NNN=5;
Center=index_all_turningpoints(NNN);
 plot(real(Signal));
 xlim([Center-10000 Center+10000]);
%% Position

Position_waveform=Wavelength_Laser/2*Phase_waveform/2/pi;   %由於起點抓的位置不同 有些時候這裡要乘上-1
Time_waveform=((1/Sampling_rate):(1/Sampling_rate):(1/Sampling_rate)*length(Position_waveform))';
plot(Voltage);

plot(Position_waveform);
plot(Voltage,Position_waveform);
xlabel('Voltage (V)');
ylabel('Position (micron)');
plot(Time_waveform,Position_waveform);
xlabel('Time (second)');
ylabel('Position (micron)');

%% Calibration
Order=3;
Velocity=1;                 %micron/sec
Position_Sampling_Resolution=Velocity/Sampling_rate;

%% Forward
Nth_Start_Turning_Point=3;  %1~2 or 2~3
Nth_End_Turning_Point=4;

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

%% Cut specific loop in original data - Roundtrip


Nth_Start_Turning_Point=3;  %1~2 or 2~3
Nth_End_Turning_Point=5;

Position_Patial_Roundtrip=Position_waveform(index_all_turningpoints(Nth_Start_Turning_Point):(index_all_turningpoints(Nth_End_Turning_Point)-1));
Voltage_Patial_Roundtrip=Voltage(index_all_turningpoints(Nth_Start_Turning_Point):(index_all_turningpoints(Nth_End_Turning_Point)-1));
Time_Patial_Roundtrip=Time_waveform(index_all_turningpoints(Nth_Start_Turning_Point):(index_all_turningpoints(Nth_End_Turning_Point)-1));

plot(Time_Patial_Roundtrip,Position_Patial_Roundtrip);
[Position_Patial_Min index_Position_Patial_Min]=min(Position_Patial_Roundtrip);
[Position_Patial_Max index_Position_Patial_Max]=max(Position_Patial_Roundtrip);


plot(Voltage_Patial_Roundtrip,Position_Patial_Roundtrip);

%% Cut specific loop in original data - Forward


Nth_Start_Turning_Point=3;  %1~2 or 2~3
Nth_End_Turning_Point=4;

Position_Patial_Forward=Position_waveform(index_all_turningpoints(Nth_Start_Turning_Point):(index_all_turningpoints(Nth_End_Turning_Point)-1));
Voltage_Patial_Forward=Voltage(index_all_turningpoints(Nth_Start_Turning_Point):(index_all_turningpoints(Nth_End_Turning_Point)-1));
Time_Patial_Forward=Time_waveform(index_all_turningpoints(Nth_Start_Turning_Point):(index_all_turningpoints(Nth_End_Turning_Point)-1));

plot(Time_Patial_Forward,Position_Patial_Forward);
[Position_Patial_Min index_Position_Patial_Min]=min(Position_Patial_Forward);
[Position_Patial_Max index_Position_Patial_Max]=max(Position_Patial_Forward);


plot(Voltage_Patial_Forward,Position_Patial_Forward);

%% Cut specific loop in original data - Backward


Nth_Start_Turning_Point=4;  %1~2 or 2~3
Nth_End_Turning_Point=5;

Position_Patial_Backward=Position_waveform(index_all_turningpoints(Nth_Start_Turning_Point):(index_all_turningpoints(Nth_End_Turning_Point)-1));
Voltage_Patial_Backward=Voltage(index_all_turningpoints(Nth_Start_Turning_Point):(index_all_turningpoints(Nth_End_Turning_Point)-1));
Time_Patial_Backward=Time_waveform(index_all_turningpoints(Nth_Start_Turning_Point):(index_all_turningpoints(Nth_End_Turning_Point)-1));

plot(Time_Patial_Backward,Position_Patial_Backward);
[Position_Patial_Min index_Position_Patial_Min]=min(Position_Patial_Backward);
[Position_Patial_Max index_Position_Patial_Max]=max(Position_Patial_Backward);


plot(Voltage_Patial_Backward,Position_Patial_Backward);
%% x'=a*x+b之a, Ratio array correction
plot(Voltage_Patial_Forward,Position_Patial_Forward,flipud(Voltage_Patial_Backward),flipud(Position_Patial_Backward));

V_f=Voltage_Patial_Forward;
X_f=Position_Patial_Forward;
V_bf=flipud(Voltage_Patial_Backward);
X_bf=flipud(Position_Patial_Backward);

if length(X_f)>length(X_bf)
    X_bf((length(X_bf)+1):length(X_f))=X_bf(length(X_bf));
else
    X_f((length(X_f)+1):length(X_bf))=X_f(length(X_f));
end
plot(V_f,X_f,V_f,X_bf);
X_diff=X_bf-X_f;

% 因為剛好array長度為偶數, 我就先不考慮奇數的情況了
%%
Offset=2245;        %Check: 令Ratio大致左右對稱; 之後驗算時Ratio_Check不要出現0
clear Ratio
for p=1:((length(X_diff))-1-Offset)
    Ratio(p)=X_diff(length(X_diff)-p-Offset)/X_diff(p);
end
Ratio_root=(Ratio).^0.5;
Ratio_root(length(Ratio_root):length(X_diff))=0;        %會有一些0, 但那是因為X_diff前幾項也為0
X_diff_Better=X_diff.*Ratio_root';
plot(1:length(Ratio_root),Ratio_root);
%
plot(X_diff_Better);

%% 驗算
Offset=2250;

for p=1:((length(X_diff_Better))-1-Offset)
    Ratio_Check(p)=X_diff_Better(length(X_diff_Better)-p-Offset)/X_diff_Better(p);
end
plot(1:length(Ratio_Check),Ratio_Check);

%plot(1:length(Ratio),Ratio,1:length(X_diff),X_diff)
%%
%%
plot(Voltage_Patial_Forward,Position_Patial_Forward,flipud(Voltage_Patial_Backward),flipud(Position_Patial_Backward));

V_f=Voltage_Patial_Forward;
X_f=Position_Patial_Forward;
V_bf=flipud(Voltage_Patial_Backward);
X_bf=flipud(Position_Patial_Backward);

if length(X_f)>length(X_bf)
    X_bf((length(X_bf)+1):length(X_f))=X_bf(length(X_bf));
else
    X_f((length(X_f)+1):length(X_bf))=X_f(length(X_f));
end
plot(V_f,X_f,V_f,X_bf);

%V_diff=V_bf-V_f;
X_diff=X_bf-X_f;

X_mean=(X_bf+X_f)/2;

X_mean_ideal=[X_mean(1):(X_mean(end)-X_mean(1))/(length(X_mean)-1):X_mean(end)]';

plot(V_f,X_diff);
plot(V_f,X_mean,V_f,X_mean_ideal);
X_f_better=X_mean_ideal-X_diff/2;
X_bf_better=X_mean_ideal+X_diff/2;


X_b_better=flipud(X_bf_better);

if length(Voltage_Patial_Backward)>length(X_b_better)
    X_b_better((length(X_b_better)+1):length(Voltage_Patial_Backward))=X_b_better(length(X_b_better));
else
    X_b_better=X_b_better(1:length(Voltage_Patial_Backward));
end


plot(V_f,X_f_better,V_f,X_bf_better,V_f,X_f,V_f,X_bf);
legend('X_f_better','X_bf_better','X_f','X_bf');
%% (16/02/01 NEW!) To Decomposing the V-X relation    w, forward, measured based on V_f, X_f_better, X_bf_better (note! not X_b_better)

Delta_vth=1;
Min_V=0;
Max_V=10;


vth_array=Min_V:Delta_vth:Max_V-(Delta_vth);%0:delta_V:max(V_measured)-delta_V;  %[0 0.1 0.2 0.3 0.4 0.5 1 2 4 6 8 9 9.5 9.6 9.7 9.8 9.9];


V_measured_forward=V_f;
X_measured_forward=X_f_better;%Position_Patial_Forward;
V_measured_backward=V_f;
X_measured_backward=X_bf_better;%Position_Patial_Backward;


plot(V_measured_backward,X_measured_backward,V_measured_forward,X_measured_forward,V_measured_forward,X_mean_ideal);
xlabel('V');
ylabel('X');


w_measured_forward=zeros(length(vth_array),1);   
w_measured_backward=zeros(length(vth_array),1);   
%forward
index_vth=1;
for p=1:length(vth_array)
    if p==1
        index_vth_next=find(V_measured_forward>vth_array(p+1),1,'first');
        w_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth));
    elseif p<length(vth_array)
        index_vth_next=find(V_measured_forward>vth_array(p+1),1,'first');
        w_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth))-sum(w_measured_forward(1:p-1));
    else
        index_vth_next=length(V_measured_forward);
        w_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth))-sum(w_measured_forward(1:p-1));
    end

    index_vth=index_vth_next;
end
%backward 16/1/27 NEW!!!!!!!!!!!!!! FORCE第一項之w與forward相同
index_vth=length(X_measured_backward);
for p=1:length(vth_array)
    if p==1
        index_vth_next=find((V_measured_backward)<(Max_V-vth_array(p+1)),1,'last');
        w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth));%w_measured_forward(p);%
    elseif p<length(vth_array)
        index_vth_next=find((V_measured_backward)<(Max_V-vth_array(p+1)),1,'last');
        w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth))-sum(w_measured_backward(1:p-1));
    else
        index_vth_next=1;
        w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth))-sum(w_measured_backward(1:p-1));
    end

    index_vth=index_vth_next;
end

%w_measured_forward(w_measured_forward<0)=0;
%w_measured_backward(w_measured_backward<0)=0;

plot(1:length(w_measured_forward),w_measured_forward,1:length(w_measured_forward),w_measured_backward);
%% Calibrate with 2nd loop of unipolar poling (也就是非initial loading curve, 求得的是w和2*vth)
V_start=0;
V_max_1=10;
V_max_2=10;
V_max_3=10;
V_min_1=0;
V_min_2=0;
V_min_3=0;
d_V=0.01;
V_forward_1=V_start:d_V:V_max_1;
V_backward_1=(V_max_1-d_V):-d_V:V_min_1;
V_forward_2=(V_min_1+d_V):d_V:V_max_2;
V_backward_2=(V_max_2-d_V):-d_V:V_min_2;
V_forward_3=(V_min_2+d_V):d_V:V_max_3;
V_backward_3=(V_max_3-d_V):-d_V:V_min_3;

V=[V_forward_1 V_backward_1];% V_forward_2 V_backward_2 V_forward_3 V_backward_3];
plot(V)

Rate=0.01; %V/sec
Time=Rate:Rate:Rate*length(V);
vth_measured_another_assumption=vth_array/2;
X_another_assumption=zeros(length(V),length(w_measured_forward));
NNNNN=10;
for s=1:(length(w_measured_forward)) %s=NNNNN:NNNNN%
    for p=1:length(V)
        if p==1
            X_another_assumption(p,s)=w_measured_backward(s)*(V(p)+vth_measured_another_assumption(s));
            %X_another_assumption(p,s)=w_measured_forward(s)*(V(p)+vth_measured_another_assumption(s));
        else
            X_another_assumption(p,s)=max(w_measured_forward(s)*(V(p)-vth_measured_another_assumption(s)),min(w_measured_backward(s)*(V(p)+vth_measured_another_assumption(s)),X_another_assumption(p-1,s)));
            %X_another_assumption(p,s)=max(w_measured_forward(s)*(V(p)-vth_measured_another_assumption(s)),min(w_measured_forward(s)*(V(p)+vth_measured_another_assumption(s)),X_another_assumption(p-1,s)));
        end
    end
end
X_total_another_assumption=sum(X_another_assumption,2);

plot(V,X_another_assumption,V,X_total_another_assumption);
xlabel('V (volt)');
ylabel('X (micron)');
legend('X_1','X_2','X_3','X_4','X_5','X_6','X_7','X_8','X_9','X_1_0','X_1_1','X_1_2','X_t_o_t_a_l');
plot(V,X_total_another_assumption,V_f,X_f_better,V_f,X_bf_better);
xlabel('V (volt)');
ylabel('X (micron)');
%%
plot(V,X_total_another_assumption-(X_total_another_assumption(1))+Position_Patial_Roundtrip(1),V_measured_forward,X_measured_forward,V_measured_backward,X_measured_backward);
plot(V,X_total_another_assumption-(X_total_another_assumption(1))+Position_Patial_Roundtrip(1),Voltage_Patial_Roundtrip,Position_Patial_Roundtrip);

        xlabel('V (volt)');
    ylabel('X (micron)');
    legend('V-X (calculated)','V-X (measured)');
    
    plot(V,X_total_another_assumption-(X_total_another_assumption(1))+X_measured_forward(1),V_measured_forward,X_measured_forward,V_measured_backward,X_measured_backward);
%%

plot(V,X_total_another_assumption-(X_total_another_assumption(1))+X_f_better(1),V_f,X_f_better,V_f,X_bf_better);

        xlabel('V (volt)');
    ylabel('X (micron)');
    legend('V-X (calculated)','V-X (measured)');

    %%
    NNN=7;
    plot(X_another_assumption(:,NNN));