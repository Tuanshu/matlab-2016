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
%% (16/01/08 NEW!) To Decomposing the V-X relation (Initial loading curve)
 index_negative_w=[];

V_measured=Voltage_Patial_New;
X_measured=Position_Patial_New;
plot(V_measured,X_measured);
xlabel('V');
ylabel('X');
delta_V=0.5;

vth_measured=[0 0.1 0.2 0.3 0.4 0.5 1 1.5 2 3 4 6 8 9 9.5 9.7 9.8 9.9];%0:delta_V:max(V_measured)-delta_V;

w_measured=zeros(length(vth_measured),1);
index_vth=1;
for p=1:length(vth_measured)
    if p==1
        index_vth_next=find(V_measured>vth_measured(p+1),1,'first');
        w_measured(p)=(X_measured(index_vth_next)-X_measured(index_vth))/(V_measured(index_vth_next)-V_measured(index_vth));
    elseif p<length(vth_measured)
        index_vth_next=find(V_measured>vth_measured(p+1),1,'first');
        w_measured(p)=(X_measured(index_vth_next)-X_measured(index_vth))/(V_measured(index_vth_next)-V_measured(index_vth))-sum(w_measured(1:p-1));
    else
        index_vth_next=length(V_measured);
        w_measured(p)=(X_measured(index_vth_next)-X_measured(index_vth))/(V_measured(index_vth_next)-V_measured(index_vth))-sum(w_measured(1:p-1));
    end

    index_vth=index_vth_next;
end
index_negative_w=find(w_measured<0,1,'first');
if isempty(index_negative_w)==0
    vth_measured=vth_measured(1:(index_negative_w-1));
    w_measured=w_measured(1:(index_negative_w-1));
end

%% (16/01/08 NEW!) to graph out this result



V_start=0;
V_max_1=10;
V_max_2=10;
V_max_3=10;
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

plot(V,X_total);%,V,X_total);
xlabel('V (volt)');
ylabel('X (micron)');
legend('X_1','X_2','X_3','X_4','X_5','X_6','X_7','X_8','X_9','X_1_0','X_1_1','X_1_2','X_t_o_t_a_l');

%
plot(V,X_total,V_measured,X_measured);
plot(V,X_total,Voltage,Position_waveform);

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
  %  Position_Segments=cell(length(index_all_turningpoints));    %2 for x and y coordinate

  %  for p=1:length(index_all_turningpoints)
  %      if p==1
  %          Position_Segments{1}=Position_waveform(1:index_all_turningpoints(1));
  %      else
  %          Position_Segments{p}=Position_waveform(index_all_turningpoints(p-1):index_all_turningpoints(p));
  %      end
  %  end

    
    %%
 %   Max_array_Size=max([index_all_turningpoints(1) diff(index_all_turningpoints)']);
 %   for p=1:length(index_all_turningpoints)
 %       if length(Position_Segments{QQQ})<Max_array_Size
 %           Position_Segments{QQQ}(length(Position_Segments{QQQ}):Max_array_Size)=0;
 %       end
 %   end
            %%
       % QQQ=1;
    %plot(Position_Segments{QQQ});
    %%
index_Voltage_Zero=find(Voltage>0,1,'first');
plot(V,X_total,Voltage,Position_waveform-Position_waveform(index_Voltage_Zero));
    xlabel('V');
    ylabel('X');
    
plot(V,X_total);
    xlabel('V');
    ylabel('X');
%

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



%% Calibrate with 2nd loop of unipolar poling (也就是非initial loading curve, 求得的是w和2*vth)
    
V_start=0;
V_max_1=10;
V_max_2=10;
V_max_3=10;
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

V=[V_forward_1 V_backward_1];% V_forward_2 V_backward_2 V_forward_3 V_backward_3];
plot(V)

Rate=0.01; %V/sec
Time=Rate:Rate:Rate*length(V);

vth_measured_another_assumption=vth_measured/2;

X_another_assumption=zeros(length(V),length(w_measured));
for s=1:length(w_measured)
    for p=1:length(V)
        if p==1
            X_another_assumption(p,s)=max(w_measured(s)*(V(p)-vth_measured_another_assumption(s)),w_measured(s)*(V(p)+vth_measured_another_assumption(s)));
        else
            X_another_assumption(p,s)=max(w_measured(s)*(V(p)-vth_measured_another_assumption(s)),min(w_measured(s)*(V(p)+vth_measured_another_assumption(s)),X_another_assumption(p-1,s)));
        end
    end
end
X_total_another_assumption=sum(X_another_assumption,2);

plot(V,X_another_assumption,V,X_total_another_assumption);
xlabel('V (volt)');
ylabel('X (micron)');
legend('X_1','X_2','X_3','X_4','X_5','X_6','X_7','X_8','X_9','X_1_0','X_1_1','X_1_2','X_t_o_t_a_l');

%
plot(V,X_total_another_assumption,V_measured,X_measured);
plot(V,X_total_another_assumption-(X_total_another_assumption(1))+Position_Patial_Roundtrip(1),Voltage_Patial_Roundtrip,Position_Patial_Roundtrip);

        xlabel('V (volt)');
    ylabel('X (micron)');
    legend('V-X (calculated)','V-X (measured)');
    
    
    plot(Voltage,Position_waveform);
    
    xlabel('V (volt)');
    ylabel('X (micron)');
    
    
    plot(Time,V);
    xlabel('Time (second)');
    ylabel('V (volt)');
    
    %plot(Time,X_total);
    %xlabel('Time (second)');
    %ylabel('Displacement (micron)');
    
    %plot(V,X_total);
    %    xlabel('V (volt)');
    %ylabel('X (micron)');
    
%Output=[Time_with_pre Voltage_with_pre];

%dlmwrite(sprintf('Waveform_SR%dHz_V%dmicronsec_Order%d_withbackwalking.txt',Sampling_rate,Velocity,Order),Output,'delimiter','\t','newline','pc','precision', '%.6f');
%% Try to find the memory-free nonlinearity
%Interpole to the same basis
Ideal_dVoltage_Forward=(max(Voltage_Patial_Forward)-min(Voltage_Patial_Forward))/(length(Voltage_Patial_Forward)-1);
Ideal_dVoltage_Backward=(min(Voltage_Patial_Backward)-max(Voltage_Patial_Backward))/(length(Voltage_Patial_Backward)-1);

Ideal_Voltage_Forward=min(Voltage_Patial_Forward):Ideal_dVoltage_Forward:max(Voltage_Patial_Forward);
Ideal_Voltage_Backward=max(Voltage_Patial_Backward):Ideal_dVoltage_Backward:min(Voltage_Patial_Backward);
Position_Patial_Backward_Forwardbasis=interp1(Ideal_Voltage_Backward(length(Voltage_Patial_Backward):-1:1),Position_Patial_Backward(length(Voltage_Patial_Backward):-1:1),Ideal_Voltage_Forward,'linear','extrap')';

plot(Voltage_Patial_Forward,Position_Patial_Forward,Voltage_Patial_Forward,Position_Patial_Backward_Forwardbasis);

Position_Patial_Mean=(Position_Patial_Forward+Position_Patial_Backward_Forwardbasis)/2;

%% (16/01/22 NEW!) To Decomposing the V-X relation (Iterative method), first to get the w arrays:   w, forward, measured

Delta_vth=1;
Min_V=0;
Max_V=10;


vth_array=Min_V:Delta_vth:Max_V-(Delta_vth);%0:delta_V:max(V_measured)-delta_V;


V_measured_forward=Voltage_Patial_Forward;
X_measured_forward=Position_Patial_Forward;
V_measured_backward=Voltage_Patial_Backward;
X_measured_backward=Position_Patial_Backward;


plot(V_measured_backward,X_measured_backward,V_measured_forward,X_measured_forward);
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
%backward
index_vth=1;
for p=1:length(vth_array)
    if p==1
        index_vth_next=find((Max_V-V_measured_backward)>vth_array(p+1),1,'first');
        w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth));
    elseif p<length(vth_measured_backward)
        index_vth_next=find((Max_V-V_measured_backward)>vth_array(p+1),1,'first');
        w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth))-sum(w_measured_backward(1:p-1));
    else
        index_vth_next=length(V_measured_backward);
        w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth))-sum(w_measured_backward(1:p-1));
    end

    index_vth=index_vth_next;
end
%% (16/01/22 NEW!) NOT To Decomposing the V-X relation, just get the slope array (X_d array)

Delta_vth=1;
Min_V=0;
Max_V=10;


vth_array=Min_V:Delta_vth:Max_V-(Delta_vth);%0:delta_V:max(V_measured)-delta_V;


V_measured_forward=Voltage_Patial_Forward;
X_measured_forward=Position_Patial_Forward;
V_measured_backward=Voltage_Patial_Backward;
X_measured_backward=Position_Patial_Backward;


plot(V_measured_backward,X_measured_backward,V_measured_forward,X_measured_forward);
xlabel('V');
ylabel('X');


sdxdv_measured_forward=zeros(length(vth_array),1);      %dxdv==dx/dv
sdxdv_measured_backward=zeros(length(vth_array),1);   
%forward
index_vth=1;
for p=1:length(vth_array)
    if p==1
        index_vth_next=find(V_measured_forward>vth_array(p+1),1,'first');
        sdxdv_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth));
    elseif p<length(vth_array)
        index_vth_next=find(V_measured_forward>vth_array(p+1),1,'first');
        sdxdv_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth));
    else
        index_vth_next=length(V_measured_forward);
        sdxdv_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth));
    end
    index_vth=index_vth_next;
end
%backward
index_vth=1;
for p=1:length(vth_array)
    if p==1
        index_vth_next=find((Max_V-V_measured_backward)>vth_array(p+1),1,'first');
        sdxdv_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth));
    elseif p<length(vth_array)
        index_vth_next=find((Max_V-V_measured_backward)>vth_array(p+1),1,'first');
        sdxdv_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth));
    else
        index_vth_next=length(V_measured_backward);
        sdxdv_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth));
    end
    index_vth=index_vth_next;
end

plot(vth_array,sdxdv_measured_forward,vth_array,(sdxdv_measured_backward));

%% Note, 以上的sdxdv_measured_backward是倒著排序的, 所以要轉回來
sdxdv_for_m=sdxdv_measured_forward;
sdxdv_bak_m=flipud(sdxdv_measured_backward);

%% Iterative method (~GNA)
%initial condition
S=ones([length(vth_array) 1]);
dxdv=sdxdv_measured_forward;

d_S=1E-10;      
d_dxdv=1E-10;
      

sdxdv_for = @ (S,dxdv) S.*dxdv;
sdxdv_bak = @ (S,dxdv) S.*flipud(dxdv); %因我們只求forward的dx/dv, 而這裡的應該要放backward的dx/dv
                   
det_array = @ (Q1,Q2,Q3,Q4) Q1.*Q4-Q2.*Q3;

%% Start the fitting    
MSE_total_OK=0.05;
Max_Number_of_Loop=1000;
Merit_Best=999999999999999999999;
Current_Loop=1;
MSE_A=99999999;
MSE_B=99999999;
MSE_total=99999999;
while (Current_Loop<Max_Number_of_Loop) && (MSE_total>MSE_total_OK)
    d_sdxdv_for_d_S = (sdxdv_for(S+d_S,dxdv)-sdxdv_for(S,dxdv))/d_S;
    d_sdxdv_bak_d_S = (sdxdv_bak(S+d_S,dxdv)-sdxdv_bak(S,dxdv))/d_S;
                       
    d_sdxdv_for_d_dxdv = (sdxdv_for(S,dxdv+d_dxdv)-sdxdv_for(S,dxdv))/d_S;
    d_sdxdv_bak_d_dxdv = (sdxdv_bak(S,dxdv+d_dxdv)-sdxdv_bak(S,dxdv))/d_S;
           
    delta_sdxdv_for=sdxdv_for_m-sdxdv_for(S,dxdv);
    delta_sdxdv_bak=sdxdv_bak_m-sdxdv_bak(S,dxdv);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        increment_S = det_array(delta_sdxdv_for,d_sdxdv_for_d_dxdv,delta_sdxdv_bak,d_sdxdv_bak_d_dxdv)./det_array(d_sdxdv_for_d_S,d_sdxdv_for_d_dxdv,d_sdxdv_bak_d_S,d_sdxdv_bak_d_dxdv);
                        increment_dxdv = det_array(d_sdxdv_for_d_S,delta_sdxdv_for,d_sdxdv_bak_d_S,delta_sdxdv_bak)./det_array(d_sdxdv_for_d_S,d_sdxdv_for_d_dxdv,d_sdxdv_bak_d_S,d_sdxdv_bak_d_dxdv);
                        %dthickness = det_array(Q1,Q2,delta_A,Q4,Q5,delta_B,Q7,Q8,delta_C)./det_array(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9);
                        %嘗試用MATLAB的左除法
                        %increment=[d_sdxdv_for_d_S d_sdxdv_for_d_dxdv;d_sdxdv_bak_d_S d_sdxdv_bak_d_dxdv].\[delta_sdxdv_for;delta_sdxdv_bak];
                        
                        S=S+0.5*increment_S;
                        dxdv=dxdv+0.1*increment_dxdv;
                        
                        S(isnan(S))=1;
                        S(isinf(S))=1;
                        dxdv(isnan(dxdv))=100;
                        dxdv(isinf(dxdv))=1;

                        
    MSE_sdxdv_for=(sum((delta_sdxdv_for./sdxdv_for_m).^2)/length(delta_sdxdv_for)).^0.5;
    MSE_sdxdv_bak=(sum((delta_sdxdv_bak./sdxdv_bak_m).^2)/length(delta_sdxdv_bak)).^0.5;
    MSE_total=(sum((delta_sdxdv_for./sdxdv_for_m).^2)/length(delta_sdxdv_for)+sum((delta_sdxdv_bak./sdxdv_bak_m).^2)/length(delta_sdxdv_bak)).^0.5;
                        
                        Current_Loop=Current_Loop+1;
                        fprintf('%d/%d     \n',Current_Loop,Max_Number_of_Loop);
                        %fprintf('%d/%d\n',MSE_A,MSE_A_OK);
                        %fprintf('%d/%d\n',MSE_B,MSE_B_OK);
                        %fprintf('%d/%d\n',MSE_C,MSE_C_OK);
                        %fprintf('%d/%d\n',MSE_total,MSE_total_OK);
                        

                        
end
subplot(2,1,1)
plot(S.*dxdv);
subplot(2,1,2)
plot(S.*flipud(dxdv));