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
X_offset=-4.84;

plot(Voltage_Patial_Forward,Position_Patial_Forward,flipud(Voltage_Patial_Backward),flipud(Position_Patial_Backward));

V_f=Voltage_Patial_Forward;
X_f=Position_Patial_Forward+X_offset;
V_bf=flipud(Voltage_Patial_Backward);
X_bf=flipud(Position_Patial_Backward+X_offset);
%% 2/15 下面幾行反而會礙事, 所以註解掉      2/17 又會用到了
if length(X_f)>length(X_bf)
    X_bf((length(X_bf)+1):length(X_f))=X_bf(length(X_bf));
else
    X_f((length(X_f)+1):length(X_bf))=X_f(length(X_f));
end
%plot(V_f,X_f,V_f,X_bf);
%X_diff=X_bf-X_f;

%% 2/15 Change to uniform X basis
X_basis_NumOfSam=1000;
X_basis_min=min([X_f; X_bf]);
X_basis_max=max([X_f; X_bf]);
X_basis=X_basis_min:(X_basis_max-X_basis_min)/(X_basis_NumOfSam-1):X_basis_max;

V_f_Xbasis=interp1(X_f,V_f,X_basis,'linear','extrap');
V_bf_Xbasis=interp1(X_bf,V_bf,X_basis,'linear','extrap');

plot(X_basis,V_f_Xbasis,X_basis,V_bf_Xbasis);

delta_V=V_f_Xbasis-V_bf_Xbasis;
V_com=0.5*(V_f_Xbasis+V_bf_Xbasis);



plot(X_basis,V_com);

%% 2016/02/16 try to find the X', 因為是要求反函數, 會用到類似interp1(V1,X1,X2)的形式(沒打錯!) > 結果好像搞錯了
%V_com_uniform_min=min([V_com]);
%V_com_uniform_max=max([V_com]);
%V_com_uniform=V_com_uniform_min:(V_com_uniform_max-V_com_uniform_min)/(X_basis_NumOfSam-1):V_com_uniform_max;
%X_basis_uniformV=interp1(V_com,X_basis,V_com_uniform);

%plot(V_com,X_basis,V_com_uniform,X_basis_uniformV);


%Coefficient=1;

%X_n=interp1(V_com/max(V_com)*max(X_basis),X_basis,X_basis);
X_n=V_com/max(V_com)*max(X_basis);


plot(X_basis,V_com,X_n,V_com);

plot(X_basis,X_n);      %一開始x的定義好像會對這個圖造成影響, 嘗試改變offset 好, 我終於有這個registration map了


plot(X_n,V_f_Xbasis,X_n,V_bf_Xbasis);

%% 2016/02/16 To get the X(X') (uniform X')

X_n_uniform_min=min([X_n]);
X_n_uniform_max=max([X_n]);
X_n_uniform=X_n_uniform_min:(X_n_uniform_max-X_n_uniform_min)/(X_basis_NumOfSam-1):X_n_uniform_max;
X_uniformXn=interp1(X_n,X_basis,X_n_uniform);

plot(X_n,X_basis,X_n_uniform,X_uniformXn);


%% 2016/02/16 To decomposite  the X(X')

Delta_r=5;
Min_X_n=0;
Max_X_n=400;    %其實最大只有到385


r_array=Min_X_n:Delta_r:Max_X_n-(Delta_r);%0:delta_V:max(V_measured)-delta_V;  %[0 0.1 0.2 0.3 0.4 0.5 1 2 4 6 8 9 9.5 9.6 9.7 9.8 9.9];

S_array=zeros(length(r_array),1);   
%forward
index_r=1;
for p=1:length(r_array)
    if p==1
        index_r_next=find(X_n_uniform>r_array(p+1),1,'first');
        S_array(p)=(X_uniformXn(index_r_next)-X_uniformXn(index_r))/(X_n_uniform(index_r_next)-X_n_uniform(index_r));
    elseif p<length(r_array)
        index_r_next=find(X_n_uniform>r_array(p+1),1,'first');
        S_array(p)=(X_uniformXn(index_r_next)-X_uniformXn(index_r))/(X_n_uniform(index_r_next)-X_n_uniform(index_r))-sum(S_array(1:p-1));
    else
        index_r_next=length(X_n_uniform);
        S_array(p)=(X_uniformXn(index_r_next)-X_uniformXn(index_r))/(X_n_uniform(index_r_next)-X_n_uniform(index_r))-sum(S_array(1:p-1));
    end

    index_r=index_r_next;
end
plot(r_array,S_array);
%% 2016/2/16 Test: Given X_n, calculate X
X_n_Given=X_n_uniform;
X_Calculated=zeros(length(X_n_Given),length(r_array));
for s=1:(length(r_array)) %s=NNNNN:NNNNN%
    if r_array(s)==0
        X_Calculated(:,s)=S_array(s).*X_n_Given;
	elseif r_array(s)>0
        X_Calculated(:,s)=S_array(s).*max(X_n_Given-r_array(s),0);
	elseif r_array(s)<0
        X_Calculated(:,s)=S_array(s).*max(X_n_Given-r_array(s),0);
	end
end

plot(X_n_Given,sum(X_Calculated,2),X_n_Given,X_uniformXn);
%% 好像會用到: To decomposite the X'(X)

Delta_r_rev=5;
Min_X=0;
Max_X=400;    %其實最大只有到385


r_rev_array=Min_X:Delta_r_rev:Max_X-(Delta_r_rev);%0:delta_V:max(V_measured)-delta_V;  %[0 0.1 0.2 0.3 0.4 0.5 1 2 4 6 8 9 9.5 9.6 9.7 9.8 9.9];

S_rev_array=zeros(length(r_rev_array),1);   
%forward
index_r_rev=1;
for p=1:length(r_rev_array)
    if p==1
        index_r_rev_next=find(X_basis>r_rev_array(p+1),1,'first');
        S_rev_array(p)=(X_n(index_r_rev_next)-X_n(index_r_rev))/(X_basis(index_r_rev_next)-X_basis(index_r_rev));
    elseif p<length(r_rev_array)
        index_r_rev_next=find(X_basis>r_rev_array(p+1),1,'first');
        S_rev_array(p)=(X_n(index_r_rev_next)-X_n(index_r_rev))/(X_basis(index_r_rev_next)-X_basis(index_r_rev))-sum(S_rev_array(1:p-1));
    else
        index_r_next=length(X_n_uniform);
        S_rev_array(p)=(X_n(index_r_rev_next)-X_n(index_r_rev))/(X_basis(index_r_rev_next)-X_basis(index_r_rev))-sum(S_rev_array(1:p-1));
    end

    index_r_rev=index_r_rev_next;
end
plot(r_rev_array,S_rev_array);


%% 好像會用到: Given X, calculate X_n (會看起來對的比較不好是因為算出來的起點一定在0附近(無initial offset), 但X_n會有, 因為是基於Vcom算出來的, 但其實只是差offset而已)
X_Given=X_basis;
X_n_Calculated=zeros(length(X_Given),length(r_rev_array));
for s=1:(length(r_rev_array)) %s=NNNNN:NNNNN%
    if r_rev_array(s)==0
        X_n_Calculated(:,s)=S_rev_array(s).*X_Given;
	elseif r_rev_array(s)>0
        X_n_Calculated(:,s)=S_rev_array(s).*max(X_Given-r_rev_array(s),0);
	elseif r_rev_array(s)<0
        X_n_Calculated(:,s)=S_rev_array(s).*max(X_Given-r_rev_array(s),0);
	end
end
%plot(X_n_Given,sum(X_Calculated,2),X_n_Given,X_uniformXn);

plot(X_Given,sum(X_n_Calculated,2),X_Given,X_n);

%% 2016/2/17 用上述function generate X_for_n和X_bak_n                    (待會順便
X_f_n_Calculated=zeros(length(X_f),length(r_rev_array));
X_bf_n_Calculated=zeros(length(X_bf),length(r_rev_array));
for s=1:(length(r_rev_array)) %s=NNNNN:NNNNN%
    if r_rev_array(s)==0
        X_f_n_Calculated(:,s)=S_rev_array(s).*X_f;
        X_bf_n_Calculated(:,s)=S_rev_array(s).*X_bf;
	elseif r_rev_array(s)>0
        X_f_n_Calculated(:,s)=S_rev_array(s).*max(X_f-r_rev_array(s),0);
        X_bf_n_Calculated(:,s)=S_rev_array(s).*max(X_bf-r_rev_array(s),0);
	elseif r_rev_array(s)<0
        X_f_n_Calculated(:,s)=S_rev_array(s).*max(X_f-r_rev_array(s),0);
        X_bf_n_Calculated(:,s)=S_rev_array(s).*max(X_bf-r_rev_array(s),0);
	end
end

plot(V_f,X_f,V_f,X_bf,V_f,sum(X_f_n_Calculated,2),V_f,sum(X_bf_n_Calculated,2));
plot(V_f,sum(X_f_n_Calculated,2),V_f,sum(X_bf_n_Calculated,2),V_f,0.5*(sum(X_f_n_Calculated,2)+sum(X_bf_n_Calculated,2)));  %%奇怪, 平均之後接近原點處有點不直?

%% 2016/2/16 接下來求基於X'(V)之w和vth
Delta_vth=1;
Min_V=0;
Max_V=10;


vth_array=Min_V:Delta_vth:Max_V-(Delta_vth);%0:delta_V:max(V_measured)-delta_V;  %[0 0.1 0.2 0.3 0.4 0.5 1 2 4 6 8 9 9.5 9.6 9.7 9.8 9.9];


V_measured_forward=V_f;
X_measured_forward=X_f_Ratioed_adjusted;%Position_Patial_Forward;
V_measured_backward=V_f;
X_measured_backward=X_bf_Ratioed_adjusted;%Position_Patial_Backward;


%plot(V_measured_backward,X_measured_backward,V_measured_forward,X_measured_forward,V_measured_forward,X_mean_ideal);
%xlabel('V');
%ylabel('X');


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



%%



Center_X=195;

Xc=abs(X_basis-Center_X);

Xc_no_abs=X_basis-Center_X;

index_Xc_went_positive=find(Xc_no_abs>0,1,'first');

Xc_1=Xc(1:index_Xc_went_positive);
Xc_2=Xc((index_Xc_went_positive+1):end);



plot(Xc,delta_V);

%% 因為剛好array長度為偶數, 我就先不考慮奇數的情況了
%
Amp=1;
Offset=2245;        %Check: 令Ratio大致左右對稱; 之後驗算時Ratio_Check不要出現0
clear Ratio
for p=1:((length(X_diff))-1-Offset)
    Ratio(p)=X_diff(length(X_diff)-p-Offset)/X_diff(p);
end
Ratio_root=Amp.*(Ratio).^0.5;
Ratio_root(length(Ratio_root):length(X_diff))=0;        %會有一些0, 但那是因為X_diff前幾項也為0
X_diff_Better=X_diff.*Ratio_root';
plot(1:length(Ratio_root),Ratio_root);
%%
plot(X_diff_Better);

% 驗算
Offset=2250;

for p=1:((length(X_diff_Better))-1-Offset)
    Ratio_Check(p)=X_diff_Better(length(X_diff_Better)-p-Offset)/X_diff_Better(p);
end
plot(1:length(Ratio_Check),Ratio_Check);
%%
%plot(1:length(Ratio),Ratio,1:length(X_diff),X_diff)
% x'=a*x+b之b, offset correction

a=Ratio_root';

X_f_Ratioed=a.*X_f;
X_bf_Ratioed=a.*X_bf;

XX=X_bf_Ratioed-X_f_Ratioed;
plot(XX);
X_mean_ratioed=a.*(X_f+X_bf)/2;
[maxvalue maxindex]=max(X_mean_ratioed);

X_mean_ratioed_ideal=[X_mean_ratioed(1):(maxvalue-X_mean_ratioed(1))/(maxindex-1):X_mean_ratioed(1)+(maxvalue-X_mean_ratioed(1))/(maxindex-1)*(length(X_mean_ratioed)-1)]';
plot(X_mean_ratioed_ideal);
X_mean_ratioed_adjust_to_ideal=X_mean_ratioed_ideal-X_mean_ratioed;
b=X_mean_ratioed_adjust_to_ideal;
X_f_Ratioed_adjusted=X_f_Ratioed+b;
X_bf_Ratioed_adjusted=X_bf_Ratioed+b;


plot(1:length(X_f_Ratioed),X_f_Ratioed_adjusted,1:length(X_bf_Ratioed),X_bf_Ratioed_adjusted);
plot(X_f,X_f_Ratioed_adjusted,X_bf,X_bf_Ratioed_adjusted);
%% (16/02/02 NEW!) To Decomposing the V-X relation    w, forward, measured based on X_f_Ratioed_adjusted and X_bf_Ratioed_adjusted)
Delta_vth=1;
Min_V=0;
Max_V=10;


vth_array=Min_V:Delta_vth:Max_V-(Delta_vth);%0:delta_V:max(V_measured)-delta_V;  %[0 0.1 0.2 0.3 0.4 0.5 1 2 4 6 8 9 9.5 9.6 9.7 9.8 9.9];


V_measured_forward=V_f;
X_measured_forward=X_f_Ratioed_adjusted;%Position_Patial_Forward;
V_measured_backward=V_f;
X_measured_backward=X_bf_Ratioed_adjusted;%Position_Patial_Backward;


%plot(V_measured_backward,X_measured_backward,V_measured_forward,X_measured_forward,V_measured_forward,X_mean_ideal);
%xlabel('V');
%ylabel('X');


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

%%

% plot(Voltage_Patial_Forward,Position_Patial_Forward,flipud(Voltage_Patial_Backward),flipud(Position_Patial_Backward));
% 
% 
% V_f=Voltage_Patial_Forward;
% X_f=Position_Patial_Forward;
% V_bf=flipud(Voltage_Patial_Backward);
% X_bf=flipud(Position_Patial_Backward);
% 
% if length(X_f)>length(X_bf)
%     X_bf((length(X_bf)+1):length(X_f))=X_bf(length(X_bf));
% else
%     X_f((length(X_f)+1):length(X_bf))=X_f(length(X_f));
% end
% plot(V_f,X_f,V_f,X_bf);
% 
% %V_diff=V_bf-V_f;
% X_diff=X_bf-X_f;
% 
% X_mean=(X_bf+X_f)/2;
% 
% X_mean_ideal=[X_mean(1):(X_mean(end)-X_mean(1))/(length(X_mean)-1):X_mean(end)]';
% 
% plot(V_f,X_diff);
% plot(V_f,X_mean,V_f,X_mean_ideal);
% X_f_better=X_mean_ideal-X_diff/2;
% X_bf_better=X_mean_ideal+X_diff/2;
% 
% 
% X_b_better=flipud(X_bf_better);
% 
% if length(Voltage_Patial_Backward)>length(X_b_better)
%     X_b_better((length(X_b_better)+1):length(Voltage_Patial_Backward))=X_b_better(length(X_b_better));
% else
%     X_b_better=X_b_better(1:length(Voltage_Patial_Backward));
% end
% 
% 
% plot(V_f,X_f_better,V_f,X_bf_better,V_f,X_f,V_f,X_bf);
% legend('X_f_better','X_bf_better','X_f','X_bf');
%% (16/02/01 NEW!) To Decomposing the V-X relation    w, forward, measured based on V_f, X_f_better, X_bf_better (note! not X_b_better)
% 
% Delta_vth=1;
% Min_V=0;
% Max_V=10;
% 
% 
% vth_array=Min_V:Delta_vth:Max_V-(Delta_vth);%0:delta_V:max(V_measured)-delta_V;  %[0 0.1 0.2 0.3 0.4 0.5 1 2 4 6 8 9 9.5 9.6 9.7 9.8 9.9];
% 
% 
% V_measured_forward=V_f;
% X_measured_forward=X_f_better;%Position_Patial_Forward;
% V_measured_backward=V_f;
% X_measured_backward=X_bf_better;%Position_Patial_Backward;
% 
% 
% plot(V_measured_backward,X_measured_backward,V_measured_forward,X_measured_forward,V_measured_forward,X_mean_ideal);
% xlabel('V');
% ylabel('X');
% 
% 
% w_measured_forward=zeros(length(vth_array),1);   
% w_measured_backward=zeros(length(vth_array),1);   
% %forward
% index_vth=1;
% for p=1:length(vth_array)
%     if p==1
%         index_vth_next=find(V_measured_forward>vth_array(p+1),1,'first');
%         w_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth));
%     elseif p<length(vth_array)
%         index_vth_next=find(V_measured_forward>vth_array(p+1),1,'first');
%         w_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth))-sum(w_measured_forward(1:p-1));
%     else
%         index_vth_next=length(V_measured_forward);
%         w_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth))-sum(w_measured_forward(1:p-1));
%     end
% 
%     index_vth=index_vth_next;
% end
% %backward 16/1/27 NEW!!!!!!!!!!!!!! FORCE第一項之w與forward相同
% index_vth=length(X_measured_backward);
% for p=1:length(vth_array)
%     if p==1
%         index_vth_next=find((V_measured_backward)<(Max_V-vth_array(p+1)),1,'last');
%         w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth));%w_measured_forward(p);%
%     elseif p<length(vth_array)
%         index_vth_next=find((V_measured_backward)<(Max_V-vth_array(p+1)),1,'last');
%         w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth))-sum(w_measured_backward(1:p-1));
%     else
%         index_vth_next=1;
%         w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth))-sum(w_measured_backward(1:p-1));
%     end
% 
%     index_vth=index_vth_next;
% end
% 
% %w_measured_forward(w_measured_forward<0)=0;
% %w_measured_backward(w_measured_backward<0)=0;
% 
% plot(1:length(w_measured_forward),w_measured_forward,1:length(w_measured_forward),w_measured_backward);
%% Calibrate with 2nd loop of unipolar poling (也就是非initial loading curve, 求得的是w和2*vth)
V_start=0.01;
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
plot(V,X_total_another_assumption,V_f,X_f_Ratioed_adjusted,V_f,X_bf_Ratioed_adjusted);
xlabel('V (volt)');
ylabel('X (micron)');
%%
plot(V,X_total_another_assumption-(X_total_another_assumption(1))+Position_Patial_Roundtrip(1),V_measured_forward,X_measured_forward,V_measured_backward,X_measured_backward);
plot(V,X_total_another_assumption-(X_total_another_assumption(1))+Position_Patial_Roundtrip(1),Voltage_Patial_Roundtrip,Position_Patial_Roundtrip);

        xlabel('V (volt)');
    ylabel('X (micron)');
    legend('V-X (calculated)','V-X (measured)');
    
    plot(V,X_total_another_assumption-(X_total_another_assumption(1))+X_measured_forward(1),V_measured_forward,X_measured_forward,V_measured_backward,X_measured_backward);
%% Trial: try to transfrom the X_total_another_assumption to the original one

%1st step: to break down the X_total_another_assumption
X_total_another_assumption_f=X_total_another_assumption(1:length(V_forward_1));
X_total_another_assumption_b=X_total_another_assumption((length(V_forward_1)+1):end);
X_total_another_assumption_bf=flipud(X_total_another_assumption_b);
V_f_Trial=V_forward_1';
V_bf_Trial=flipud(V_backward_1');
plot(X_total_another_assumption_bf);
[maxvalue maxindex]=max(V_f);
V_f_for_interpolation=[V_f(1):(maxvalue-V_f(1))/(maxindex-1):V_f(1)+(maxvalue-V_f(1))/(maxindex-1)*(length(V_f)-1)]';
a_Trial=interp1(V_f_for_interpolation,a,V_f_Trial);
b_Trial=interp1(V_f_for_interpolation,b,V_f_Trial);
plot(V_f_for_interpolation,a,V_f_Trial,a_Trial);
X_f_Derived_Original=(X_total_another_assumption_f-(X_total_another_assumption_f(1))+X_f(1)-b_Trial)./a_Trial;
X_bf_Derived_Original=(X_total_another_assumption_bf-(X_total_another_assumption_bf(1))+X_bf(1)-b_Trial)./a_Trial;

%plot(V_f_Trial,X_total_another_assumption_f-X_total_another_assumption_f(1)+X_f_Ratioed_adjusted(1),V_f_Trial,X_total_another_assumption_bf-X_total_another_assumption_bf(1)+X_bf_Ratioed_adjusted(1),V_f,X_f_Ratioed_adjusted,V_f,X_bf_Ratioed_adjusted);%,V_f,X_f,V_f,X_bf);
plot(V_f_Trial,X_total_another_assumption_f,V_f_Trial,X_total_another_assumption_bf,V_f,X_f_Ratioed_adjusted,V_f,X_bf_Ratioed_adjusted);%,V_f,X_f,V_f,X_bf);
X_Test_f=(X_f_Ratioed_adjusted-b)./a;
X_Test_bf=(X_bf_Ratioed_adjusted-b)./a;
X_Trial_f=(X_total_another_assumption_f-X_total_another_assumption_f(1)+X_f_Ratioed_adjusted(1)-b_Trial)./a_Trial;
X_Trial_bf=(X_total_another_assumption_bf-X_total_another_assumption_bf(1)+X_bf_Ratioed_adjusted(1)-b_Trial)./a_Trial;
plot(V_f,X_Test_f,V_f,X_Test_bf);

plot(V_f_Trial,X_Trial_f,V_f_Trial,X_Trial_bf,V_f,X_Test_f,V_f,X_Test_bf);
% 所以結論是: OK, 但必須注意一開始的offset, 因為是unipolar
%% 然後要breakdown a and b array, to match the form of vth etc.
a_array=zeros([length(vth_array) 1]);
b_array=zeros([length(vth_array) 1]);
index_vth=1;
for p=1:length(vth_array)
    if p<length(vth_array)
        index_vth_next=find(V_f>vth_array(p+1),1,'first');
        a_array(p)=mean(a(index_vth:index_vth_next));
        b_array(p)=min(b(index_vth:index_vth_next));
    else
        index_vth_next=length(V_f);
        a_array(p)=mean(a(index_vth:index_vth_next));
        b_array(p)=min(b(index_vth:index_vth_next));
        
    end
    index_vth=index_vth_next;
end
%% Note: 由於x'=ax+b, 這個transform是apply在x上, 故會有一張x versus x'的作圖, 在該圖上, 可以用類似找w和vth的方式, 找出S function的參數
%% (2/5 突然發現, 我一開始的a和b的定義可能就不太好) > 這代表根本不是transform嗎?
plot(X_f_Derived_Original-X_f_Derived_Original(1),X_total_another_assumption_f-X_total_another_assumption_f(1),X_bf_Derived_Original-X_bf_Derived_Original(1),X_total_another_assumption_bf-X_total_another_assumption_bf(1));
for s=1:(length(a_array)) %s=NNNNN:NNNNN%   %每一個
    for p=1:length(V)
        if p==1
            X_another_assumption_back_to_Original(p,s)=w_measured_backward(s)*(V(p)+vth_measured_another_assumption(s));
            %X_another_assumption(p,s)=w_measured_forward(s)*(V(p)+vth_measured_another_assumption(s));
        else
            X_another_assumption_back_to_Original(p,s)=max(w_measured_forward(s)*(V(p)-vth_measured_another_assumption(s)),min(w_measured_backward(s)*(V(p)+vth_measured_another_assumption(s)),X_another_assumption(p-1,s)));
            %X_another_assumption(p,s)=max(w_measured_forward(s)*(V(p)-vth_measured_another_assumption(s)),min(w_measured_forward(s)*(V(p)+vth_measured_another_assumption(s)),X_another_assumption(p-1,s)));
        end
    end
end





        xlabel('V (volt)');
    ylabel('X (micron)');
    legend('V-X (calculated)','V-X (measured)');

    %%
    NNN=2;
    plot(X_another_assumption(:,NNN));