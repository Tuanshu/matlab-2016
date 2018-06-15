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
    elseif p<length(vth_array)
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
sx_for_m=X_measured_forward(1:min(length(X_measured_forward),length(X_measured_backward)));
sx_bak_m=flipud(X_measured_backward(1:min(length(X_measured_forward),length(X_measured_backward))));
% simplified case
%sdxdv_for_m(:)=sdxdv_measured_forward.*sdxdv_measured_backward;
%sdxdv_bak_m(:)=sdxdv_measured_forward.*(sdxdv_measured_backward).^-2;

%% Iterative method (~LMA)
%initial condition
S=ones([length(min(length(X_measured_forward),length(X_measured_backward))) 1]);
%S(1:5)=1.5;
%S=rand(length(vth_array),1)*0.1+ones([length(vth_array) 1])-0.05;
x=X_measured_forward(1:min(length(X_measured_forward),length(X_measured_backward)));%./(rand(length(vth_array),1)*0.1+ones([length(vth_array) 1])-0.05);

d_S=1E-10;      
d_x=1E-10;
      

sx_for = @ (S,x) S.*x;
sx_bak = @ (S,x) S.*flipud(x); %因我們只求forward的dx/dv, 而這裡的應該要放backward的dx/dv
                   
% simplified case
%sdxdv_for = @ (S,dxdv) dxdv.*S;
%sdxdv_bak = @ (S,dxdv) dxdv.*(S).^-2; %因我們只求forward的dx/dv, 而這裡的應該要放backward的dx/dv
%sdxdv_for = @ (S,dxdv) dxdv+S;
%sdxdv_bak = @ (S,dxdv) dxdv+2*S; %因我們只求forward的dx/dv, 而這裡的應該要放backward的dx/dv


det_array = @ (Q1,Q2,Q3,Q4) Q1.*Q4-Q2.*Q3;

%% Start the fitting    

v=1.5;
lambda=0.01;

MSE_total_OK=0.00001;
Max_Number_of_Loop=1000;
MSE_best=999999999999999999999;
Current_Loop=1;
MSE_A=99999999;
MSE_B=99999999;
MSE_total=99999999;
clear MSE_record MSE_sx_for_record MSE_sx_bak_record
while (Current_Loop<Max_Number_of_Loop) && (MSE_total>MSE_total_OK)
    S_pre=S;
    x_pre=x;
    d_sx_for_d_S = (sx_for(S+d_S,x)-sx_for(S,x))/d_S;
    d_sx_bak_d_S = (sx_bak(S+d_S,x)-sx_bak(S,x))/d_S;
                       
    d_sx_for_d_x = (sx_for(S,x+d_x)-sx_for(S,x))/d_x;
    d_sx_bak_d_x = (sx_bak(S,x+d_x)-sx_bak(S,x))/d_x;
           
    delta_sx_for=sx_for_m-sx_for(S,x);
    delta_sx_bak=sx_bak_m-sx_bak(S,x);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GO=1;
    w=1;
    while (GO>0)
        if GO==1    %for lambda
            Q1=d_sx_for_d_S+lambda/v*d_sx_for_d_S;
            Q2=d_sx_for_d_x;
            Q3=d_sx_bak_d_S;
            Q4=d_sx_bak_d_x+lambda/v*d_sx_bak_d_x;
            increment_S = det_array(delta_sx_for,Q2,delta_sx_bak,Q4)./det_array(Q1,Q2,Q3,Q4);
            increment_x = det_array(Q1,delta_sx_for,Q3,delta_sx_bak)./det_array(Q1,Q2,Q3,Q4);
                        %dthickness = det_array(Q1,Q2,delta_A,Q4,Q5,delta_B,Q7,Q8,delta_C)./det_array(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9);
                        %嘗試用MATLAB的左除法
                        %increment=[d_sdxdv_for_d_S d_sdxdv_for_d_dxdv;d_sdxdv_bak_d_S d_sdxdv_bak_d_dxdv].\[delta_sdxdv_for;delta_sdxdv_bak];
            S_temp=S+0.1*increment_S;
            x_temp=x+0.1*increment_x;
                        
            S_temp(isnan(S_temp))=S_pre(isnan(S_temp));
            S_temp(isinf(S_temp))=S_pre(isinf(S_temp));
            x_temp(isnan(x_temp))=x_pre(isnan(x_temp));
            x_temp(isinf(x_temp))=x_pre(isinf(x_temp));
            
            
            delta_sx_for_temp=sx_for_m-sx_for(S_temp,x_temp);
            delta_sx_bak_temp=sx_bak_m-sx_bak(S_temp,x_temp);
            
            MSE_total_temp(1)=(sum((delta_sx_for_temp./sx_for_m).^2)/length(delta_sx_for_temp)+sum((delta_sx_bak_temp./sx_bak_m).^2)/length(delta_sx_bak_temp)).^0.5;
            
            if MSE_total_temp(1)<=MSE_best
                S=S_temp;
                x=x_temp;
                lambda=lambda/v;
                MSE_best=MSE_total_temp(1);
                MSE_record(Current_Loop)=MSE_total_temp(1);
                GO_record(Current_Loop)=GO;
                GO=-9999999;
            end
        elseif GO==2    %for lambda/v
            Q1=d_sx_for_d_S+lambda*d_sx_for_d_S;
            Q2=d_sx_for_d_x;
            Q3=d_sx_bak_d_S;
            Q4=d_sx_bak_d_x+lambda*d_sx_bak_d_x;
            increment_S = det_array(delta_sx_for,Q2,delta_sx_bak,Q4)./det_array(Q1,Q2,Q3,Q4);
            increment_x = det_array(Q1,delta_sx_for,Q3,delta_sx_bak)./det_array(Q1,Q2,Q3,Q4);
                        %dthickness = det_array(Q1,Q2,delta_A,Q4,Q5,delta_B,Q7,Q8,delta_C)./det_array(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9);
                        %嘗試用MATLAB的左除法
                        %increment=[d_sdxdv_for_d_S d_sdxdv_for_d_dxdv;d_sdxdv_bak_d_S d_sdxdv_bak_d_dxdv].\[delta_sdxdv_for;delta_sdxdv_bak];
            S_temp=S+0.1*increment_S;
            x_temp=x_temp+0.1*increment_x;
                        
            S_temp(isnan(S_temp))=S_pre(isnan(S_temp));
            S_temp(isinf(S_temp))=S_pre(isinf(S_temp));
            x_temp(isnan(x_temp))=x_pre(isnan(x_temp));
            x_temp(isinf(x_temp))=x_pre(isinf(x_temp));

            delta_sx_for_temp=sx_for_m-sx_for(S_temp,x_temp);
            delta_sx_bak_temp=sx_bak_m-sx_bak(S_temp,x_temp);
            
            MSE_total_temp(2)=(sum((delta_sx_for_temp./sx_for_m).^2)/length(delta_sx_for_temp)+sum((delta_sx_bak_temp./sx_bak_m).^2)/length(delta_sx_bak_temp)).^0.5;
                
            if MSE_total_temp(2)<=MSE_best
                S=S_temp;
                x=x_temp;
                MSE_best=MSE_total_temp(2);
                MSE_record(Current_Loop)=MSE_total_temp(2);
                GO_record(Current_Loop)=GO;
                GO=-9999999;
            end
            
        else
            Q1=d_sx_for_d_S+lambda*(v^w)*d_sx_for_d_S;
            Q2=d_sx_for_d_x;
            Q3=d_sx_bak_d_S;
            Q4=d_sx_bak_d_x+lambda*(v^w)*d_sx_bak_d_x;
            increment_S = det_array(delta_sx_for,Q2,delta_sx_bak,Q4)./det_array(Q1,Q2,Q3,Q4);
            increment_x = det_array(Q1,delta_sx_for,Q3,delta_sx_bak)./det_array(Q1,Q2,Q3,Q4);
                        %dthickness = det_array(Q1,Q2,delta_A,Q4,Q5,delta_B,Q7,Q8,delta_C)./det_array(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9);
                        %嘗試用MATLAB的左除法
                        %increment=[d_sdxdv_for_d_S d_sdxdv_for_d_dxdv;d_sdxdv_bak_d_S d_sdxdv_bak_d_dxdv].\[delta_sdxdv_for;delta_sdxdv_bak];
            S_temp=S+0.1*increment_S;
            x_temp=x_temp+0.1*increment_x;
                        
            S_temp(isnan(S_temp))=S_pre(isnan(S_temp));
            S_temp(isinf(S_temp))=S_pre(isinf(S_temp));
            x_temp(isnan(x_temp))=x_pre(isnan(x_temp));
            x_temp(isinf(x_temp))=x_pre(isinf(x_temp));

            delta_sx_for_temp=sx_for_m-sx_for(S_temp,x_temp);
            delta_sx_bak_temp=sx_bak_m-sx_bak(S_temp,x_temp);
            
            MSE_total_temp(GO)=(sum((delta_sx_for_temp./sx_for_m).^2)/length(delta_sx_for_temp)+sum((delta_sx_bak_temp./sx_bak_m).^2)/length(delta_sx_bak_temp)).^0.5;
                

            if lambda*(v^w)==inf
                S=S_temp;
                x=x_temp;
                lambda=99999999999999999;
                MSE_best=MSE_total_temp(GO);
                MSE_record(Current_Loop)=MSE_total_temp(GO);
                GO_record(Current_Loop)=GO;
                GO=-9999999;      
            elseif MSE_total_temp(GO)<=MSE_best
                S=S_temp;
                x=x_temp;
                lambda=lambda*(v^w);
                MSE_best=MSE_total_temp(GO);
                MSE_record(Current_Loop)=MSE_total_temp(GO);
                GO_record(Current_Loop)=GO;
                GO=-9999999;  
            else
                w=w+1;
            end
            
            
        end
        GO=GO+1;
    end
           
    Current_Loop=Current_Loop+1;
    fprintf('%d/%d     \n',Current_Loop,Max_Number_of_Loop);
                        %fprintf('%d/%d\n',MSE_A,MSE_A_OK);
                        %fprintf('%d/%d\n',MSE_B,MSE_B_OK);
                        %fprintf('%d/%d\n',MSE_C,MSE_C_OK);
                        %fprintf('%d/%d\n',MSE_total,MSE_total_OK);
                        

    %MSE_sdxdv_for=(sum((delta_sdxdv_for./sdxdv_for_m).^2)/length(delta_sdxdv_for)).^0.5;
    %MSE_sdxdv_bak=(sum((delta_sdxdv_bak./sdxdv_bak_m).^2)/length(delta_sdxdv_bak)).^0.5;           
end
subplot(2,1,1)
plot(sx_for(S,x));
subplot(2,1,2)
plot(sx_bak(S,x));

plot(MSE_record(10:end));
%subplot(2,1,1)
%plot(dxdv+S);
%subplot(2,1,2)
%plot(dxdv+2*S);