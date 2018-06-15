clear all;

%%

Sampling_rate=100;
d_t=1/Sampling_rate;

Predicted_X_V_Ratio=430/10; %micron/volt
Predicted_PZT_Velocity=0.8;   %micron/sec
Calculated_V_Rate=Predicted_PZT_Velocity/Predicted_X_V_Ratio;   %volt/sec
d_V=d_t*Calculated_V_Rate;


V_start=0;
V_max_1=10;
V_max_2=10;
V_max_3=10;
V_max_4=8;
V_max_5=6;

V_min_1=0;
V_min_2=0;
V_min_3=0;
V_min_4=0;
V_min_5=0;

V_forward_1=V_start:d_V:V_max_1;
V_backward_1=(V_max_1-d_V):-d_V:V_min_1;
V_forward_2=(V_min_1+d_V):d_V:V_max_2;
V_backward_2=(V_max_2-d_V):-d_V:V_min_2;
V_forward_3=(V_min_2+d_V):d_V:V_max_3;
V_backward_3=(V_max_3-d_V):-d_V:V_min_3;
V_forward_4=(V_min_3+d_V):d_V:V_max_4;
V_backward_4=(V_max_4-d_V):-d_V:V_min_4;
V_forward_5=(V_min_4+d_V):d_V:V_max_5;
V_backward_5=(V_max_5-d_V):-d_V:V_min_5;


V=[V_forward_1 V_backward_1 V_forward_2 V_backward_2 V_forward_3 V_backward_3  V_forward_4 V_backward_4  V_forward_5 V_backward_5]';
T=[d_t:d_t:d_t*length(V)]';

plot(T,V);

Output=[T V];
cd('D:\160112_PZT calibration');
dlmwrite(sprintf('Waveform_SR%dHz_V%dmicronsec_10V.txt',Sampling_rate,Predicted_PZT_Velocity),Output,'delimiter','\t','newline','pc','precision', '%.6f');