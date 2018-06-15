clear all

D1=49;
D1_2=51;

F1=50;
F2=9;

%M0=[1 D1;0 1];
M0=[1 50+D1;0 1];
M0_2=[1 D1_2;0 1];
M1=[1 0;-1/F1 1];
Min=[1 140;0 1];
M2=[1 0;-1/F2 1];
%M3=[1 D2;0 1];

Vin=[0; 1];


Vout=M2*Min*M1*M0*Vin;
Vout_2=M2*Min*M1*M0_2*Vin;

D2=Vout(1)/Vout(2)*-1;
D2_2=Vout_2(1)/Vout_2(2)*-1;