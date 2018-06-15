clear all

If_Selective=1;

R1=1;

R2=0.15;

T1=1-R1;
T2=1-R2;

if If_Selective ==0
    T_Sample=T1*T2*T2*T1
    T_Reference=T1*R2*R1*R2*T1
else
    T_Sample=T2*T2
    T_Reference=R2*R1*R2
end

Ratio=T_Sample/T_Reference