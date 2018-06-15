clear all

A=dlmread('D:\MATLAB\141211_PZT waveform file generation\Waveform_f0.01_V10_O1.txt');

finalindex=round(length(A)/2);

Time=A(1:finalindex,1);
Voltage=A(1:finalindex,2);

Time_ms=Time*1000;

plot(Time_ms,Voltage);


groot = fittype( @(a, b, c, x) a*x+b*x.^0.5+c);
%g4 = fittype( @(a, b, c, d, e, x) a*x.^4+b*x.^3+c*x.^2+d*x+e);
g2 = fittype( @(a, b, c, x) a*x.^2+b*x+c);

COEFroot=fit(Time_ms,Voltage,groot);
%COEF4=fit(Time_ms,Voltage,g4);
COEF2=fit(Time_ms,Voltage,g2);

Voltage_fitroot=COEFroot.a*Time_ms+COEFroot.b*Time_ms.^0.5+COEFroot.c;

%Voltage_fit4=COEF4.a*Time_ms.^4+COEF4.b*Time_ms.^3+COEF4.c*Time_ms.^2+COEF4.d*Time_ms+COEF4.e;
Voltage_fit2=COEF2.a*Time_ms.^2+COEF2.b*Time_ms+COEF2.c;

plot(Time_ms,Voltage,Time_ms,Voltage_fit2,Time_ms,Voltage_fitroot);

legend('Original','2nd order fitting', '0.5 fitting');

Residual=Voltage_fitroot-Voltage;

plot(Time_ms,Residual./Voltage*100);
xlabel('Time (ms)');
ylabel('Deviation Ratio (%)');
ylim([-5 5]);

legend('Original','2nd order fitting', '0.5 fitting');

legend('Original','2nd order fitting', '0.5 fitting');