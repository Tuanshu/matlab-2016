clear all

Folder_Path='I:\160908_Photon Transfer Curves\';

Parameter='Power';

Max_ADU=4096;
QE=0.5;
FWC_PCO_Spec=36000;
FWC_Imperx_Spec=20000;

Readout_Noise_PCO=23;
Readout_Noise_Imperx=16;


Output_Imperx=dlmread([Folder_Path sprintf('Output_Different %s_Imperx.txt',Parameter)]);

Output_PCO=dlmread([Folder_Path sprintf('Output_Different %s_PCO.txt',Parameter)]);

plot(Output_Imperx(:,1),Output_Imperx(:,2),'-o',Output_PCO(:,1),Output_PCO(:,2),'-o');
legend('Imperx B0620','PCO HS4','Location','NorthWest');
xlim([0 3500]);
ylim([0 20]);
xlabel('Signal (ADU)');
ylabel('Noise (ADU)');
grid on


g = fittype( @(a,b,c,x) (a*x.^2+b*x+c)); %a: excess, b: shot, c: dark
Fit_Imperx = fit(Output_Imperx(:,1), Output_Imperx(:,2).^2, g, 'StartPoint', [0.015, 0.0569, -0]);
Fit_Imperx.a
Fit_Imperx.b
Fit_Imperx.c

Fit_PCO = fit(Output_PCO(:,1), Output_PCO(:,2).^2, g, 'StartPoint', [0.015, 0.0569, -0]);
Fit_PCO.a
Fit_PCO.b
Fit_PCO.c

Ratio=Fit_Imperx.b/Fit_PCO.b;
Ratio_Expected=Readout_Noise_Imperx/Readout_Noise_PCO;
%FWC_estimated_Imperx=Max_ADU/Fit_Imperx.b;
%FWC_estimated_PCO=Max_ADU/Fit_PCO.b;

plot(Output_Imperx(:,1),Output_Imperx(:,2).^2,'-o',Output_PCO(:,1),Output_PCO(:,2).^2,'-o');
legend('Imperx B0620','PCO HS4','Location','NorthWest');
xlim([0 3500]);
ylim([0 400]);
xlabel('Signal (ADU)');
ylabel('Noise^2 (ADU)');
grid on


Ratio_Noise=Fit_Imperx.c/Fit_PCO.c;
Ratio_Noise_Expected=FWC_PCO_Spec/FWC_Imperx_Spec;


plot(FWC_PCO_Spec/1000,Fit_PCO.b,'o',FWC_Imperx_Spec/1000,Fit_Imperx.b,'o');
xlim([0 60]);
ylim([0 0.3]);
xlabel('Full Well Capacity (ke)');
ylabel('2nd Order Noise Coefficient');
grid on

Fit_FWC=fit([FWC_PCO_Spec/1000 FWC_Imperx_Spec/1000]',[Fit_PCO.b Fit_Imperx.b]','poly1');

FWC_Array=0:1:60;

Coef_Line_Est=Fit_FWC.p1*FWC_Array+Fit_FWC.p2;


plot(FWC_PCO_Spec/1000,Fit_PCO.b,'o',FWC_Imperx_Spec/1000,Fit_Imperx.b,'o',FWC_Array,Coef_Line_Est);
xlim([0 60]);
ylim([0 0.3]);
xlabel('Full Well Capacity (ke)');
ylabel('2nd Order Noise Coefficient');
grid on