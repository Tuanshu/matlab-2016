clear all

%B1

N_B1=648^2;

QE_B1=0.5;

Frame_Rate_B1=1000;    %at ROI 648*76

Light_Power_B1=14;     

FWC_B1=20000;

CE_B1=1;       %normalized

Ave_Factor=64;

Corresponding_Thickness=0.78/1.4/2;    %micron

AAA=N_B1*Frame_Rate_B1*FWC_B1/(QE_B1*CE_B1*Light_Power_B1);

%%

if_set_FOV=1;
if if_set_FOV == 1
    Sampling_Resolution=0.45;

    FOV_Length=0.5;   %mm
    N_B2=(FOV_Length/Sampling_Resolution*1000)^2;
else
N_B2=648^2;
end
Case=3;     %1 for best, 2 for worse, 3 for normal

if Case == 1

    QE_B2=0.8;

    Light_Power_B2=40;   

    FWC_B2=20000;

    CE_B2=3;       %normalized
    
elseif Case == 2
    QE_B2=0.3;

    Light_Power_B2=10;   

    FWC_B2=20000;

    CE_B2=1.5;       %normalized
elseif Case == 3
    QE_B2=0.5;

    Light_Power_B2=20;   

    FWC_B2=20000;

    CE_B2=2;       %normalized
    
end
Max_Frame_Rate_B2=AAA/N_B2*QE_B2*Light_Power_B2/FWC_B2*CE_B2;
Max_Scanning_Speed_B2=Max_Frame_Rate_B2/Ave_Factor*Corresponding_Thickness;

