clear all

Target_Speed=0.1;   %mm/sec

Max_Moving_Speed=9;    %mm/sec

Avaliable_Time_per_mm=1/Target_Speed;   %sec
Time_for_Moving_per_mm=1/Max_Moving_Speed;  %sec

Avaliable_Time_for_Settle_per_mm=Avaliable_Time_per_mm-Time_for_Moving_per_mm;  %sec

Step_Size=3/1000;%0.78/2/1.4/4/1000;    %mm

Number_of_Step_per_mm=1/Step_Size;  %#

Max_Settle_Time=Avaliable_Time_for_Settle_per_mm/Number_of_Step_per_mm; %sec

Max_Settle_Time_ms=Max_Settle_Time*1000;