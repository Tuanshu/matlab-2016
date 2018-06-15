clear all

Step_Size=1;    %mm

Target_Velocity=1;  %mm/sec

Acc=50;             %Acceleration, 我猜是mm/sec^2
Dec=50;             %Deceleration, 我猜是mm/sec^2

Temporal_Resolution=0.0001;            %sec

Delta=(2*Step_Size/(1/Acc+1/Dec))^0.5-Target_Velocity;

if Delta >= 0        %Case 1, reach Target velocity
    t1=Target_Velocity/Acc;
    t2=(Step_Size/Target_Velocity)-0.5*(1/Acc+1/Dec)*Target_Velocity;
    t3=Target_Velocity/Dec;
elseif Delta < 0
    t1=(Delta+Target_Velocity)/Acc;
    t2=0;
    t3=(Delta+Target_Velocity)/Dec;
end

Total_Time=t1+t2+t3;


Time=0:Temporal_Resolution:(ceil(Total_Time/Temporal_Resolution)*Temporal_Resolution);

Position=zeros([length(Time) 1]);

for p=1:length(Time)
    if Time(p) <= t1
        Position(p)=0.5*Acc*Time(p)^2;
    elseif (Delta >= 0) && (Time(p)<=(t1+t2))
        Position(p)=0.5*Acc*t1^2+Target_Velocity*(Time(p)-t1);
    else Time(p) <= (t1+t2+t3)
        if Delta >= 0
            Position(p)=0.5*Acc*t1^2+Target_Velocity*(Time(p)-t1)-0.5*Dec*(Time(p)-t1-t2)^2;
        elseif Delta < 0
            Position(p)=0.5*Acc*t1^2+(Delta+Target_Velocity)*(Time(p)-t1)-0.5*Dec*(Time(p)-t1-t2)^2;
        end

    end
end

plot(Time,Position);
% Time_Required_To_Reach_Maximum_Velocity=(2*Step_Size/(Acc+(Acc^2)/Dec))^0.5;
% 
% Maximum_Velocity=Acc*Time_Required_To_Reach_Maximum_Velocity;