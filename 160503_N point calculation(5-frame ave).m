clear all

cd('D:\Users\TuanShu\130228_FF-OCT data\');
Data=dlmread('filtered.txt'); 

Position_UpperLimit=110;

Position=Data(:,1);

Signal=Data(:,2);

Signal=Signal(Position<Position_UpperLimit);
Position=Position(Position<Position_UpperLimit);

Number_of_Point_per_Carrier_Before_Averaging=140;

N_Array=[3 4 6 9 18 36 144];

Number_of_Points_per_Period=150;
 
Averaging_2=1;
Window_array_2=ones(Averaging_2,1);
NN=1;

for qq=1:length(N_Array)

    N=N_Array(qq);
    %if N==3
    %    Averaging=48;
    %elseif N==4
    %    Averaging=36;
    %elseif N==6
    %    Averaging=24;    
    %elseif N==9
    %    Averaging=16;  
    %elseif N==18
    %    Averaging=8;  
    %elseif N==36
    %    Averaging=4;  
    %elseif N==144
    %    Averaging=1;
    %end
    Averaging=fix(Number_of_Points_per_Period/N);

    Window_array_1=ones(Averaging,1);

    %N=Number_of_Point_per_Carrier_Before_Averaging/Averaging;

    %Window_array=ones(Averaging,1);

    %Signal_ave=conv(Signal,Window_array,'same')/Averaging;

    %C= @ (n,m) factorial(n)/factorial(m)/factorial(n-m);

    clear Signal_New Position_New

    Signal_New(1:fix(length(Signal)/Averaging))=0;
    Position_New(1:fix(length(Signal)/Averaging))=0;
    Signal_New=Signal_New';
    Position_New=Position_New';

    for p=1:fix(length(Signal)/Averaging)

        Signal_New(p)=mean(Signal((1+(p-1)*Averaging):(p*Averaging)));

        Position_New(p)=Position(1+(p-1)*Averaging);
        disp(p);
    end
    
    
    %Signal_New(1:fix(length(Signal)))=0;
    %Position_New(1:fix(length(Signal)))=0;
    %Signal_New=Signal_New';
    %Position_New=Position_New';

    %for p=1:fix(length(Signal))
    %    if p+Averaging < length(Signal)
    %        Signal_New(p)=mean(Signal(p:(p+Averaging)));
    %    else
    %        Signal_New(p)=mean(Signal(p:end));
    %    end        
    %    Position_New(p)=Position(p);
    %    disp(p);
    %end

    %for p=1:fix(length(Signal))
    %    if p+Averaging < length(Signal)
    %        Signal_New(p)=mean(Signal(p:(p+Averaging)));
    %    else
    %        Signal_New(p)=mean(Signal(p:end));
    %    end
    %    Position_New(p)=Position(p);
    %    disp(p);
    %end
    %Signal_New=conv(Signal,Window_array_1,'same')/Averaging;
    
    
    %Better_Array_Length=fix(length(Signal_New)/N)*N;
    %Signal_New=Signal_New(1:Better_Array_Length);
    %Position_New=Position_New(1:Better_Array_Length);
    %Position_New=Position;
    
    clear Signal_New_Env Position_New_Env
    Signal_New_Env(1:fix(length(Signal_New)))=0;
    Position_New_Env(1:fix(length(Position_New)))=0;

    %for p=1:fix(length(Signal_New)/N)
    %    Temp=0;
    %    for m=1:N
    %       for n=1:N
    %           Temp=Temp+(Signal_New((p-1)*N+m)-Signal_New((p-1)*N+n))^2;
    %       end
    %    end
    %    Signal_New_Env(p)=(Temp^0.5)/N;
    %    Position_New_Env(p)=Position_New(1+(p-1)*N);
    %end
    for p=1:(fix((length(Signal_New)-N)/NN))
        Temp=0;
        for m=1:N
           for n=1:N
               Temp=Temp+(Signal_New(1+(p-1)*NN+m)-Signal_New(1+(p-1)*NN+n))^2;
           end
        end
        Signal_New_Env(p)=(Temp^0.5)/N;
        Position_New_Env(p)=Position_New(1+(p-1)*NN);
        disp(p);
    end

    Signal_New_Env=conv(Signal_New_Env,Window_array_2,'same')/Averaging_2;

    if qq==1
        Signal_New_Env_1=Signal_New_Env;
        Position_New_Env_1=Position_New_Env;
    elseif qq==2
        Signal_New_Env_2=Signal_New_Env;
        Position_New_Env_2=Position_New_Env;
    elseif qq==3
        Signal_New_Env_3=Signal_New_Env;
        Position_New_Env_3=Position_New_Env;
    elseif qq==4
        Signal_New_Env_4=Signal_New_Env;
        Position_New_Env_4=Position_New_Env;
    elseif qq==5
        Signal_New_Env_5=Signal_New_Env;
        Position_New_Env_5=Position_New_Env;   
    elseif qq==6
        Signal_New_Env_6=Signal_New_Env;
        Position_New_Env_6=Position_New_Env;
    elseif qq==7
        Signal_New_Env_7=Signal_New_Env;
        Position_New_Env_7=Position_New_Env;
    end
end
plot(Position_New_Env_1(Position_New_Env_1<20),Signal_New_Env_1(Position_New_Env_1<20));
xlabel('Position (micron)');
ylabel('Amplitude (a.u.)');
%plot(Position_New,Signal_New,Position_New_Env_1,(Signal_New_Env_1));
plot(Position_New_Env_1,(Signal_New_Env_1),Position_New_Env_2,(Signal_New_Env_2),Position_New_Env_3,(Signal_New_Env_3),Position_New_Env_4,(Signal_New_Env_4),Position_New_Env_5,(Signal_New_Env_5),Position_New_Env_6,(Signal_New_Env_6),Position_New_Env_7,(Signal_New_Env_7));

plot(Position_New_Env_1,log10(Signal_New_Env_1)*10,Position_New_Env_2,log10(Signal_New_Env_2)*10,Position_New_Env_3,log10(Signal_New_Env_3)*10,Position_New_Env_4,log10(Signal_New_Env_4)*10,Position_New_Env_5,log10(Signal_New_Env_5)*10,Position_New_Env_6,log10(Signal_New_Env_6)*10,Position_New_Env_7,log10(Signal_New_Env_7)*10);

%plot(Position_New_Env_1,log10(Signal_New_Env_1-Signal_New_Env_1(end))*10,Position_New_Env_2,log10(Signal_New_Env_2-Signal_New_Env_2(end))*10,Position_New_Env_3,log10(Signal_New_Env_3-Signal_New_Env_3(end))*10);
plot(Position_New_Env_1,log10(Signal_New_Env_1)*10,Position_New_Env_2,log10(Signal_New_Env_2)*10,Position_New_Env_3,log10(Signal_New_Env_3)*10,Position_New_Env_4,log10(Signal_New_Env_4)*10,Position_New_Env_5,log10(Signal_New_Env_5)*10,Position_New_Env_6,log10(Signal_New_Env_6)*10,Position_New_Env_7,log10(Signal_New_Env_7)*10);
xlabel('Position (micron)');
ylabel('Amplitude (dB)');
legend('3-point','4-point','6-point','9-point','18-point','36-point','144-point');

plot(Position_New_Env_1,Signal_New_Env_1,Position_New_Env_2,Signal_New_Env_2,Position_New_Env_3,Signal_New_Env_3,Position_New_Env_4,Signal_New_Env_4,Position_New_Env_5,Signal_New_Env_5,Position_New_Env_6,Signal_New_Env_6,Position_New_Env_7,Signal_New_Env_7);
xlabel('Position (micron)');
ylabel('Amplitude (a.u.)');
legend('3-point','4-point','6-point','9-point','18-point','36-point','144-point');

dlmwrite('Signal_filtered.txt',Signal_New_Env_1','delimiter','\t','newline','pc');

dlmwrite('Position.txt',Position_New_Env_1','delimiter','\t','newline','pc');