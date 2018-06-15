clear all

Number_of_Camera=3;

Number_of_Electron_provided=20000*[1:32];

Exposure_Time=1;

%%
Specification_Matrix=zeros(Number_of_Camera,2); %2 specs, FWC and RN

for p=1:3
    Camera=p;
    
    if Camera == 1  %Imperx B0620
    FWC=20000;
    RN=14;
    %DC=1000;
    elseif Camera == 2 %PCO dimax S1
    FWC=36000;
    RN=23;
    %DC=530;

    elseif Camera == 3 %CMOSIS CSI2100
    FWC=2000000;
    RN=945;
    %DC=2600;
    end 
    
    Specification_Matrix(p,1)=FWC;
    Specification_Matrix(p,2)=RN;

end

%% To calculate the number of frame required
SNR=zeros(Number_of_Camera,length(Number_of_Electron_provided));
for p=1:length(Number_of_Electron_provided)
    Number_of_Frame_Required=ones(Number_of_Camera,1);

    Number_of_Frame_Required=ceil(Number_of_Electron_provided(p)./Specification_Matrix(:,1));

    %% Shot noise
    Shot_noise=ones(Number_of_Camera,1);
    Shot_noise=Number_of_Electron_provided(p)^0.5;

    %% Accumulated readout noise

    Accumulated_Readout_Noise=(Number_of_Frame_Required.*Specification_Matrix(:,2).^2).^0.5;
    
  %% Total noise
    Total_noise=(Shot_noise^2+Accumulated_Readout_Noise.^2).^0.5;
    
 %% SNR
    SNR(:,p)=20*log10(Number_of_Electron_provided(p)./Total_noise);
end

plot(Number_of_Electron_provided/1000,SNR');
xlabel('Number of Input Electron (ke-)');
ylabel('SNR before 4-point calc (dB)');
legend('Imperx','PCO','CMOSIS');