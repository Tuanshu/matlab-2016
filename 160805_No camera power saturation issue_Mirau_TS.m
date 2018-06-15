clear all
%close all


N=1E10;        %入射之光子數(單位時間內, 在特定波長)
IE=0.5;           %interference efficiency

Ratio=10;     %the ratio between sample signal and sample stray light 1:ratio

R_sample_Total=0.01;   %The TOTAL Signal reflectance from sample

R1_1D=0:0.001:1;        %the reference glass
R2_1D=1:-0.001:0;%0:0.002:1;        %the BS

%R2_1D=R2_1D(end:-1:1);
%
R1=repmat(R1_1D,[length(R2_1D) 1]);
R2=repmat(R2_1D',[1 length(R1_1D)]);
T1=1-R1;
T2=1-R2;

R_stray=R_sample_Total*Ratio/(Ratio+1);   %assumed
R_sample=R_sample_Total/(Ratio+1);

% Transmittance

T_Sample=T1.*T2.*R_sample.*T2.*T1;  %the overall T of Mirau for Sample path
T_Reference=T1.*R2.*R1.*R2.*T1;     %the overall T of Mirau for Reference path
T_StrayLight=T1.*T2.*R_stray.*T2.*T1;  %the overall T of Mirau for Sample path
T_Total=T_Sample+T_Reference+T_StrayLight;


N_Sample=T_Sample*N;
N_Reference=T_Reference*N;
N_StrayLight=T_StrayLight*N;


N_Total=N_Sample+N_Reference+N_StrayLight;

N_Interference=2*IE.*(N_Sample.*N_Reference).^0.5;

N_Noise=(N_Total).^0.5;

SNR=(N_Interference./N_Noise).^2;

Equivalent_Interference_Efficiency=N_Interference./N_Total;  %when N_StrayLight=0 and N_Sample=N_Reference, IE=1, this value is equal to 1


% Try to consider the CCD saturation by NORMALIZING the system
%簡單說就是:
%1. 先找系統穿透率(T_Sample+T_Reference+T_StrayLight
%2. 決定要在Detector收到多少電子 (N_Detector_Set)
%3. 計算因此對應的入射電子數N_Set (注意! 這個N_Set就會是2D array了)
%4. 用這個N_Set來進行SNR_Set計算

N_Detector_Set=1E9;

N_Set=N_Detector_Set./T_Total;


N_Sample_Set=T_Sample.*N_Set;
N_Reference_Set=T_Reference.*N_Set;
N_StrayLight_Set=T_StrayLight.*N_Set;

N_Total_Set=N_Sample_Set+N_Reference_Set+N_StrayLight_Set;




N_Interference_Set=2*IE.*(N_Sample_Set.*N_Reference_Set).^0.5;

N_Noise_Set=(N_Total_Set).^0.5;

SNR_Set=(N_Interference_Set./N_Noise_Set).^2;

imagesc((SNR_Set),'xdata',R1_1D,'ydata',R2_1D);
xlabel('R_1');
ylabel('R_2');
axis equal

colormap(jet);
colorbar
axis equal
set(gca,'YDir','normal');

xlim([0 max(R1_1D)]);
ylim([0 max(R2_1D)]);

x_index_array=1:length(R1_1D);
for p=1:(length(R2_1D))
    [value y_index_array(p)]=max(SNR_Set(:,p));
end
hold on
plot(R1_1D(x_index_array),R2_1D(y_index_array),'w');
hold off

MAX=0;
for p=1:(length(R2_1D)-1)
    if T_Total(x_index_array(p),y_index_array(p))>MAX
        p_record=p;
        MAX=T_Total(x_index_array(p),y_index_array(p));
    end
end
R1_MAX=R1_1D(x_index_array(p_record));
R2_MAX=R2_1D(y_index_array(p_record));

hold on
plot(R1_MAX,R2_MAX,'bo');
hold off

fprintf('R1_Best=%.3f\n',R1_MAX);
fprintf('R2_Best=%.3f\n',R2_MAX);

%% to find max T_Total


%figure;imagesc(N_Reference_Set./N_Sample_Set,'xdata',R1_1D,'ydata',R2_1D);
%xlabel('1st R');
%ylabel('BS R');
%axis equal
%colormap(gray);
%caxis([0 1000])
%figure;imagesc(N_Set,'xdata',R1_1D,'ydata',R2_1D);
%xlabel('1st R');
%ylabel('BS R');
%axis equal
%colormap(gray);
%caxis([0 N_Detector_Set*100])