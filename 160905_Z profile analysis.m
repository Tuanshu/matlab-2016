Z_Porfile=dlmread('C:\Users\Owner\Desktop\160901_Images\160905_Imperx Z profile for comparison.txt');
Z_Porfile=dlmread('F:\Processed Data\AVG_Reslice of AVG_160905_Imperx_648x3_2070fps_4x_Ave_Factor_32.txt');
Z_Porfile=dlmread('F:\Processed Data\AVG_Reslice of AVG_160905_Imperx_648x3_2070fps_4x_Upper.txt');
Z_Porfile=dlmread('F:\Processed Data\160905_Imperx_648x3_2070fps_4x_forearm_2_Ave_Factor_32.txt');

%X_Porfile=dlmread('C:\Users\Owner\Desktop\160901_Images\SFR_X profile.txt');

% plot(X_Porfile(:,1)*0.5,X_Porfile(:,2)/16,'-o');
% 
%  xlabel('Lateral Position (micron)');
%  ylabel('Signal (ADU)');

Depth=Z_Porfile(:,1)*0.2;
LogSig=log10(Z_Porfile(:,2));
plot(Depth,LogSig);

 xlabel('Depth (micron)');
 ylabel('log(Signal (DN))');
xlim([0 100]);
ylim([0 2]);

FitResult=fit(Depth(abs(Depth-70)<10),LogSig(abs(Depth-70)<10),'poly1');

Fit_LogSig=FitResult.p1*Depth+FitResult.p2;

plot(Depth,LogSig,Depth,Fit_LogSig);

 xlabel('Depth (micron)');
 ylabel('log(Signal (DN))');
xlim([0 100]);
ylim([0 20]);


plot(Depth,20*(LogSig-Fit_LogSig(Depth==30)),Depth,20*(Fit_LogSig-Fit_LogSig(Depth==30)));
xlim([0 140]);

ylim([-8 0]);
 xlabel('Depth (micron)');
 ylabel('Signal (dB)');