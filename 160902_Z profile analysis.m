Z_Porfile=dlmread('C:\Users\Owner\Desktop\160901_Images\AVG_Reslice of AVG_160830_PCO_4160fps_4x_200microsec_200x200_forearm_Ave_Factor_64_depth 200 micron_z_profile.txt');
X_Porfile=dlmread('C:\Users\Owner\Desktop\160901_Images\SFR_X profile.txt');

plot(X_Porfile(:,1)*0.5,X_Porfile(:,2)/16,'-o');

 xlabel('Lateral Position (micron)');
 ylabel('Signal (ADU)');

Depth=Z_Porfile(:,1)*0.2;
LogSig=log10(Z_Porfile(:,2));
plot(Depth,LogSig);

 xlabel('Depth (micron)');
 ylabel('log(Signal (DN))');
xlim([0 100]);
ylim([0 20]);

FitResult=fit(Depth(abs(Depth-70)<30),LogSig(abs(Depth-70)<30),'poly1');

Fit_LogSig=FitResult.p1*Depth+FitResult.p2;

plot(Depth,LogSig,Depth,Fit_LogSig);

 xlabel('Depth (micron)');
 ylabel('log(Signal (DN))');
xlim([0 100]);
ylim([0 20]);


plot(Depth,20*(LogSig-Fit_LogSig(Depth==20)),Depth,20*(Fit_LogSig-Fit_LogSig(Depth==20)));
xlim([0 200]);

ylim([-20 0]);
 xlabel('Depth (micron)');
 ylabel('Signal (dB)');