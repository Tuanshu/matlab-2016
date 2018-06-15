clear all


Xmax_micron=10;        %micron

deltaX_micron=0.01;    %micron
deltaF=1E4;


Xmax=Xmax_micron*1E-6;
deltaX=deltaX_micron*1E-6;
N_X=Xmax/deltaX;



X=0:deltaX_micron:Xmax_micron;    %micron

Center=5;

Axial_Resolution=1.5;     %micron

PSF=gaussmf(X,[Lateral_Resolution/2/((2*log(2))^0.5) Center]);

plot(X,PSF);


Fmax=1/(deltaX*2);
F=0:deltaF:Fmax;
N_F=Fmax/deltaF+1;


S=(fft(PSF,N_F));



plot(F,abs(S));

