clear all

Xmax_mm=1;

ZP_ratio=2;

Simulated_Wavelength=0.01;

deltaX_mm=0.001;
N_X=Xmax_mm/deltaX_mm;

X=0:deltaX_mm:Xmax_mm;    %micron

S=sin(2*pi*X/Simulated_Wavelength);

plot(X,S);

FFT_S=fft(S,N_X*ZP_ratio);

F=0:(1/Xmax_mm)/ZP_ratio:(1/Xmax_mm)/ZP_ratio*(length(FFT_S)-1);          %unit:cycle/mm

plot(F,abs(FFT_S));
xlabel('Cycle/mm');

