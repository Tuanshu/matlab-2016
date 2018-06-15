clear all

Xmax_mm=1;

Simulated_Wavelength=0.02;

deltaX_mm=0.001;
N_X=Xmax_mm/deltaX_mm;

X=0:deltaX_mm:Xmax_mm;    %micron

S_1mm=sin(2*pi*X/Simulated_Wavelength);

plot(X,S_1mm);

FFT_S=fft(S_1mm);

F=0:(1/Xmax_mm):(1/Xmax_mm)*(length(FFT_S)-1);          %unit:cycle/mm

plot(F,abs(FFT_S));
xlabel('Cycle/mm');

