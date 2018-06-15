clear all

Lateral_Resolution_Array=0.01:0.01:1.5;
MTF_Half_Sampling_Array=zeros([1 length(Lateral_Resolution_Array)]);
for p=1:length(Lateral_Resolution_Array)

    Lateral_Resolution=Lateral_Resolution_Array(p);         %micron

    Pixel_Size=0.4950;        %micron

    Xmax_mm=1;
    Pixel_Size_mm=Pixel_Size/1000;  %mm
    ZP_ratio=1;

    Simulated_Wavelength=0.01;

    deltaX_mm=0.00001;
    N_X=Xmax_mm/deltaX_mm;
    Lateral_Resolution_mm=Lateral_Resolution/1000;
    X=[0:deltaX_mm:Xmax_mm]-Xmax_mm/2;    %mm
    X_micron=X*1000;    %mm

    S=gaussmf(X,[Lateral_Resolution_mm/2/((2*log(2))^0.5) 0]);

    FWHM_Spacial=X_micron(find(S>0.5,1,'last'))-X_micron(find(S>0.5,1,'first'));    %micron


    plot(X_micron,S);
    xlim([-3 3]);
    xlabel('micron');
    disp(FWHM_Spacial);


    FFT_S=fft(S,N_X*ZP_ratio);

    MTF=abs(FFT_S)/max(abs(FFT_S));
    F=0:(1/Xmax_mm)/ZP_ratio:(1/Xmax_mm)/ZP_ratio*(length(FFT_S)-1);          %unit:cycle/mm

    F_Cycle_Pixel=F*Pixel_Size_mm/Xmax_mm;

    Index_Half_Sampling=find(F_Cycle_Pixel>0.5,1,'first');

    MTF_Half_Sampling=MTF(Index_Half_Sampling);
    MTF_Half_Sampling_Array(p)=MTF_Half_Sampling;

end
plot(F,MTF);
xlabel('Cycle/mm');


plot(F_Cycle_Pixel,MTF);
xlabel('Cycle/pixel');
xlim([0 1]);
disp(MTF_Half_Sampling);


plot(Lateral_Resolution_Array,MTF_Half_Sampling_Array,'linewidth',2);
xlabel(sprintf('Optical Resolution (\\mum) (Sampling Resolution = %.02f \\mum)',Pixel_Size));
ylabel('SFR_H_S');
xlim([Lateral_Resolution_Array(1) Lateral_Resolution_Array(end)]);
ylim([0 1]);

disp(MTF_Half_Sampling);
set(gca,'XTick', 0:0.1:100,'GridLineStyle',':'	,'MinorGridLineStyle',':', 'XColor', 'black');

grid on
