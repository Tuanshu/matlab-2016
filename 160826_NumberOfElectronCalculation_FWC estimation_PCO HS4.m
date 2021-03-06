clear all

Wavelength=0.56;
h=6.62606957E-34;
C=3E8;


Row=2000;
Colomn=2000;
Byte_Skip=1024;


Pixel_Size=11;
It=2.144;%3.872;  %micronS
Max_Electron=27750;
Max_ADU=4096;
QE=0.5;

Power_1=0;
Power_2=0.31;%0.153;   %mW

Power_Diff=Power_2-Power_1;

folder_path='G:\P1.2 Test\160826_Focus_for FWC estimation\';
Image_1_file_name=sprintf('%gmicrosec_%gmW_260fps.b16',It,Power_1);
Image_2_file_name=sprintf('%gmicrosec_%gmW_260fps.b16',It,Power_2);


file_path_1=[folder_path Image_1_file_name];
file_path_2=[folder_path Image_2_file_name];

fin=fopen(file_path_1);
fseek(fin, Byte_Skip, 'bof');
Image_1=fread(fin,[Row,inf],'uint16')/16; %*Frame   不知為何, 看起來就像是要除16


imagesc(Image_1);

fin=fopen(file_path_2);
fseek(fin, Byte_Skip, 'bof');
Image_2=fread(fin,[Row,inf],'uint16')/16; %*Frame   不知為何, 看起來就像是要除16

imagesc(Image_2);

Image_Diff=Image_2-Image_1;

imagesc(Image_Diff);
colormap(gray);
axis equal
axis off
xlim([0 size(Image_Diff,2)]);
ylim([0 size(Image_Diff,1)]);

%%

Integrated_ADU=sum(Image_Diff(:));

Integrated_Electron=Integrated_ADU*Max_Electron/Max_ADU;

Estimated_Photon_Flux=Integrated_Electron/(It*1E-6)/QE;
Estimated_Optical_Power=Estimated_Photon_Flux*h*C/(Wavelength*1E-6)

Estimated_FWC=(Power_Diff*1E-3)/(Integrated_ADU/Max_ADU)/(h*C/(Wavelength*1E-6))*QE*(It*1E-6)
