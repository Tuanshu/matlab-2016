clear all

Max_Electron=20000;
Max_ADU=4096;
Wavelength=0.56;
QE=0.5;
h=6.62606957E-34;
C=3E8;
Pixel_Size=0.45;
Image_Number=2;
if Image_Number == 1
    It=50;  %micronS
    folder_path='D:\AMO\160809_in vivo OCT\160810_P1.1 Current Status_it50\';
    file_name='13_no inter.raw';
elseif Image_Number == 2
    It=500;
    folder_path='D:\AMO\160809_in vivo OCT\160804_P1.1 status record\';
    file_name='160804_both sample and ref_input ~14 mW_it 500.raw';

elseif Image_Number == 3
    It=50;
    folder_path='I:\Everday Experiements\160811_New P1.1 Test\2016_0811_1327_Car_C392_it50_inter\';
    file_name='00000000';
    
end
file_path=[folder_path file_name];

fin=fopen(file_path);
Image=fread(fin,[648,488],'uint16');
imagesc(Image);
colormap(gray);
axis equal
axis off
xlim([0 size(Image,2)]);
ylim([0 size(Image,1)]);

%%
Center=[337 248];

if Image_Number == 1
    Radius=157.14286;

elseif Image_Number == 2
    Radius=999;
    
elseif Image_Number == 3
    Radius=154+999;
    
end
X=repmat([1:648]',[1 488]);
Y=repmat([1:488],[648 1]);

Mask=ones(648,488);
Mask(((X-Center(1)).^2+(Y-Center(2)).^2)>Radius^2)=0;

Masked_Image=Image.*Mask;
imagesc(Masked_Image);
colormap(gray);
axis equal
axis off
xlim([0 size(Image,2)]);
ylim([0 size(Image,1)]);

Integrated_ADU=sum(sum(Masked_Image));

Integrated_Electron=Integrated_ADU*Max_Electron/Max_ADU;

Estimated_Photon_Flux=Integrated_Electron/(It*1E-6)/QE;
Estimated_Optical_Power=Estimated_Photon_Flux*h*C/(Wavelength*1E-6)
Diameter=Radius*2*Pixel_Size