clear all

Wavelength=0.56;
h=6.62606957E-34;
C=3E8;


Row=2000;
Colomn=2000;
Byte_Skip=1024;


Pixel_Size=11;
It=10;  %micronS
Max_Electron=36000;
Max_ADU=4096;
QE=0.5;
folder_path='F:\P1.2 Test\160826_Focus_for FWC estimation\';
Image_1_file_name='1.28microsec_0mW_260fps.b16';
Image_2_file_name='1.28microsec_0.31mW_260fps.b16';


file_path=[folder_path file_name];

fin=fopen(file_path);
if Image_Number == 1
    Image=fread(fin,[648,488],'uint16');
elseif Image_Number == 2
    Image=fread(fin,[1008,1008],'uint16','b');
end
imagesc(Image);
colormap(gray);
axis equal
axis off
xlim([0 size(Image,2)]);
ylim([0 size(Image,1)]);

%%

if Image_Number == 1
    Center=[307 274];

    Radius=179;

    X=repmat([1:648]',[1 488]);
    Y=repmat([1:488],[648 1]);
    Mask=ones(648,488);

elseif Image_Number == 2
    Center=[570 449];

    Radius=126;

    X=repmat([1:1008]',[1 1008]);
    Y=repmat([1:1008],[1008 1]);
    Mask=ones(1008,1008);

elseif Image_Number == 3
    Radius=154+999;
    
end

Mask(((X-Center(1)).^2+(Y-Center(2)).^2)>Radius^2)=0;

Masked_Image=Image.*Mask;
BND=mean(Image(Mask==0));
Masked_Image=Masked_Image-BND;

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