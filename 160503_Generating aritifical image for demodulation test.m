clear all

N=4;

Output_Raw_Path='D:\160503_Output Raw\';
Name='Test';

Image_Size=[648 488];

Number_of_Frame=250*4;


Averaging_Factor=16;

Center_Frame=Number_of_Frame/2;

Frame_Rate=250;


Lateral_Resolution=1.35;

Frame_Spacing=0.2/64;


Signal_Wavelength=0.56/1.4; %micron

Axial_Resolution=1.5/1.4;


Frame_Position_Array=(Frame_Spacing:Frame_Spacing:Frame_Spacing*Number_of_Frame)-Center_Frame*Frame_Spacing;

%% PSF generation (固定把PSF的中心放在Number_of_Frame的中心


Carrier=cos(4*pi*Frame_Position_Array/Signal_Wavelength);
Envelope=exp(-1*log(2)*(Frame_Position_Array/(Axial_Resolution/2)).^2);
PSF=Carrier.*Envelope;
plot(Frame_Position_Array,Carrier,Frame_Position_Array,Envelope);
plot(Frame_Position_Array,PSF);

%%

Signal_Circle_Lateral_Position=[150 150];
Signal_Circle_Size=145;


Noise_Circle_Lateral_Position_Original=[Image_Size(1)-Signal_Circle_Lateral_Position(1) Image_Size(2)-Signal_Circle_Lateral_Position(2)];
Noise_Circle_Size=145;

Noise_Circle_Lateral_Max_Movement=[1 1]/2;       %pixel

Noise_Circle_Lateral_Movement_Period=20;        %frame #



%%
X_Grid=repmat(1:Image_Size(1),[Image_Size(2) 1])';
Y_Grid=repmat([1:Image_Size(2)]',[1 Image_Size(1)])';


Single_Frame=zeros([Image_Size(1) Image_Size(2)]);

Single_Frame(((X_Grid-Signal_Circle_Lateral_Position(1)).^2+(Y_Grid-Signal_Circle_Lateral_Position(2)).^2)<Signal_Circle_Size^2)=1;

Data_Volume=zeros([Image_Size(1) Image_Size(2) Number_of_Frame]);

Current_Noise_Circle_Lateral_Position=Noise_Circle_Lateral_Position_Original;
Blank_Frame=zeros([Image_Size(1) Image_Size(2)]);

for p=1:Number_of_Frame
    Current_Noise_Frame=Blank_Frame;
    Current_Noise_Frame((X_Grid-Current_Noise_Circle_Lateral_Position(1)).^2+(Y_Grid-Current_Noise_Circle_Lateral_Position(2)).^2<Noise_Circle_Size^2)=1;
    Data_Volume(:,:,p)=Single_Frame.*PSF(p)+Current_Noise_Frame;
    Current_Noise_Circle_Lateral_Position=Current_Noise_Circle_Lateral_Position+Noise_Circle_Lateral_Max_Movement*cos(2*pi*p/Noise_Circle_Lateral_Movement_Period);
    disp(p);
end

%% Ave
Data_Volume_Ave=zeros(Image_Size(1),Image_Size(2),floor(Number_of_Frame/Averaging_Factor));
Frame_Position_Array_Ave=zeros([floor(Number_of_Frame/Averaging_Factor) 1 1]);
for p=1:floor(Number_of_Frame/Averaging_Factor)
    Data_Volume_Ave(:,:,p)=mean(Data_Volume(:,:,(1+(p-1)*Averaging_Factor):p*Averaging_Factor),3);
    Frame_Position_Array_Ave(p)=mean(Frame_Position_Array((1+(p-1)*Averaging_Factor):p*Averaging_Factor));
end

clear Data_Volume

%%
NNN=50;

imagesc(Data_Volume_Ave(:,:,NNN));

axis equal
axis off
caxis([0 1]);
colormap(gray);

%%
fout=fopen([Output_Raw_Path,sprintf('%s %d-%d-%d',Name,Image_Size(1),Image_Size(2),Number_of_Frame)],'w+');
fwrite(fout,Data_Volume_Ave,'float32','b');

%%
Index=1:N;
Index_Choice=nchoosek(Index,2);

%%
Processed_Data_Volume=zeros([Image_Size(1) Image_Size(2) floor(Number_of_Frame/4/Averaging_Factor)]);

for p=1:size(Processed_Data_Volume,3)
    Raw=Data_Volume_Ave(:,:,(1+(p-1)*N):p*N);
    Temp=0;
    for q=1:length(Index_Choice)
        Temp=Temp+(Raw(:,:,Index_Choice(q,1))-Raw(:,:,Index_Choice(q,2))).^2;
    end
    Processed_Data_Volume(:,:,p)=Temp.^0.5;
    disp(p);
end

%%


fout=fopen([Output_Raw_Path,sprintf('Processed_%s %d-%d-%d',Name,Image_Size(1),Image_Size(2),Number_of_Frame)],'w+');
fwrite(fout,Processed_Data_Volume,'float32','b');

%% Record max vibration noise
L_Processed_Data_Volume=Processed_Data_Volume(1:floor(Image_Size(1)/2),:,:);
R_Processed_Data_Volume=Processed_Data_Volume((floor(Image_Size(1)/2)+1):end,:,:);
Max_L_Processed=max(L_Processed_Data_Volume(:));
Max_R_Processed=max(R_Processed_Data_Volume(:));
Mean_L_Processed=mean(L_Processed_Data_Volume(:));
Mean_R_Processed=mean(R_Processed_Data_Volume(:));
clear L_Processed_Data_Volume R_Processed_Data_Volume
%%
SB=11;
LB=21;

Test_Array=squeeze(Data_Volume_Ave(Signal_Circle_Lateral_Position(1),Signal_Circle_Lateral_Position(2),:));
plot(Test_Array);

FFT_Test_Array=fft(Test_Array,length(Test_Array));
plot(real(FFT_Test_Array));

FFT_Test_Array(1:SB)=0;
FFT_Test_Array(LB:end)=0;

Test_Array_Filtered=abs(ifft(FFT_Test_Array));


plot(Test_Array_Filtered);

%%
clear Processed_Data_Volume 
FFT_Data_Volume=fft(Data_Volume_Ave,[],3);

FFT_Data_Volume(:,:,1:SB)=0;
FFT_Data_Volume(:,:,LB:end)=0;

Data_Volume_Filtered=abs(ifft(FFT_Data_Volume,[],3));
%%
fout=fopen([Output_Raw_Path,sprintf('Filtered_%s %d-%d-%d',Name,Image_Size(1),Image_Size(2),Number_of_Frame)],'w+');
fwrite(fout,Data_Volume_Filtered,'float32','b');

%%

L_Filtered_Data_Volume=Data_Volume_Filtered(1:floor(Image_Size(1)/2),:,:);
R_Filtered_Data_Volume=Data_Volume_Filtered((floor(Image_Size(1)/2)+1):end,:,:);
Max_L_Filtered=max(L_Filtered_Data_Volume(:));
Max_R_Filtered=max(R_Filtered_Data_Volume(:));
Mean_L_Filtered=mean(L_Filtered_Data_Volume(:));
Mean_R_Filtered=mean(R_Filtered_Data_Volume(:));
clear L_Filtered_Data_Volume R_Filtered_Data_Volume

SNR_Processed=Max_L_Processed/Max_R_Processed;
SNR_Filtered=Max_L_Filtered/Max_R_Filtered;


SNR_Processed_Mean=Mean_L_Processed/Mean_R_Processed;
SNR_Filtered_Mean=Mean_L_Filtered/Mean_R_Filtered;

plot(Frame_Position_Array,PSF/max(PSF),Frame_Position_Array_Ave,Test_Array_Filtered/max(Test_Array_Filtered));
