clear all

frame_width=648;
frame_height=488;
number_of_frame=2000;
Averaging_Factor=18;

folder_path='I:\Everday Experiements\';
file_path=[folder_path '160712_Glass Interface Carrier 648-488-2000'];
fin=fopen(file_path);
A=fread(fin,[frame_width,frame_height*number_of_frame],'uint16','b');
A(isnan(A))=0;

Data_Volume=reshape(A,frame_width,frame_height,number_of_frame);

clear A
fclose('all');

%%
Q=1001;
Image(:,:)=Data_Volume(:,:,Q);
imagesc(Image);

%%
Averaged_Image_Array=zeros([number_of_frame 1 1]);
for p=1:number_of_frame
    Averaged_Image_Array(p)=mean(mean(Data_Volume(:,:,p),1),2);
end

plot(Averaged_Image_Array);
%% Image Averaging

Averaged_Data_Volume=zeros(frame_width,frame_height,floor(number_of_frame/Averaging_Factor));

for p=1:floor(number_of_frame/Averaging_Factor)
    Averaged_Data_Volume(:,:,p)=mean(Data_Volume(:,:,(1+(p-1)*Averaging_Factor):p*Averaging_Factor),3);
end
%%

Q=40;
Image(:,:)=Averaged_Data_Volume(:,:,Q);
imagesc(Image);
%%
Averaged_Averaged_Image_Array=zeros([floor(number_of_frame/Averaging_Factor) 1 1]);
for p=1:floor(number_of_frame/Averaging_Factor)
    Averaged_Averaged_Image_Array(p)=mean(mean(Averaged_Data_Volume(:,:,p),1),2);
end

plot(Averaged_Averaged_Image_Array);
%%

Amp=zeros(frame_width,frame_height,floor(number_of_frame/Averaging_Factor));
Phase=zeros(frame_width,frame_height,floor(number_of_frame/Averaging_Factor));


for p=1:floor(number_of_frame/Averaging_Factor)
    E1=Averaged_Data_Volume(:,:,p);
    
    if p+1>floor(number_of_frame/Averaging_Factor)
        E2=E1;
        E3=E1;
        E4=E1;

    elseif p+2>floor(number_of_frame/Averaging_Factor)
        E2=Averaged_Data_Volume(:,:,p+1);
        E3=E2;
        E4=E2;
    elseif p+3>floor(number_of_frame/Averaging_Factor)
        E2=Averaged_Data_Volume(:,:,p+1);
        E3=Averaged_Data_Volume(:,:,p+2);
        E4=E3;
    else
        E2=Averaged_Data_Volume(:,:,p+1);
        E3=Averaged_Data_Volume(:,:,p+2);
        E4=Averaged_Data_Volume(:,:,p+3);     
    end
    
    %Amp(:,:,p)=((E1-E2-E3+E4).^2+(E1-E2+E3-E4).^2).^0.5;   paper
    %Amp(:,:,p)=((E1-E2+E3+E4).^2+(E1-E2-E3-E4).^2).^0.5;
    Amp(:,:,p)=((4*(E1.^2+E2.^2+E3.^2+E4.^2)-(E1+E2+E3+E4).^2).^0.5)*(2^0.5)/4;
    Phase(:,:,p)=atan((E1-E2+E3+E4)./(E1-E2-E3-E4));

end
%zeros(frame_width,frame_height,number_of_frame);
% for q=1:number_of_frame
%     Data_Volume(:,:,q)=mean(mean(A(:,(frame_height*(q-1)+1):frame_height*q)));
% end
%%

Q=54;

Image1(:,:)=Amp(:,:,Q);
Image2(:,:)=Phase(:,:,Q);

subplot(1,2,1)
imagesc(Image1);

subplot(1,2,2)
imagesc(Image2);

%% Axial plot
Base=1:length(Amp_Axial);
X=504;
Y=341;

Amp_Axial=squeeze(Amp(X,Y,:));
Phase_Axial=squeeze(Phase(X,Y,:));

subplot(2,1,1)
plot(Base,Averaged_Averaged_Image_Array,Base,Amp_Axial);
subplot(2,1,2)

plot(Phase_Axial);
