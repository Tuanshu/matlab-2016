clear all


Axial_ave_Factor=2;
Maximum_Axial_Frame=400;
Temp=0;
Axial_Length_Original=size(After_Npoint_Image_Stack,3);
Axial_Length_Used=floor(Axial_Length_Original/Axial_ave_Factor)*Axial_ave_Factor;
Reduced_Length=Axial_Length_Used/Axial_ave_Factor;
for p=1:Axial_ave_Factor
   Temp=Temp+After_Npoint_Image_Stack(:,:,(Axial_ave_Factor-(p-1)):Axial_ave_Factor:(Axial_ave_Factor*Reduced_Length)-(p-1));
end
Reduced_Stack=Temp/Axial_ave_Factor;
Reduced_Image=squeeze(mean(Reduced_Stack,2))';
Reduced_Image=Reduced_Image(1:min(size(Reduced_Image,1),Maximum_Axial_Frame),:);