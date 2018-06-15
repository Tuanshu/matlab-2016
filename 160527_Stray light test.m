clear all

% Data reading & stacking
Folder_Path='D:\160527_PD stray light test\';

withoutPDwPDlens_File_Name='withoutPDwPDlens-3500';
withoutPDwithoutPDlens_File_Name='withoutPDwithoutPDlens-3500';
wPDwPDlens_File_Name='wPDwPDlens-3500';
wPDwithoutPDlens_File_Name='wPDwithoutPDlens-3500_ADC47';


Input_Image_withoutPDwPDlens_File=imread([Folder_Path withoutPDwPDlens_File_Name],'jpg')';
Input_Image_withoutPDwithoutPDlens=imread([Folder_Path withoutPDwithoutPDlens_File_Name],'jpg')';
Input_Image_wPDwPDlens=imread([Folder_Path wPDwPDlens_File_Name],'jpg')';
Input_Image_wPDwithoutPDlens=imread([Folder_Path wPDwithoutPDlens_File_Name],'jpg')';

wPDlens=Input_Image_wPDwPDlens-Input_Image_withoutPDwPDlens_File;
withoutPDlens=Input_Image_wPDwithoutPDlens-Input_Image_withoutPDwithoutPDlens;

imagesc(wPDlens);
axis equal
axis off
caxis([0 10]);

imagesc(withoutPDlens);
axis equal
axis off
caxis([0 10]);
