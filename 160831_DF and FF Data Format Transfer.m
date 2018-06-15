clear all

Folder_Path='D:\MATLAB\160831_DF and FF Related\';

DF_Name='Dark_frame.txt';       %iio DF names: 0.bin, 1.bin, 2.bin, 3.bin, 4.bin

FF_Name='FF.txt';               %iio FF names: flat field.bin

Output_DF_Name='0.bin';
Output_FF_Name='flat field.bin';

DF=dlmread([Folder_Path DF_Name]);
FF=dlmread([Folder_Path FF_Name]);
%DF=zeros(size(FF,1),size(FF,2));

%FF=ones(size(FF,1),size(FF,2));
Repeat_Number=768;

%DF_Dummy=zeros(648,488);
%DF_Dummy(end,:)=1;

DF_Reshape=reshape(DF,[size(DF,1)*size(DF,2) 1]);
FF_Reshape=reshape(FF,[size(DF,1)*size(DF,2) 1]);



fid = fopen([Folder_Path Output_DF_Name], 'a+');
for p=1:Repeat_Number
    fwrite(fid, DF_Reshape, 'single');
end
fclose(fid);

fid = fopen([Folder_Path Output_FF_Name], 'a+');
for p=1:Repeat_Number
    fwrite(fid, FF_Reshape, 'single');
end
fclose(fid);
