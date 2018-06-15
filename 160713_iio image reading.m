clear all

Header_Size=1024;
folder_path='I:\160713_iio image data\';
Data=2;

if Data == 1        %Deep
    file_path=[folder_path '011_015.bin'];
elseif Data == 2    %Wide
    file_path=[folder_path '20160713171718.bin'];
end

fin=fopen(file_path);
Header=fread(fin,Header_Size,'ulong');

if Header(2) == 0
    Colomn=648;
    Row=488;
else
    Row=Header(2);
    Colomn=Header(1);
end
    

fin=fopen(file_path);
fseek(fin, Header_Size*4, 'bof');
if Data == 2
    Image=fread(fin,[Row,Colomn],'ulong');
elseif Data == 1
    Image_Temp=fread(fin,[Row,Colomn],'ulong'); %*Frame
    for p=1:1
        Image(:,:,p)=Image_Temp(:,(1+(p-1)*Colomn):(p*Colomn));
    end
end


%%
imagesc(Image(:,:,1));
colormap(gray);

Min=110;
Max=113;

Image_Norm=(Image-Min*1E7)/((Max-Min)*1E7);
Image_Norm(Image_Norm>1)=1;

Image_Norm(Image_Norm<0)=0;


imagesc(Image_Norm);
axis equal
fclose('all');