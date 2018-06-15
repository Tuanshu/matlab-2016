clear all

Header_Size=1024;
folder_path='I:\160713_iio image data\';
Data=1;

if Data == 1        %Deep
    file_path=[folder_path '011_015.bin'];
elseif Data == 2    %Wide
    file_path=[folder_path '20160713171718.bin'];
end

fin=fopen(file_path);
Header=fread(fin,Header_Size,'ulong');

if Header(2) == 0
    Colomn=488;
    Row=648;
else
    Row=Header(2);
    Colomn=Header(1);
end
    

fin=fopen(file_path);
if Header(257)>0
    fseek(fin, 256*4, 'bof');
else
    fseek(fin, Header_Size*4, 'bof');
end

if Data == 2
    Image=fread(fin,[Row,Colomn],'ulong');
elseif Data == 1
    Image_Temp=fread(fin,[Row,inf],'ulong'); %*Frame
    Frame=size(Image_Temp,2)/Colomn;
    for p=1:Frame
        Image(:,:,p)=Image_Temp(:,(1+(p-1)*Colomn):(p*Colomn));
    end
end


%%
Q=65;
Image_Read=Image(:,:,Q);
Array=zeros([Frame 1]);
for p=1:Frame
    Array(p)=mean(mean(Image(:,:,p),1),2);
end


imagesc(Image_Read);
colormap(gray);

Min=110;
Max=113;

Image_Norm=(Image_Read-Min*1E7)/((Max-Min)*1E7);
Image_Norm(Image_Norm>1)=1;

Image_Norm(Image_Norm<0)=0;


imagesc(Image_Norm);
axis equal
fclose('all');