clear all

tic
Header_Size=256;

Depth_Sampling_Resolution=0.2;  %micron

First_index_range=[6 17];
Second_index_range=[8 23];

First_Length=First_index_range(2)-First_index_range(1)+1;
Second_Length=Second_index_range(2)-Second_index_range(1)+1;

Number_of_FOV=First_Length*Second_Length;

Colomn=488;
Row=648;

folder_path='J:\IIO Data\FOV Data\20160714170851_After 1st Z-compensation_S40_E80_1X_Again\';
file_list=dir(folder_path);


for p=1:First_Length
    for q=1:Second_Length
        file_path=[folder_path sprintf('%03d_%03d.bin',First_index_range(1)+p-1,Second_index_range(1)+q-1)];

        fin=fopen(file_path);
        Header=fread(fin,Header_Size,'ulong');


        fin=fopen(file_path);
        fseek(fin, 256*4, 'bof');

        Image_Temp=fread(fin,[Row,inf],'ulong'); %*Frame
        Frame=size(Image_Temp,2)/Colomn;
        
        if (p == 1) && (q == 1)
            Array_Map=zeros(First_Length,Second_Length,Frame);
        end
        
        Image_Stack=zeros(Row,Colomn,Frame);
        for r=1:Frame
            Image_Stack(:,:,r)=Image_Temp(:,(1+(r-1)*Colomn):(r*Colomn));
        end

        %%
        for r=1:Frame
            Array_Map(p,q,r)=mean(mean(Image_Stack(:,:,r),1),2);
        end
        disp(((p-1)*Second_Length+q)/Number_of_FOV);
    end
end
toc

dlmwrite([folder_path 'Array_Map.txt'],Array_Map,'delimiter','\t','newline','pc','precision', '%.6f');
%%
%Interface_Map=dlmread(['Interface_Map.txt']);

Relative_First_index_read=11;
Relative_Second_index_read=7;

plot(squeeze(Array_Map(Relative_First_index_read,Relative_Second_index_read,:)));

%% Z-map

[Max_Value_Map Max_Index_Map]=max(Array_Map,[],3);
imagesc(Max_Index_Map);

Interface_Map=Depth_Sampling_Resolution*Max_Index_Map;

First_index_Array=(First_index_range(1):First_index_range(2))';
First_index_Map=repmat(First_index_Array,[1 Second_Length]);

Second_index_Array=(Second_index_range(1):Second_index_range(2));
Second_index_Map=repmat(Second_index_Array,[First_Length 1]);


First_index_Map_Array=First_index_Map(:);
Second_index_Map_Array=Second_index_Map(:);
Interface_Map_Array=Interface_Map(:);

%% Math

XY=[First_index_Map_Array';Second_index_Map_Array';ones([1 length(First_index_Map_Array)])]';
Z=Interface_Map_Array;%-min(Z_max_Pos);

beta_matrix=XY\Z;
a=beta_matrix(1);
b=beta_matrix(2);
c=beta_matrix(3);

%%

Z_Map_Plane=a*First_index_Map+b*Second_index_Map+c;
Recidual_Map=Interface_Map-Z_Map_Plane;
All_First_index_Array=(0:23)';
All_First_index_Map=repmat(All_First_index_Array,[1 32]);

All_Second_index_Array=(0:31);
All_Second_index_Map=repmat(All_Second_index_Array,[24 1]);

All_Z_Map_Plane=a*All_First_index_Map+b*All_Second_index_Map+c;


All_Z_Map_Plane_Min=All_Z_Map_Plane-min(All_Z_Map_Plane(:));

imagesc(All_Z_Map_Plane_Min);

dlmwrite([folder_path 'Interface_Map.txt'],Interface_Map,'delimiter','\t','newline','pc','precision', '%.6f');
csvwrite([folder_path 'Interface_Map.csv'],Interface_Map);

csvwrite([folder_path 'compensate.csv'],All_Z_Map_Plane_Min);

Test=dlmread('compensate.csv');
Recidual_Map=Test-All_Z_Map_Plane_Min;
imagesc(Recidual_Map);
%%
Q=225;
Image_Read=Image_Stack(:,:,Q);


imagesc(Image_Read);
colormap(gray);

Min=106;
Max=110;

Image_Norm=(Image_Read-Min*1E7)/((Max-Min)*1E7);
Image_Norm(Image_Norm>1)=1;

Image_Norm(Image_Norm<0)=0;


imagesc(Image_Norm);
axis equal
fclose('all');

%% To add previous Map