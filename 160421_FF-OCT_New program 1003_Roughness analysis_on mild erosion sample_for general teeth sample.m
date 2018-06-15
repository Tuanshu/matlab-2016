clear all
tic
% Data reading & stacking
set=1;
If_schematic=0;
If_select_ROI=0;
If_generate_hex_map=0;
If_Analysis=0;

Number_of_Superpixel=100;


Training_Data_Storage_Folder_Path='D:\160421_Teeth_Trianing\\';

if set == 1     %30 sec, DEJ, good, important
    folder_path='G:\Everday Experiements\151229_MMH\151229_Temp_28th scan\';
    division_number=2:30;
    NNN=136;
    
    Start_Searching_Index=60-2*16;
    End_Searching_Index=200-2*16;
    
    First_index=15;
    Second_index=15;
    Third_index=15;
    
    
    ROI=[1 648 1 488];
    Class_index_array=[1 2];    %1: 0sec, 2: 30sec, 3: 60sec

    
elseif set == 2     %30 sec, detine only, good result, but less interest
    folder_path='G:\Everday Experiements\151229_MMH\151229_Temp_29th scan\';
    division_number=2:19;
    NNN=136;
    
    Start_Searching_Index=60-2*16;
    End_Searching_Index=200-2*16;
    
    First_index=15;
    Second_index=15;
    Third_index=15;
    

    
    ROI=[1 648 1 488];
    Class_index_array=[1 2];    %1: 0sec, 2: 30sec, 3: 60sec

    
    
    
elseif set == 3 %30 sec, enamel only, not so obvious, near boundary
    folder_path='G:\Everday Experiements\151229_MMH\151229_Temp_30th scan\';
    division_number=2:20;
    NNN=436;
    
    Start_Searching_Index=70-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=15;
    Second_index=15;
    Third_index=15;
    
    
    ROI=[1 648 1 488];
    Class_index_array=[1 2];    %1: 0sec, 2: 30sec, 3: 60sec

    

elseif set == 4 %30 sec, DEJ, only enamel were cut, good
    folder_path='G:\Everday Experiements\151229_MMH\151229_Temp_31th scan\';
    division_number=2:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=200;
    
    First_index=15;
    Second_index=15;
    Third_index=15;
    
    ROI=[1 648 1 488];
    Class_index_array=[1 2];    %1: 0sec, 2: 30sec, 3: 60sec

    

elseif set == 5 %30 sec, DEJ, bad
    folder_path='G:\Everday Experiements\151228_MMH\151228_2nd scan_Sampl 3_Temp_Deep100_DEJ\';
    division_number=2:15;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=15;
    Second_index=15;
    Third_index=15;
    
    ROI=[1 648 1 488];
    Class_index_array=[1 2];    %1: 0sec, 2: 30sec, 3: 60sec

    
elseif set == 6 %上一個轉90度, 也看不出來, 但好像有一些膠帶痕跡, 我猜是沒貼好
    folder_path='G:\Everday Experiements\151228_MMH\151228_3rd scan_Sampl 3_Temp_Deep100_DEJ rotate 90\';
    division_number=2:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=15;
    Second_index=15;
    Third_index=15;
    

    ROI=[1 648 1 488];
    Class_index_array=[1 2];    %1: 0sec, 2: 30sec, 3: 60sec

    
elseif set == 7  %在6的旁邊, 好像蠻明顯是沒貼好的
    folder_path='G:\Everday Experiements\151228_MMH\151228_4th scan_Sampl 3_Temp_Deep100_Dentine_rorate 90\';
    division_number=2:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=15;
    Second_index=15;
    Third_index=15;
    
    
    ROI=[1 648 1 488];
    Class_index_array=[1 2];    %1: 0sec, 2: 30sec, 3: 60sec

elseif set == 8  %60sec, DEJ, good, 但好像是內層
    folder_path='G:\Everday Experiements\151228_MMH\151228_5th scan_Sampl 2_Temp_Deep100_DEJ\';
    division_number=2:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=15;
    Second_index=15;
    Third_index=15;
    

    ROI=[1 648 1 488];
    Class_index_array=[1 3];    %1: 0sec, 2: 30sec, 3: 60sec

    
    
elseif set == 9  %60sec, DEJ, 好像D比E清楚
    folder_path='G:\Everday Experiements\151228_MMH\151228_6th scan_Sampl 2_Temp_Deep100_DEJ\';
    division_number=2:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=15;
    Second_index=15;
    Third_index=15;
    
    
    ROI=[1 648 1 488];

    
    ROI=[1 648 1 488];
    Class_index_array=[1 3];    %1: 0sec, 2: 30sec, 3: 60sec

    
elseif set == 10  %60sec, repeat previous, better, DEJ
    folder_path='G:\Everday Experiements\151228_MMH\151228_7th scan_Sampl 2_Temp_Deep100_DEJ\';
    division_number=2:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=15;
    Second_index=15;
    Third_index=15;
    

    
    
    ROI=[1 648 1 488];
    Class_index_array=[1 3];    %1: 0sec, 2: 30sec, 3: 60sec

    
    
elseif set == 11  %60 sec, DEJ, good
    folder_path='G:\Everday Experiements\151228_MMH\151228_8th scan_Sampl 2_Temp_Deep100_DEJ\';
    division_number=5:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=15;
    Second_index=15;
    Third_index=15;
    
    
    ROI=[1 648 1 488];
    Class_index_array=[1 3];    %1: 0sec, 2: 30sec, 3: 60sec

    
    
elseif set == 12  %60 sec, DEJ, good
    folder_path='G:\Everday Experiements\151228_MMH\151228_10th scan_Sampl 2_Temp_Deep100_DEJ\';
    division_number=5:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=15;
    Second_index=15;
    Third_index=15;
    

    ROI=[1 648 1 488];
    Class_index_array=[1 3];    %1: 0sec, 2: 30sec, 3: 60sec

    

    
elseif set == 13  %240 sec, 但反而沒那麼清楚, 應該是全畫面皆enamel
    folder_path='G:\Everday Experiements\151228_MMH\151228_11th scan_Sampl 1_Temp_Deep100_0-4\';
    division_number=0:20;
    NNN=86;
    Start_Searching_Index=64;
    End_Searching_Index=180-2*16;
    
    First_index=15;
    Second_index=15;
    Third_index=15;
    

    ROI=[1 648 1 488];
    Class_index_array=[1 4];    %1: 0sec, 2: 30sec, 3: 60sec, 4: 240 sec

    
    
    
elseif set == 14  %240 sec, DEJ, 但也不清楚
    folder_path='G:\Everday Experiements\151228_MMH\151228_12th scan_Sampl 1_Temp_Deep100_DEJ_4-0\';
    division_number=0:20;
    NNN=86;
    Start_Searching_Index=59;
    End_Searching_Index=180-2*16;
    
    First_index=15;
    Second_index=15;
    Third_index=15;
    

    
    ROI=[1 648 1 488];
    Class_index_array=[1 4];    %1: 0sec, 2: 30sec, 3: 60sec, 4: 240sec

    
    
elseif set == 15  %120 sec, 他說有NaN, 那算了
    folder_path='G:\Everday Experiements\151228_MMH\151228_13th scan_Sampl 1_Temp_Deep100_2-0\';
    division_number=0:20;
    NNN=86;
    Start_Searching_Index=64;
    End_Searching_Index=180-2*16;
    
    First_index=15;
    Second_index=15;
    Third_index=15;
    

end
%parts=regexp(folder_path_without_index, '\','split');
%Parent_folder=folder_path_without_index(1:(length(folder_path_without_index)-length(parts{end})));

cmin=12;%12;%11
cmax_set=16;%16;%60;%120;%20;
X_offset=10;
Y_offset=30;

frame_width=648;
frame_height=488;


num_of_frame_per_division=16;

normalization_factor=50;


%%
temp_frame_volume=zeros(frame_width,frame_height,num_of_frame_per_division*length(division_number));
cd(folder_path);
for NN=1:length(division_number)
    file_path=[folder_path sprintf('%08d',division_number(NN))];
    fin=fopen(file_path);
    A=fread(fin,[frame_width,frame_height*num_of_frame_per_division],'float32','b');
    if fin ==-1
        k=k+1;
        fclose('all');
    else
        for q=1:num_of_frame_per_division
        temp_frame_volume(:,:,(NN-1)*num_of_frame_per_division+q)=A(:,(frame_height*(q-1)+1):frame_height*q);
        end               
        fclose('all');
    end
    disp(NN);
end

%% B-scan
cmin=8;%12;%11
cmax_set=30;%16;%60;%120;%20;

B_scan(:,:)=temp_frame_volume(:,NNN,:);
B_scan_Adjusted=(B_scan'-cmin)/cmax_set;

subplot (2,1,1)

imagesc(B_scan_Adjusted);
    colormap('gray');
    caxis([0 1]);
    axis equal
    xlim([0 size(B_scan_Adjusted,2)]);
    ylim([0 size(B_scan_Adjusted,1)]);

%% Volume Reducing
Reduced_Volume=permute(temp_frame_volume(:,:,Start_Searching_Index:End_Searching_Index),[2 1 3]);

    %%
NNN=153;
B_scan_Reduced(:,:)=Reduced_Volume(:,NNN,:);
B_scan_Reduced_Adjusted=(B_scan_Reduced'-cmin)/cmax_set;
    

%% 3D volume opening
Sphere_Radius=2;
Element=ones(Sphere_Radius*2+1,Sphere_Radius*2+1,Sphere_Radius*2+1);
Xgrid=ones(size(Element,1),size(Element,2),size(Element,3));
Ygrid=ones(size(Element,1),size(Element,2),size(Element,3));
Zgrid=ones(size(Element,1),size(Element,2),size(Element,3));

for p=1:size(Element,2)
    Xgrid(p,:,:)=p;
end
for q=1:size(Element,1)
    Ygrid(:,q,:)=q;
end
for w=1:size(Element,3)
    Zgrid(:,:,w)=w;
end

Element((((Xgrid-(Sphere_Radius+1)).^2+(Ygrid-(Sphere_Radius+1)).^2+(Zgrid-(Sphere_Radius+1)).^2).^0.5)>Sphere_Radius)=0;
    %%
    %Q=5;
    %imagesc(Element(:,:,Q));
    
%%
Opened_Volume=imopen(Reduced_Volume,Element);
Closed_Volume=imclose(Reduced_Volume,Element);

%% B-scan Opened
cmin=8;%12;%11
cmax_set=80;%16;%60;%120;%20;
    
    
%% 3D volume opening along Z only
Sphere_Radius=2;
Element_Z=ones(1,1,Sphere_Radius*2+1);
Xgrid=ones(size(Element_Z,1),size(Element_Z,2),size(Element_Z,3));
Ygrid=ones(size(Element_Z,1),size(Element_Z,2),size(Element_Z,3));
Zgrid=ones(size(Element_Z,1),size(Element_Z,2),size(Element_Z,3));

for p=1:size(Element_Z,2)
    Xgrid(p,:,:)=p;
end
for q=1:size(Element_Z,1)
    Ygrid(:,q,:)=q;
end
for w=1:size(Element_Z,2)
    Zgrid(:,:,w)=w;
end

Element_Z((((Zgrid-(Sphere_Radius+1)).^2).^0.5)>Sphere_Radius)=0;
    %%
    %Q=10;
    %imagesc(Element_Z(:,:,Q));
    
%%
Opened_Z_Volume=imopen(Reduced_Volume,Element_Z);
%% Peak searching 
[max_value max_index_map_Reduced]=max(Reduced_Volume(:,:,:),[],3);
[max_value max_index_map_Opened_Reduced]=max(Opened_Volume(:,:,:),[],3);
[max_value max_index_map_Closed_Reduced]=max(Closed_Volume(:,:,:),[],3);
[max_value max_index_map_Opened_Z_Reduced]=max(Opened_Z_Volume(:,:,:),[],3);

max_index_map=max_index_map_Reduced+Start_Searching_Index-1;
max_index_map_Opened=max_index_map_Opened_Reduced+Start_Searching_Index-1;
max_index_map_Closed=max_index_map_Closed_Reduced+Start_Searching_Index-1;
max_index_map_Opened_Z=max_index_map_Opened_Z_Reduced+Start_Searching_Index-1;

    %%
Height_map=max_index_map*0.2;

Height_map_Closed=max_index_map_Closed*0.2;

%% Z separation
First_Zone=zeros(size(Reduced_Volume,1),size(Reduced_Volume,2),First_index);
Second_Zone=zeros(size(Reduced_Volume,1),size(Reduced_Volume,2),Second_index);
Third_Zone=zeros(size(Reduced_Volume,1),size(Reduced_Volume,2),Third_index);

for p=1:size(Reduced_Volume,1)
    for q=1:size(Reduced_Volume,2)
        if max_index_map_Closed_Reduced(p,q)+First_index+Second_index+Third_index<size(Reduced_Volume,3)
            First_Zone(p,q,:)=Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+1):(max_index_map_Closed_Reduced(p,q)+First_index));
            Second_Zone(p,q,:)=Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+First_index+1):(max_index_map_Closed_Reduced(p,q)+First_index+Second_index));
            Third_Zone(p,q,:)=Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+First_index+Second_index+1):(max_index_map_Closed_Reduced(p,q)+First_index+Second_index+Third_index));
        elseif max_index_map_Closed_Reduced(p,q)+First_index+Second_index<size(Reduced_Volume,3)
            First_Zone(p,q,:)=Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+1):(max_index_map_Closed_Reduced(p,q)+First_index));
            Second_Zone(p,q,:)=Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+First_index+1):(max_index_map_Closed_Reduced(p,q)+First_index+Second_index));
            Third_Zone(p,q,1:size(Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+First_index+Second_index+1):end),3))=Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+First_index+Second_index+1):end);
        elseif max_index_map_Closed_Reduced(p,q)+First_index<size(Reduced_Volume,3)
            First_Zone(p,q,:)=Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+1):(max_index_map_Closed_Reduced(p,q)+First_index));
            Second_Zone(p,q,1:size(Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+First_index+1):end),3))=Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+First_index+1):end);
        else
            First_Zone(p,q,1:size(Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+1):end),3))=Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+1):end);
    
        end
    end
end
First_Map=mean(First_Zone,3);
Second_Map=mean(Second_Zone,3);
Third_Map=mean(Third_Zone,3);
Total_Map=(First_Map*First_index+Second_Map*Second_index+Third_Map*Third_index)/(First_index+Second_index+Third_index);


Ratio_Map=mean(Second_Zone,3)./mean(First_Zone,3);
Ratio_Map_2=mean(Third_Zone,3)./mean(Second_Zone,3);

Ratio_Map=Ratio_Map(ROI(3):ROI(4),ROI(1):ROI(2));
Ratio_Map_2=Ratio_Map_2(ROI(3):ROI(4),ROI(1):ROI(2));


subplot (2,2,1)


imagesc(First_Map(:,:));
%caxis([0 2]);
colormap('gray');
axis equal
xlim([0 size(First_Map,2)]);
ylim([0 size(First_Map,1)]);
axis off


subplot (2,2,2)


imagesc(Second_Map(:,:));
%caxis([0 100]);
colormap('gray');
axis equal
xlim([0 size(Second_Map,2)]);
ylim([0 size(Second_Map,1)]);
axis off

subplot (2,2,3)


imagesc(Third_Map(:,:));
%caxis([0 2]);
colormap('gray');
axis equal
xlim([0 size(Third_Map,2)]);
ylim([0 size(Third_Map,1)]);
axis off


subplot (2,2,4)


imagesc(Total_Map(:,:));
%caxis([0 100]);
colormap('gray');
axis equal
xlim([0 size(Total_Map,2)]);
ylim([0 size(Total_Map,1)]);
axis off


%% 生成hex grid
if If_generate_hex_map == 1
    rows=size(Ratio_Map,1);       %image size x
    cols=size(Ratio_Map,2);       %image size y
    k=Number_of_Superpixel;          %number of superpixel

    S = sqrt(rows*cols / (k * sqrt(3)/2));
    nodeCols = round(cols/S - 0.5);
    S = cols/(nodeCols + 0.5); 

    nodeRows = round(rows/(sqrt(3)/2*S));
    vSpacing = rows/nodeRows;

    k = nodeRows * nodeCols;

    C = zeros(2,k);          % Cluster centre data  1:3 is mean Lab value,
                             % 4:5 is row, col of centre, 6 is No of pixels
    l = -ones(rows, cols);   % Pixel labels.
    d = inf(rows, cols);     % Pixel distances from cluster centres.

    kk = 1;
    r = vSpacing/2;
    %% 
    for ri = 1:nodeRows

        if mod(ri,2), c = S/2; else, c = S;  end

        for ci = 1:nodeCols
            cc = round(c); rr = round(r);
            C(:, kk) = [cc; rr];
            c = c+S;
            kk = kk+1;
        end

        r = r+vSpacing;
    end

    S = round(S);  % We need S to be an integer from now on

    for kk = 1:k  % for each cluster

       % Get subimage around cluster
       rmin = max(C(2,kk)-S, 1);   rmax = min(C(2,kk)+S, rows); 
       cmin = max(C(1,kk)-S, 1);   cmax = min(C(1,kk)+S, cols); 

        [x,y] = meshgrid(cmin:cmax, rmin:rmax);
        x = x-C(1,kk);  % x and y dist from cluster centre
        y = y-C(2,kk);
        D = sqrt((x.^2 + y.^2)/S^2);

       subd =  d(rmin:rmax, cmin:cmax);
       subl =  l(rmin:rmax, cmin:cmax);
       updateMask = D < subd;
       subd(updateMask) = D(updateMask);
       subl(updateMask) = kk;

       d(rmin:rmax, cmin:cmax) = subd;
       l(rmin:rmax, cmin:cmax) = subl;   
    end
    imagesc(l);
    Hex_Grid=l;
    dlmwrite([Training_Data_Storage_Folder_Path,sprintf('Number_%d_Hex_Grid.txt',Number_of_Superpixel)],Hex_Grid,'delimiter','\t','newline','pc','precision', '%.6f');

else
    Hex_Grid=dlmread([Training_Data_Storage_Folder_Path,sprintf('Number_%d_Hex_Grid.txt',Number_of_Superpixel)]);
    kk=max(Hex_Grid(:));
end    
%% 選擇ROI
if If_select_ROI == 1
    ROI_Map=zeros(size(Hex_Grid,1),size(Hex_Grid,2),length(Class_index_array));
    ROI_Map_Hex=ROI_Map;
    Fig_ROI=figure;
    for p=1:length(Class_index_array)
       ROI_Map(:,:,p)=roipoly(Ratio_Map);
    end
    close(Fig_ROI);
    %%
    for p=1:length(Class_index_array)
        ROI_Map_Temp=ROI_Map(:,:,p);
        ROI_Map_Hex_Temp=zeros(size(Hex_Grid,1),size(Hex_Grid,2));
        for q=1:kk
            Current_Superpixel_BW=ROI_Map_Temp(Hex_Grid==q);

            ROI_Map_Hex_Temp(Hex_Grid==q)=min(Current_Superpixel_BW(:));
        end
        ROI_Map_Hex(:,:,p)=ROI_Map_Hex_Temp*Class_index_array(p);
        disp(p);
    end
    
    ROI_Map_Hex=sum(ROI_Map_Hex,3);

    imagesc(ROI_Map_Hex);



    dlmwrite([Training_Data_Storage_Folder_Path,sprintf('Set_%d_Class_index_array.txt',set)],Class_index_array,'delimiter','\t','newline','pc','precision', '%.6f');
    dlmwrite([Training_Data_Storage_Folder_Path,sprintf('Number_%d_Set_%d_ROI_Map_Hex.txt',Number_of_Superpixel,set)],ROI_Map_Hex,'delimiter','\t','newline','pc','precision', '%.6f');

else
    Class_index_array=dlmread([Training_Data_Storage_Folder_Path,sprintf('Set_%d_Class_index_array.txt',set)]);
    ROI_Map_Hex=dlmread([Training_Data_Storage_Folder_Path,sprintf('Number_%d_Set_%d_ROI_Map_Hex.txt',Number_of_Superpixel,set)]);

end
%% 生成上述ROI map hex之邊界(for 示意圖)
Line_thickness=2;

Element=ones(Line_thickness*2+1,Line_thickness*2+1);
ROI_Map_Boundary=imdilate(ROI_Map_Hex,Element)-ROI_Map_Hex;

imagesc(ROI_Map_Boundary);


%% 生成每個superpixel之ratio以及X與Y map, 同時, 檢查該superpixel是否屬於某個selected ROI (group), 
%% 若是, 則將該superpixel之index, First value, Second value紀錄並存到cell中

[Y_grid X_grid]=meshgrid(1:size(Ratio_Map,2),1:size(Ratio_Map,1));
Ratio_Map_Hex=Ratio_Map;
First_Map_Hex=First_Map;
Second_Map_Hex=Second_Map;
X_grid_Hex=X_grid;
Y_grid_Hex=Y_grid;

Data_Array=cell([length(Class_index_array) 1]);


for p=1:kk
    
    %Current_Superpixel=Ratio_Map(Hex_Grid==p);
    Current_Superpixel_First=First_Map(Hex_Grid==p);
    Current_Superpixel_Second=Second_Map(Hex_Grid==p);
    
    Current_Superpixel_ROI=max(ROI_Map_Hex(Hex_Grid==p));

    %Current_Superpixel(isnan(Current_Superpixel))=[];
    %Current_Superpixel(isinf(Current_Superpixel))=[];
    
    Current_Superpixel_First(isnan(Current_Superpixel_First))=[];
    Current_Superpixel_First(isinf(Current_Superpixel_First))=[];
    
    Current_Superpixel_Second(isnan(Current_Superpixel_Second))=[];
    Current_Superpixel_Second(isinf(Current_Superpixel_Second))=[];

    
    %Current_Data=[mean(Current_Superpixel_First) mean(Current_Superpixel_Second)];
    Current_Data=[max(Current_Superpixel_First) max(Current_Superpixel_Second)];

    
    First_Map_Hex(Hex_Grid==p)=Current_Data(1);
    Second_Map_Hex(Hex_Grid==p)=Current_Data(2);
    %Ratio_Map_Hex(Hex_Grid==p)=Current_Data(3);

    if Current_Superpixel_ROI>0
        for q=1:length(Class_index_array)
            if Class_index_array(q)==Current_Superpixel_ROI
                Data_Array{q}=[Data_Array{q}; Current_Data];
            end
        end
    end
    

    disp(p);
end
subplot(1,2,1)
imagesc(First_Map_Hex);
colormap('gray');
axis equal
xlim([0 size(First_Map_Hex,2)]);
ylim([0 size(First_Map_Hex,1)]);
axis off

subplot(1,2,2)
imagesc(Second_Map_Hex);
colormap('gray');
axis equal
xlim([0 size(Second_Map_Hex,2)]);
ylim([0 size(Second_Map_Hex,1)]);
axis off
%% 存Data
for p=1:length(Class_index_array)
    dlmwrite([Training_Data_Storage_Folder_Path,sprintf('Number_%d_Set_%d_Class_%d_Data_Array.txt',Number_of_Superpixel,set,Class_index_array(p))],Data_Array{p},'delimiter','\t','newline','pc','precision', '%.6f');
end

%% show含邊界之First與Second map
First_Map_Hex_with_Boundary=First_Map_Hex;
Second_Map_Hex_with_Boundary=Second_Map_Hex;

First_Map_Hex_with_Boundary(ROI_Map_Boundary>0)=0;
Second_Map_Hex_with_Boundary(ROI_Map_Boundary>0)=0;

First_Map_Hex_Color=repmat(First_Map_Hex_with_Boundary,[1 1 3])/max(First_Map_Hex(:));
Second_Map_Hex_Color=repmat(Second_Map_Hex_with_Boundary,[1 1 3])/max(Second_Map_Hex(:));

Color_index_Set=[1 0 0; 0 1 0]';

for p=1:length(Class_index_array)
    ROI_Map_Boundary_Temp=ROI_Map_Boundary;
    ROI_Map_Boundary_Temp(ROI_Map_Boundary==Class_index_array(p))=1;
    ROI_Map_Boundary_Temp(ROI_Map_Boundary~=Class_index_array(p))=0;
    ROI_Map_Boundary_Temp=repmat(ROI_Map_Boundary_Temp,[1 1 3]);
    ROI_Map_Boundary_Temp(:,:,1)=ROI_Map_Boundary_Temp(:,:,1)*Color_index_Set(1,p);
    ROI_Map_Boundary_Temp(:,:,2)=ROI_Map_Boundary_Temp(:,:,2)*Color_index_Set(2,p);
    ROI_Map_Boundary_Temp(:,:,3)=ROI_Map_Boundary_Temp(:,:,3)*Color_index_Set(3,p);

    First_Map_Hex_Color=First_Map_Hex_Color+ROI_Map_Boundary_Temp;
    Second_Map_Hex_Color=Second_Map_Hex_Color+ROI_Map_Boundary_Temp;
end

subplot(1,2,1)
imagesc(First_Map_Hex_Color);
colormap('gray');
axis equal
xlim([0 size(First_Map_Hex_Color,2)]);
ylim([0 size(First_Map_Hex_Color,1)]);
axis off

subplot(1,2,2)
imagesc(Second_Map_Hex_Color);
colormap('gray');
axis equal
xlim([0 size(Second_Map_Hex_Color,2)]);
ylim([0 size(Second_Map_Hex_Color,1)]);
axis off


%% record values to file

imwrite(First_Map_Hex_Color,[Training_Data_Storage_Folder_Path,sprintf('Number_%d_Set_%d_Image_First_Map_Hex.png',Number_of_Superpixel,set)],'png');
imwrite(Second_Map_Hex_Color,[Training_Data_Storage_Folder_Path,sprintf('Number_%d_Set_%d_Image_Second_Map_Hex.png',Number_of_Superpixel,set)],'png');

imwrite(First_Map_Hex/max(First_Map_Hex(:)),[Training_Data_Storage_Folder_Path,sprintf('Number_%d_Set_%d_Image_First_Map_Hex_no_Frame.png',Number_of_Superpixel,set)],'png');
imwrite(Second_Map_Hex/max(Second_Map_Hex(:)),[Training_Data_Storage_Folder_Path,sprintf('Number_%d_Set_%d_Image_Second_Map_Hex_no_Frame.png',Number_of_Superpixel,set)],'png');


%% Try to do k-mean
if If_Analysis == 1

Data=[reshape(Ratio_Map_Hex,[],1)];% reshape(Std_map_Hex,[],1)];% reshape(X_grid_Hex,[],1) reshape(Y_grid_Hex,[],1)];

idx=kmeans(Data,2);

idx_map=reshape(idx,size(Ratio_Map_Hex,1),size(Ratio_Map_Hex,2));

imagesc(idx_map);
end
%% 示意圖
if If_schematic == 1
     QQ=136;
    % 
    % subplot (2,1,1)
    % 
    % imagesc(B_scan_Reduced_Adjusted);
    %     colormap('gray');
    %     caxis([0 1]);
    %     axis equal
    %     xlim([0 size(B_scan_Reduced_Adjusted,2)]);
    %     ylim([0 size(B_scan_Reduced_Adjusted,1)]);
    %     axis off
    %     
    %     
    % subplot (2,1,2)

    Forced_max=1;
    Forced_min=0;
    B_scan_Reduced_Adjusted_Note=repmat(B_scan_Reduced_Adjusted,[1 1 3]);
    B_scan_Reduced_Adjusted_Note=(B_scan_Reduced_Adjusted_Note-Forced_min)/(Forced_max-Forced_min);
    B_scan_Reduced_Adjusted_Note(B_scan_Reduced_Adjusted_Note>1)=1;
    B_scan_Reduced_Adjusted_Note(B_scan_Reduced_Adjusted_Note<0)=0;

    for p=1:size(B_scan_Reduced_Adjusted_Note,2)
        B_scan_Reduced_Adjusted_Note((max_index_map_Closed_Reduced(p,QQ)+1):(max_index_map_Closed_Reduced(p,QQ)+First_index),p,1)=1;
        B_scan_Reduced_Adjusted_Note((max_index_map_Closed_Reduced(p,QQ)+1):(max_index_map_Closed_Reduced(p,QQ)+First_index),p,2)=0;
        B_scan_Reduced_Adjusted_Note((max_index_map_Closed_Reduced(p,QQ)+1):(max_index_map_Closed_Reduced(p,QQ)+First_index),p,3)=0;

        B_scan_Reduced_Adjusted_Note((max_index_map_Closed_Reduced(p,QQ)+First_index+1):(max_index_map_Closed_Reduced(p,QQ)+First_index+Second_index),p,1)=0;
        B_scan_Reduced_Adjusted_Note((max_index_map_Closed_Reduced(p,QQ)+First_index+1):(max_index_map_Closed_Reduced(p,QQ)+First_index+Second_index),p,2)=1;
        B_scan_Reduced_Adjusted_Note((max_index_map_Closed_Reduced(p,QQ)+First_index+1):(max_index_map_Closed_Reduced(p,QQ)+First_index+Second_index),p,3)=0;
    end

% 
% image(B_scan_Reduced_Adjusted_Note);
%     axis equal
%     xlim([0 size(B_scan_Reduced_Adjusted,2)]);
%     ylim([0 size(B_scan_Reduced_Adjusted,1)]);
%     axis off
end

%% Running window removing outliner (outliner = point away from local mean by N std
% if If_std_map == 1
%     Sample_Size=10;  %size of the matrix, must be odd number
% 
%     %Weight_mask=ones(Sample_Size);
%     Weight_mask=fspecial('gaussian', Sample_Size,Sample_Size/4);    %gaussian mask
%     % imagesc(Weight_mask);
%     % colormap('gray');
%     % axis equal
%     % xlim([0 size(Weight_mask,2)]);
%     % ylim([0 size(Weight_mask,1)]);
%     % axis off
% 
%     % std
%     Std_map=zeros(size(Height_map,1),size(Height_map,2));
%     for p=1:size(Height_map,1)
%         for q=1:size(Height_map,2)
%             Local_Sample=Height_map(p:min(size(Height_map,1),(p+Sample_Size-1)),q:min(size(Height_map,2),(q+Sample_Size-1)));
%             Mean=mean(Local_Sample(:));
%             Weight_mask_Local=Weight_mask(1:min(Sample_Size,Sample_Size-((p+Sample_Size-1)-size(Height_map,1))),1:min(Sample_Size,Sample_Size-((q+Sample_Size-1)-size(Height_map,2))));
%             Weight_mask_Local=Weight_mask_Local/sum(Weight_mask_Local(:));
%             Std_map(p,q)=(sum(sum(((Local_Sample-Mean).^2).*Weight_mask(1:min(Sample_Size,Sample_Size-((p+Sample_Size-1)-size(Height_map,1))),1:min(Sample_Size,Sample_Size-((q+Sample_Size-1)-size(Height_map,2))))))).^0.5;
% 
%         end
%     end
% Std_map=Std_map(ROI(3):ROI(4),ROI(1):ROI(2));
% % %%
% % imagesc(Std_map);
% % colormap('gray');
% % axis equal
% % xlim([0 size(Std_map,2)]);
% % ylim([0 size(Std_map,1)]);
% % axis off
% end
