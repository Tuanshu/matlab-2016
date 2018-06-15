clear all
tic
% Data reading & stacking
set=1;
If_std_map=0;
If_schematic=0;
if set == 1     %30 sec, DEJ, good, important
    folder_path='G:\Everday Experiements\151229_MMH\151229_Temp_28th scan\';
    division_number=2:30;
    NNN=136;
    
    Start_Searching_Index=60-2*16;
    End_Searching_Index=200-2*16;
    
    First_index=5;
    Second_index=10;
elseif set == 2     %30 sec, detine only, good result, but less interest
    folder_path='G:\Everday Experiements\151229_MMH\151229_Temp_29th scan\';
    division_number=2:19;
    NNN=136;
    
    Start_Searching_Index=60-2*16;
    End_Searching_Index=200-2*16;
    
    First_index=5;
    Second_index=10;
elseif set == 3 %30 sec, enamel only, not so obvious, near boundary
    folder_path='G:\Everday Experiements\151229_MMH\151229_Temp_30th scan\';
    division_number=2:20;
    NNN=436;
    
    Start_Searching_Index=70-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=5;
    Second_index=10;
    
elseif set == 4 %30 sec, DEJ, only enamel were cut, good
    folder_path='G:\Everday Experiements\151229_MMH\151229_Temp_31th scan\';
    division_number=2:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=10;
    Second_index=10;
elseif set == 5 %30 sec, DEJ, bad
    folder_path='G:\Everday Experiements\151228_MMH\151228_2nd scan_Sampl 3_Temp_Deep100_DEJ\';
    division_number=2:15;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=10;
    Second_index=10;
elseif set == 6 %上一個轉90度, 也看不出來, 但好像有一些膠帶痕跡, 我猜是沒貼好
    folder_path='G:\Everday Experiements\151228_MMH\151228_3rd scan_Sampl 3_Temp_Deep100_DEJ rotate 90\';
    division_number=2:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=10;
    Second_index=10;
elseif set == 7  %在6的旁邊, 好像蠻明顯是沒貼好的
    folder_path='G:\Everday Experiements\151228_MMH\151228_4th scan_Sampl 3_Temp_Deep100_Dentine_rorate 90\';
    division_number=2:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=10;
    Second_index=10;
    
elseif set == 8  %60sec, DEJ, good, 但好像是內層
    folder_path='G:\Everday Experiements\151228_MMH\151228_5th scan_Sampl 2_Temp_Deep100_DEJ\';
    division_number=2:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=10;
    Second_index=10;
    
elseif set == 9  %DEJ, 好像D比E清楚
    folder_path='G:\Everday Experiements\151228_MMH\151228_6th scan_Sampl 2_Temp_Deep100_DEJ\';
    division_number=2:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=10;
    Second_index=10;
    
elseif set == 10  %60sec, repeat previous, better, DEJ
    folder_path='G:\Everday Experiements\151228_MMH\151228_7th scan_Sampl 2_Temp_Deep100_DEJ\';
    division_number=2:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=10;
    Second_index=10;
    
elseif set == 11  %60 sec, DEJ, good
    folder_path='G:\Everday Experiements\151228_MMH\151228_8th scan_Sampl 2_Temp_Deep100_DEJ\';
    division_number=5:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=10;
    Second_index=10;
    
elseif set == 12  %60 sec, DEJ, good
    folder_path='G:\Everday Experiements\151228_MMH\151228_10th scan_Sampl 2_Temp_Deep100_DEJ\';
    division_number=5:20;
    NNN=86;
    Start_Searching_Index=50-2*16;
    End_Searching_Index=180-2*16;
    
    First_index=10;
    Second_index=10;
    
elseif set == 13  %240 sec, 但反而沒那麼清楚, 應該是全畫面皆enamel
    folder_path='G:\Everday Experiements\151228_MMH\151228_11th scan_Sampl 1_Temp_Deep100_0-4\';
    division_number=0:20;
    NNN=86;
    Start_Searching_Index=64;
    End_Searching_Index=180-2*16;
    
    First_index=10;
    Second_index=10;
    
elseif set == 14  %240 sec, DEJ, 但也不清楚
    folder_path='G:\Everday Experiements\151228_MMH\151228_12th scan_Sampl 1_Temp_Deep100_DEJ_4-0\';
    division_number=0:20;
    NNN=86;
    Start_Searching_Index=59;
    End_Searching_Index=180-2*16;
    
    First_index=10;
    Second_index=10;
    
elseif set == 15  %120 sec, 他說有NaN, 那算了
    folder_path='G:\Everday Experiements\151228_MMH\151228_13th scan_Sampl 1_Temp_Deep100_2-0\';
    division_number=0:20;
    NNN=86;
    Start_Searching_Index=64;
    End_Searching_Index=180-2*16;
    
    First_index=10;
    Second_index=10;
    
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
%%

%Adjusted_Volume=(temp_frame_volume-cmin)/cmax_set;  %因為無從做max intensity判斷, 只好用固定值
%Adjusted_Volume(Adjusted_Volume<0)=0; 
%Adjusted_Volume(Adjusted_Volume>1)=1; 


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
% imagesc(Reduced_Volume(:,:,35));
%     caxis([0 50]);
% 
%     colormap('gray');
%     axis equal
    %%
NNN=153;
B_scan_Reduced(:,:)=Reduced_Volume(:,NNN,:);
B_scan_Reduced_Adjusted=(B_scan_Reduced'-cmin)/cmax_set;
% imagesc(B_scan_Reduced_Adjusted);
%     colormap('gray');
%     caxis([0 1]);
%     axis equal
%     xlim([0 size(B_scan_Reduced_Adjusted,2)]);
%     ylim([0 size(B_scan_Reduced_Adjusted,1)]);

% %%
%     
% Q=132;
% B_scan_Temp(:,:)=Reduced_Volume(:,Q,:);
% B_scan_Temp_Adjusted=(B_scan_Temp'-cmin)/cmax_set;
% 
% imagesc(B_scan_Temp_Adjusted); 
% colormap('gray');
% caxis([0 1]);
% axis equal
% xlim([0 size(B_scan_Reduced_Adjusted,2)]);
% ylim([0 size(B_scan_Reduced_Adjusted,1)]);

    %%
    %Q=5;
    %imagesc(Element(:,:,Q));
    
%%


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
% 
% NNN=136;
% B_scan_Opened(:,:)=Opened_Volume(:,NNN,:);
% B_scan_Opened_Adjusted=(B_scan_Opened'-cmin)/cmax_set;
% imagesc(B_scan_Opened_Adjusted);
%     colormap('gray');
%     caxis([0 1]);
%     axis equal
%     xlim([0 size(B_scan_Opened_Adjusted,2)]);
%     ylim([0 size(B_scan_Opened_Adjusted,1)]);

    
    
    
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
%% B-scan Opened along Z
%cmin=8;%12;%11
%cmax_set=80;%16;%60;%120;%20;
% 
% NNN=136;
% B_scan_Opened_Z(:,:)=Opened_Z_Volume(:,NNN,:);
% B_scan_Opened_Z_Adjusted=(B_scan_Opened_Z'-cmin)/cmax_set;
% 
% subplot (3,1,1)
% 
% imagesc(B_scan_Reduced_Adjusted);
%     colormap('gray');
%     caxis([0 0.5]);
%     axis equal
%     xlim([0 size(B_scan_Reduced_Adjusted,2)]);
%     ylim([0 size(B_scan_Reduced_Adjusted,1)]);
%     axis off
% subplot (3,1,2)
% 
% imagesc(B_scan_Opened_Adjusted);
%     colormap('gray');
%     caxis([0 0.5]);
%     axis equal
%     xlim([0 size(B_scan_Opened_Adjusted,2)]);
%     ylim([0 size(B_scan_Opened_Adjusted,1)]);
%     axis off
% subplot (3,1,3)
% 
% imagesc(B_scan_Opened_Z_Adjusted);
%     colormap('gray');
%     caxis([0 0.5]);
%     axis equal
%     xlim([0 size(B_scan_Opened_Z_Adjusted,2)]);
%     ylim([0 size(B_scan_Opened_Z_Adjusted,1)]);
%     axis off

%%
%% Peak searching 
[max_value max_index_map_Reduced]=max(Reduced_Volume(:,:,:),[],3);
[max_value max_index_map_Opened_Reduced]=max(Opened_Volume(:,:,:),[],3);
[max_value max_index_map_Closed_Reduced]=max(Closed_Volume(:,:,:),[],3);
[max_value max_index_map_Opened_Z_Reduced]=max(Opened_Z_Volume(:,:,:),[],3);

max_index_map=max_index_map_Reduced+Start_Searching_Index-1;
max_index_map_Opened=max_index_map_Opened_Reduced+Start_Searching_Index-1;
max_index_map_Closed=max_index_map_Closed_Reduced+Start_Searching_Index-1;
max_index_map_Opened_Z=max_index_map_Opened_Z_Reduced+Start_Searching_Index-1;
    % %%
    % subplot (3,1,1)
    % 
    % imagesc(max_index_map);
    %     colormap('gray');
    %     axis equal
    %     xlim([0 size(max_index_map,2)]);
    %     ylim([0 size(max_index_map,1)]);
    % subplot (3,1,2)
    % 
    % imagesc(max_index_map_Opened);
    %     colormap('gray');
    %     axis equal
    %     xlim([0 size(max_index_map_Opened,2)]);
    %     ylim([0 size(max_index_map_Opened,1)]);
    % subplot (3,1,3)
    % 
    % imagesc(max_index_map_Opened_Z);
    %     colormap('gray');
    %     axis equal
    %     xlim([0 size(max_index_map_Opened_Z,2)]);
    %     ylim([0 size(max_index_map_Opened_Z,1)]);
    % 


    %%
Height_map=max_index_map*0.2;

Height_map_Closed=max_index_map_Closed*0.2;
    % surf(30-Height_map_Closed(:,end:-1:1),'LineWidth',0.01,'EdgeAlpha',0.1);
    % axis equal
    % xlim([0 size(Height_map_Closed,2)]);
    % ylim([0 size(Height_map_Closed,1)]);
    % axis off
    % view(70,10);

    %% Reversed Height map

    %Reversed_Height_map=30-Height_map;
    % imagesc(Reversed_Height_map);
    % axis equal
    % 
    % xlim([0 size(Height_map,2)]);
    % ylim([0 size(Height_map,1)]);
    % axis off

%% Running window removing outliner (outliner = point away from local mean by N std
if If_std_map == 1
    Sample_Size=10;  %size of the matrix, must be odd number

    %Weight_mask=ones(Sample_Size);
    Weight_mask=fspecial('gaussian', Sample_Size,Sample_Size/4);    %gaussian mask
    % imagesc(Weight_mask);
    % colormap('gray');
    % axis equal
    % xlim([0 size(Weight_mask,2)]);
    % ylim([0 size(Weight_mask,1)]);
    % axis off

    % std
    Std_map=zeros(size(Height_map,1)-Sample_Size,size(Height_map,2)-Sample_Size);
    for p=1:(size(Height_map,1)-Sample_Size+1)
        for q=1:(size(Height_map,2)-Sample_Size+1)
            Local_Sample=Height_map(p:(p+Sample_Size-1),q:(q+Sample_Size-1));
            Mean=mean(Local_Sample(:));
            Std_map(p,q)=(sum(sum(((Local_Sample-Mean).^2).*Weight_mask))).^0.5;

        end
    end

% %%
% imagesc(Std_map);
% colormap('gray');
% axis equal
% xlim([0 size(Std_map,2)]);
% ylim([0 size(Std_map,1)]);
% axis off
end
%% Z separation
First_Zone=zeros(size(Reduced_Volume,1),size(Reduced_Volume,2),First_index);
Second_Zone=zeros(size(Reduced_Volume,1),size(Reduced_Volume,2),Second_index);

for p=1:size(Reduced_Volume,1)
    for q=1:size(Reduced_Volume,2)
        if max_index_map_Closed_Reduced(p,q)+First_index+Second_index<size(Reduced_Volume,3)
            First_Zone(p,q,:)=Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+1):(max_index_map_Closed_Reduced(p,q)+First_index));
            Second_Zone(p,q,:)=Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+First_index+1):(max_index_map_Closed_Reduced(p,q)+First_index+Second_index));
        elseif max_index_map_Closed_Reduced(p,q)+First_index<size(Reduced_Volume,3)
            First_Zone(p,q,:)=Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+1):(max_index_map_Closed_Reduced(p,q)+First_index));
            Second_Zone(p,q,1:size(Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+First_index+1):end),3))=Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+First_index+1):end);
        else
            First_Zone(p,q,1:size(Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+1):end),3))=Reduced_Volume(p,q,(max_index_map_Closed_Reduced(p,q)+1):end);
    
        end
    end
end

Ratio_Map=mean(Second_Zone,3)./mean(First_Zone,3);

subplot (2,1,2)


imagesc(Ratio_Map(:,:));
caxis([0 2]);
colormap('gray');
axis equal
xlim([0 size(Second_Zone,2)]);
ylim([0 size(Second_Zone,1)]);
axis off

%%
Rad3Over2 = sqrt(3) / 2;
[X Y] = meshgrid(0:1:41);
n = size(X,1);
X = Rad3Over2 * X;
Y = Y + repmat([0 0.5],[n,n/2]);
[XV YV] = voronoi(X(:),Y(:)); plot(XV,YV,'b-')
axis equal, axis([10 20 10 20]), zoom on

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

%%