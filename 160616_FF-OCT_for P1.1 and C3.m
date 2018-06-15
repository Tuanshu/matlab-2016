clear all
tic
% Data reading & stacking
%cd('F:\WideScan_20160615_MeatC4_B3_S60_R50_4X_it1862\');

TH=5;
Offset_increment=0; %for different block
Min_Min_increment=40;   %for the 1st FOV, hysteresis issue
Initial_state=3;
correction_set=0;   %5很爛, 7抓錯example
%If_update_DF=0;
%If_update_FF=0;
If_universal_FF=0;
If_universal_DF=0;
If_single_frame=0;
If_patial=0;
If_outputraw=0;
If_autotilt=1;
If_Zmap_EveryPosition=1;
frame_average_factor=20;    %averaged to stack
stack_sampling_spacing=20;
total_stack_number=10;%10;%20;%20;%40;


if correction_set == 0  %15-21005
    System=2;       %1 for P1.1, 2 for B1
    if_glass_interface_searching=1;
    folder_path_without_index='F:\WideScan_20160615_MeatC4_B3_S60_R50_4X_it1862\';


    X_incre=2.5;%5.714286;%-4.833016;%-45.7000;%19.0844
    Y_incre=13.357143;%7.400571;%-6.6429;%-33.4857
    
    cmin_2=0;   %B1:0.025   P1.1_vib:0   C3_not blacken: 0
    cmax_2=20;%30;  %without Dark frame 0.25;     %B1:0.55   P1.1_vib:0.4  C3_not blacken: 0.4

    Dark_frame_Ratio=1.2;
    
    X_mosaic_starting_index_min=6;%12;
    Y_mosaic_starting_index_min=16;
    X_mosaic_number_max=6;%6;%6;%6;%12;%12;%3;%6;%3;%-9;%6;%3;%12;%12;
    Y_mosaic_number_max=8;%16;%8;%16;%8;%16;%4;%8;%16;%4;%-12;%8;%4;%16;%16;   這4個是用來定義glass interface matrix size的, 不然相同斜度, 但不同FOV, 影像深度會不同, 很不方便

    X_mosaic_starting_index=6;%12;%+3;%+2;%+6;%+6;%+6+3;%+6;%6;%13;%+6;%6;%1;  %start from 0
    Y_mosaic_starting_index=16;%+4;%+2;%+8+4;%+6;%+8;%+8+4;%12;%8;%17;%+8;%8;%1;

    X_mosaic_number=6;%6;%6;%6;%12;%12;%3;%6;%3;%-9;%6;%3;%12;%12;
    Y_mosaic_number=8;%16;%8;%16;%8;%16;%4;%8;%16;%4;%-12;%8;%4;%16;%16;
    
    start_index_offset_from_max=0;%+8*8;%300;%+9*4;%16*4;
    
    angle=185;

    X_start_array=[735 1105 1045 895 785 1135 1735 1860 2335];
    Y_start_array=[1615 635 1810 2610 2965 2965 2115 540 1290];
    
    
    
end

cd(folder_path_without_index);


frame_width=648;
frame_height=488;

X_overlapping=25;
Y_overlapping=21;




frame_width_eff=frame_width-X_overlapping;
frame_height_eff=frame_height-Y_overlapping;


num_of_frame_per_division=16;


%% to find the total number of frame (roughly, floor)
if If_single_frame == 0
    if (X_mosaic_starting_index_min<10)&&(Y_mosaic_starting_index_min<10)
        folder_path_min_min=sprintf('%s_% d_% d\\',folder_path_without_index,Y_mosaic_starting_index_min,X_mosaic_starting_index_min);
    elseif (X_mosaic_starting_index_min>9)&&(Y_mosaic_starting_index_min<10)
        folder_path_min_min=sprintf('%s_ %d_%d\\',folder_path_without_index,Y_mosaic_starting_index_min,X_mosaic_starting_index_min);
    elseif (Y_mosaic_starting_index_min>9)&&(X_mosaic_starting_index_min<10)
        folder_path_min_min=sprintf('%s_%d_ %d\\',folder_path_without_index,Y_mosaic_starting_index_min,X_mosaic_starting_index_min);
    else
        folder_path_min_min=sprintf('%s_%d_%d\\',folder_path_without_index,Y_mosaic_starting_index_min,X_mosaic_starting_index_min);
    end
elseif If_single_frame == 1
    folder_path_min_min=folder_path_without_index;
end
    
    

file_list=dir(folder_path_min_min);
Division_array=[];
for p=1:length(file_list)
    if ~isempty(str2double(file_list(p).name)) && ~isnan(str2double(file_list(p).name))
        Division_array(length(Division_array)+1)=str2double(file_list(p).name);
    end
end
[value index]=max(Division_array);
Division_array((index-1):index)=[];  %刪掉最後一個, 16/5/11 delete another one
max_frame_index=(max(Division_array)+1)*num_of_frame_per_division;  
%%

total_FOV_number=X_mosaic_number*Y_mosaic_number;
k=0;

correction_A=ones(frame_width,frame_height);
Overlapping_Mask=zeros(frame_width,frame_height);
%% revised 2016/3/2 A (overlapping)
for tt=1:X_overlapping
    if X_mosaic_number_max>1
        correction_A(tt,:)=correction_A(tt,:)*((tt-1)/((X_overlapping-1)));
        correction_A(frame_width-tt+1,:)=correction_A(frame_width-tt+1,:)*((tt-1)/((X_overlapping-1))); 
    end 
    %correction_A(tt,:)=correction_A(tt,:)*0.5;
    %correction_A(frame_width-tt+1,:)=correction_A(frame_width-tt+1,:)*0.5;   
    Overlapping_Mask(tt,:)=Overlapping_Mask(tt,:)+1;
    Overlapping_Mask(frame_width-tt+1,:)=Overlapping_Mask(frame_width-tt+1,:)+1;

end
for tt=1:Y_overlapping
    if Y_mosaic_number_max>1
        correction_A(:,tt)=correction_A(:,tt)*((tt-1)/(Y_overlapping-1));
        correction_A(:,frame_height-tt+1)=correction_A(:,frame_height-tt+1)*((tt-1)/((Y_overlapping-1))); 
    end
    %correction_A(:,tt)=correction_A(:,tt)*0.5;
    %correction_A(:,frame_height-tt+1)=correction_A(:,frame_height-tt+1)*0.5; 
    Overlapping_Mask(:,tt)=Overlapping_Mask(:,tt)+1;
    Overlapping_Mask(:,frame_height-tt+1)=Overlapping_Mask(:,frame_height-tt+1)+1;
end

Overlapping_Mask_Array=Overlapping_Mask;

Overlapping_Mask_Array=Overlapping_Mask_Array(:);

%% Interblock offset
% First, to generate folder index map
Interblock_offset=zeros(X_mosaic_number,Y_mosaic_number);
for p=1:Y_mosaic_number
    Offset_factor=floor((Y_mosaic_starting_index+Y_mosaic_number-p)/8);     %for P1.1
    Interblock_offset(:,p)=Offset_factor*Offset_increment;
end
Interblock_offset=Interblock_offset-min(Interblock_offset(:));
%% Search for 5 point's z
if If_autotilt

% identify FOV index
tic
%Note: not using the min_min, because there can be hysteresis issue
if If_Zmap_EveryPosition == 0
    Corner_X_array=[X_mosaic_starting_index_min+1 X_mosaic_starting_index_min X_mosaic_starting_index_min+X_mosaic_number_max-1 X_mosaic_starting_index_min+X_mosaic_number_max-1 X_mosaic_starting_index_min+round(X_mosaic_number_max/2)];
    Corner_Y_array=[Y_mosaic_starting_index_min+1 Y_mosaic_starting_index_min+Y_mosaic_number_max-1 Y_mosaic_starting_index_min Y_mosaic_starting_index_min+Y_mosaic_number_max-1 Y_mosaic_starting_index_min+round(Y_mosaic_number_max/2)];
elseif If_Zmap_EveryPosition == 1
    for p=1:X_mosaic_number_max
        for q=1:Y_mosaic_number_max
            Corner_X_array((p-1)*Y_mosaic_number_max+q)=X_mosaic_starting_index_min+p-1;
            Corner_Y_array((p-1)*Y_mosaic_number_max+q)=Y_mosaic_starting_index_min+q-1;
        end
    end
end
frame_mean=zeros(num_of_frame_per_division*length(Division_array),length(Corner_X_array));
if (If_Zmap_EveryPosition == 1)  && (exist([sprintf('%s_correction_file\\',folder_path_without_index),'Measured_glass_interface_Index_map.txt'], 'file') == 2)

    Measured_glass_interface_Index_map=dlmread([sprintf('%s_correction_file\\',folder_path_without_index),'Measured_glass_interface_Index_map.txt']);

    for p=1:X_mosaic_number_max
        for q=1:Y_mosaic_number_max
            Z_max_index((p-1)*Y_mosaic_number_max+q)=Measured_glass_interface_Index_map(p,q);
        end
    end
    Z_max_index=Z_max_index';
else
    for p=1:length(Corner_X_array)

        if If_single_frame == 0

            if (Corner_X_array(p)<10)&&(Corner_Y_array(p)<10)
                folder_path=sprintf('%s_% d_% d\\',folder_path_without_index,Corner_Y_array(p),Corner_X_array(p));
            elseif (Corner_X_array(p)>9)&&(Corner_Y_array(p)<10)
                folder_path=sprintf('%s_ %d_%d\\',folder_path_without_index,Corner_Y_array(p),Corner_X_array(p));
            elseif (Corner_Y_array(p)>9)&&(Corner_X_array(p)<10)
                folder_path=sprintf('%s_%d_ %d\\',folder_path_without_index,Corner_Y_array(p),Corner_X_array(p));
            else
                folder_path=sprintf('%s_%d_%d\\',folder_path_without_index,Corner_Y_array(p),Corner_X_array(p));
            end
        elseif If_single_frame == 1
            folder_path=folder_path_without_index;
        end

        for NN=1:length(Division_array)
            file_path=[folder_path sprintf('%08d',Division_array(NN))];
            fin=fopen(file_path);
            A=fread(fin,[frame_width,frame_height*num_of_frame_per_division],'float32','b');
            A(isnan(A))=0;
            if fin ==-1
                k=k+1;
            else

                for q=1:num_of_frame_per_division
                    frame_mean((NN-1)*num_of_frame_per_division+q,p)=mean(mean(A(:,(frame_height*(q-1)+1):frame_height*q)));
                end


            end
            fclose('all');
        end
       disp(p); 
    end
    toc
    frame_mean=frame_mean-min(frame_mean(:));
    Depth_micron=0.2*[1:length(frame_mean)]';

    Z_max_index=ones([length(Corner_X_array) 1]);
    Z_max_value=ones([length(Corner_X_array) 1]);
    % calculate max
    % for p=1:5
    %     [Z_max_value(p) Z_max_index(p)]=max(frame_mean(:,p));
    % end
    % Z_max_Pos=Z_max_index*0.2;
    % %

    %find last peak
    for p=1:length(Corner_X_array)
        %[Z_max_value_array Z_max_index_array]=findpeaks(frame_mean(:,p),'MINPEAKDISTANCE',1,'MINPEAKHEIGHT',10);
        [Z_max_value_now Z_max_index_now]=max(frame_mean(:,p));

        if isempty(Z_max_index_now)~=1
            Z_max_index(p)=Z_max_index_now;
            Z_max_value(p)=Z_max_value_now;
        else
            Z_max_index(p)=1;
            Z_max_value(p)=0;
        end

    end


    plot(Depth_micron,frame_mean);
    xlabel('Depth(\mum)');
    ylabel('Intensity (a.u.)');

    hold on


    plot(Depth_micron(Z_max_index),Z_max_value,'r.');
    xlabel('Depth(\mum)');
    ylabel('Intensity (a.u.)');

    legend('Corner 1','Corner 2','Corner 3','Corner 4','Corner 5','Location','BestOutside');

    hold off
end
%% Calc Xinc and Yinc
    X=[Corner_X_array;Corner_Y_array;ones([1 length(Corner_X_array)])]';
    Y=[Z_max_index*0.2];%-min(Z_max_Pos);

    beta_matrix=zeros(size(X,2),size(X,1));
    error_array=zeros(size(X,1),1);
    for p=1:size(X,1)
        X_used=X;
        Y_used=Y;
        X_used(p,:)=[];
        Y_used(p)=[];
        beta_matrix(:,p)=X_used\Y_used;
        error_array(p)=sum((Y_used-X_used*beta_matrix(:,p)).^2);
        disp(p);
    end
    [value index]=min(error_array);

    beta=beta_matrix(:,index);
    if System == 1
        X_incre=beta(1)*(-5);
        Y_incre=beta(2)*(-5);
    elseif System == 2
        X_incre=beta(1)*(5);
        Y_incre=beta(2)*(5);
    end
    disp(index);
    disp(sprintf('Xinc=%3f\n',X_incre))
    disp(sprintf('Yinc=%3f\n',Y_incre))
    disp(sprintf('Total=%3f micron',(abs(X_incre)*24+abs(Y_incre)*32)/5))
end

%% Reconstruct Z map in case of complete Z map
if If_Zmap_EveryPosition
    if exist([sprintf('%s_correction_file\\',folder_path_without_index),'Measured_glass_interface_Index_map.txt'], 'file') == 2
        Measured_glass_interface_Index_map=dlmread([sprintf('%s_correction_file\\',folder_path_without_index),'Measured_glass_interface_Index_map.txt']);
    else
        Measured_glass_interface_index_map=zeros(X_mosaic_number_max,Y_mosaic_number_max);
        for p=1:X_mosaic_number_max
            for q=1:Y_mosaic_number_max
                Measured_glass_interface_Index_map(p,q)=Z_max_index((p-1)*Y_mosaic_number_max+q);
            end
        end
        Measured_glass_interface_Micron_map=Measured_glass_interface_Index_map*0.2;
        mkdir(sprintf('%s_correction_file\\',folder_path_without_index));
        dlmwrite([sprintf('%s_correction_file\\',folder_path_without_index),'Measured_glass_interface_Index_map.txt'],Measured_glass_interface_Index_map,'delimiter','\t','newline','pc','precision', '%.6f');
    end
    
    
end


%% glass interface defination
%Value_1_1=280;
glass_interface_index_map_set_all=zeros(X_mosaic_number_max,Y_mosaic_number_max);
    
for p=1:size(glass_interface_index_map_set_all,1)
    for q=1:size(glass_interface_index_map_set_all,2)
        glass_interface_index_map_set_all(p,q)=glass_interface_index_map_set_all(p,q)+X_incre*(p-1)+Y_incre*(q-1);
    end
end
glass_interface_index_map_set_all=glass_interface_index_map_set_all-min(glass_interface_index_map_set_all(:));

glass_interface_index_map_set_all(1,1)=glass_interface_index_map_set_all(1,1)+Min_Min_increment;

if System == 1 %P1.1
    glass_interface_index_map_set=glass_interface_index_map_set_all((X_mosaic_number_max-X_mosaic_number+1-X_mosaic_starting_index+X_mosaic_starting_index_min):(size(glass_interface_index_map_set_all,1)-X_mosaic_starting_index+X_mosaic_starting_index_min),(Y_mosaic_number_max-Y_mosaic_number+1-Y_mosaic_starting_index+Y_mosaic_starting_index_min):(size(glass_interface_index_map_set_all,2)-Y_mosaic_starting_index+Y_mosaic_starting_index_min))+Interblock_offset;   %這一行很重要, for P1.1
elseif System == 2 %B1 %先用跟P1.1一樣的, 有問題再說
    glass_interface_index_map_set=glass_interface_index_map_set_all((X_mosaic_number_max-X_mosaic_number+1-X_mosaic_starting_index+X_mosaic_starting_index_min):(size(glass_interface_index_map_set_all,1)-X_mosaic_starting_index+X_mosaic_starting_index_min),(Y_mosaic_number_max-Y_mosaic_number+1-Y_mosaic_starting_index+Y_mosaic_starting_index_min):(size(glass_interface_index_map_set_all,2)-Y_mosaic_starting_index+Y_mosaic_starting_index_min))+Interblock_offset;   %這一行很重要, for P1.1
end


%%
total_ave_frame=zeros(648,488);
total_ave_frame_after_correction=zeros(648,488);
TH_map=ones(648,488)*TH;
DF_map=ones(648,488)*999;
DF_map_after=ones(648,488)*999;
state=Initial_state;
while 1
    if state == 1
        Dark_frame=0;
        correction_image=correction_A;
    elseif state == 2
        if If_universal_DF == 1
            Dark_frame=dlmread(['D:\MATLAB\FF-OCT_P1.1_Post Data Processing\_correction_file\\','dark_frame.txt']);
        else
            Dark_frame=dlmread([sprintf('%s_correction_file\\',folder_path_without_index),'dark_frame.txt']);
        end
        correction_image=correction_A;
    elseif state == 3
        if If_universal_DF == 1
            Dark_frame=dlmread(['D:\MATLAB\FF-OCT_P1.1_Post Data Processing\_correction_file\\','dark_frame.txt']);
        else
            Dark_frame=dlmread([sprintf('%s_correction_file\\',folder_path_without_index),'dark_frame.txt']);
        end
        if If_universal_FF == 1
            correction_B_Raw_loaded=dlmread(['D:\MATLAB\FF-OCT_P1.1_Post Data Processing\_correction_file\\','correction_file.txt']);
        else
            correction_B_Raw_loaded=dlmread([sprintf('%s_correction_file\\',folder_path_without_index),'correction_file.txt']);
        end
        
        correction_B=1./(correction_B_Raw_loaded/mean(correction_B_Raw_loaded(:)));
        correction_B(isinf(correction_B))=1;
        correction_image=correction_A.*correction_B;

    end
        
    for NNN=1:total_stack_number
        if if_glass_interface_searching==1
            glass_interface_index_map_set_now=glass_interface_index_map_set;
            glass_interface_index_map_set_now(glass_interface_index_map_set_now>(max_frame_index-frame_average_factor-1-start_index_offset_from_max-(NNN-1)*stack_sampling_spacing))=max_frame_index-frame_average_factor-1-start_index_offset_from_max-(NNN-1)*stack_sampling_spacing;
        end

        stiched_image=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping);

        for N=0:(total_FOV_number-1)
            X_number=rem(N,X_mosaic_number)+X_mosaic_starting_index; %0~2
            Y_number=floor(N/X_mosaic_number)+Y_mosaic_starting_index; %0~2
            if System == 1
                X_FOV_number=X_mosaic_number-1-X_number+X_mosaic_starting_index;
                Y_FOV_number=Y_mosaic_number-1-Y_number+Y_mosaic_starting_index;
            elseif System == 2
                X_FOV_number=X_number-X_mosaic_starting_index;
                Y_FOV_number=Y_number-Y_mosaic_starting_index;
            end
            if if_glass_interface_searching==1
                division_starting_index=floor((start_index_offset_from_max+glass_interface_index_map_set_now(X_FOV_number+1,Y_FOV_number+1)+(NNN-1)*stack_sampling_spacing)/num_of_frame_per_division);   %start from zero
                the_starting_frame_index_in_the_first_division=start_index_offset_from_max+glass_interface_index_map_set_now(X_FOV_number+1,Y_FOV_number+1)+(NNN-1)*stack_sampling_spacing-division_starting_index*num_of_frame_per_division+1; %corrected, add 1
                division_end_index=     ceil((start_index_offset_from_max+glass_interface_index_map_set_now(X_FOV_number+1,Y_FOV_number+1)+(NNN-1)*stack_sampling_spacing+frame_average_factor)/num_of_frame_per_division)-1;
            else
                division_starting_index=floor((start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)/num_of_frame_per_division);   %start from zero
                the_starting_frame_index_in_the_first_division=start_index_offset_from_max+(NNN-1)*stack_sampling_spacing-division_starting_index*num_of_frame_per_division+1; %corrected, add 1
                division_end_index=     ceil((start_index_offset_from_max+(NNN-1)*stack_sampling_spacing+frame_average_factor)/num_of_frame_per_division)-1;    %應該是ceil +1 而非floor
            end
            division_number=division_starting_index:division_end_index;
            temp_frame_volume=zeros(frame_width,frame_height,num_of_frame_per_division*length(division_number));
            
            if If_single_frame == 0
                if (X_number<10)&&(Y_number<10)
                    folder_path=sprintf('%s_% d_% d\\',folder_path_without_index,Y_number,X_number);
                elseif (X_number>9)&&(Y_number<10)
                    folder_path=sprintf('%s_ %d_%d\\',folder_path_without_index,Y_number,X_number);
                elseif (Y_number>9)&&(X_number<10)
                    folder_path=sprintf('%s_%d_ %d\\',folder_path_without_index,Y_number,X_number);
                else
                    folder_path=sprintf('%s_%d_%d\\',folder_path_without_index,Y_number,X_number);
                end
            elseif If_single_frame == 1
                folder_path=folder_path_without_index;
            end

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
            end
            Averaged_frame=mean(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),3);
            if state == 1
                Averaged_frame_Temp=Averaged_frame;
                Averaged_frame_Temp(Averaged_frame<TH)=999;
                DF_map=max(TH_map,min(DF_map,Averaged_frame_Temp));
            end

            Averaged_frame=Averaged_frame-Dark_frame*Dark_frame_Ratio;

            Averaged_frame(Averaged_frame<0)=0; 

            if (state == 2) && (max(Averaged_frame(:))~=0) && (max(isnan(Averaged_frame(:)))==0)% && max(isnan(Averaged_frame(:)))~=1
                total_ave_frame=total_ave_frame+Averaged_frame;
            end

            Averaged_frame=Averaged_frame.*correction_image;   %這次試把]normalization放在corre後

            stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))=stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))+flipud(fliplr(Averaged_frame));

            disp(N);

        end
        Hist_max=min(50,max(stiched_image(:)));
        Intensity_Step_Array=0:0.01:Hist_max;
        Stitched_Histogram=hist(stiched_image(:),Intensity_Step_Array);
        Stitched_Histogram=Stitched_Histogram(2:(length(Stitched_Histogram)-1));
        Intensity_Step_Array=Intensity_Step_Array(2:(length(Intensity_Step_Array)-1));
        plot(Intensity_Step_Array,Stitched_Histogram);
        Ratio=sum(Stitched_Histogram(Intensity_Step_Array<=cmax_2))/sum(Stitched_Histogram);
        disp(Ratio);
        stiched_image=(stiched_image-cmin_2)/(cmax_2-cmin_2);  %因為無從做max intensity判斷, 只好用固定值
        %stiched_image(stiched_image<0)=0; 
        %stiched_image(stiched_image>1)=1; 
        %%
        mkdir(sprintf('%s_stiched_image\\',folder_path_without_index));
        if state==1
            imwrite(stiched_image,[sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_offset%.1f micron_X%d_Y%d_original',(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index),'.png']);
        elseif state == 2
            imwrite(stiched_image,[sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_offset%.1f micron_X%d_Y%d_with DF',(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index),'.png']);

        elseif state == 3
            imwrite(stiched_image,[sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_offset%.1f micron_X%d_Y%d_with DF_with corr',(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index),'.png']);
        
            
            if If_outputraw==1
                fout=fopen([sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_offset%.1f micron_X%d_Y%d',(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index)],'w+');
                fwrite(fout,stiched_image,'float32','b');
            end

            if If_patial==1
              %% Rotate   
                %angle=4;
                stiched_image_rotated=imrotate(stiched_image,angle);

                imagesc(stiched_image_rotated);
                colormap('gray');
                caxis([0 1]);
                axis equal
                xlim([0 size(stiched_image_rotated,2)]);
                ylim([0 size(stiched_image_rotated,1)]);

                %% Patial
                %X_start_array=[225 515 295 1465 1340 1515 1965 1990 2365];   %12: 1435   11: 835    10: 1660   9: 2035 8: 1885    7: 2085    6: 2135 5:2310 1: 750 2: 975    3: 2150     4: 2325
                %Y_start_array=[4600 4320 3620 3320 2420 2095 1145 845 220];   %12: 2600   11: 1700   10: 2250   9: 2250 8: 2075    7: 2000    6: 1750 5:1225 1: 550 2: 525     3: 250     4: 700
                X_size=768;
                Y_size=X_size/3*4;
                for p=1:length(X_start_array)
                    if NNN == 1
                        mkdir(sprintf('%s_stiched_image\\%d\\',folder_path_without_index,p));
                    end

                    stiched_image_rotated_patial=stiched_image_rotated(X_start_array(p):(X_start_array(p)+X_size-1),Y_start_array(p):(Y_start_array(p)+Y_size-1));

                    imwrite(stiched_image_rotated_patial,[sprintf('%s_stiched_image\\%d\\',folder_path_without_index,p),sprintf('patial_image_Number%d_D%.1f',p,(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2),'.png']);
                end
            end

        end

    end

    %% Dark-frame processing (Search for supremum of the Dark frame)
    if state ==1

        Times=4;
        Search_Window_Half_Size=10;
        Search_Window_Weight=fspecial('gaussian', Search_Window_Half_Size*2+1,Search_Window_Half_Size/2);    %gaussian mask
        Search_Window_Weight=Search_Window_Weight/max(Search_Window_Weight(:)); %因為是要取max而非平均, 所以filter的最大值應該是1
        Dark_frame_Temp=DF_map;
        Dark_frame_Out=DF_map;
        for w=1:Times
            for p=1:size(Dark_frame_Temp,1)
                for q=1:size(Dark_frame_Temp,2)
                    x_min=max((p-Search_Window_Half_Size+1),1);
                    x_max=min((p+Search_Window_Half_Size),size(Dark_frame_Temp,1));
                    y_min=max((q-Search_Window_Half_Size+1),1);
                    y_max=min((q+Search_Window_Half_Size),size(Dark_frame_Temp,2));
                    Window_Temp=Search_Window_Weight((x_min-p+Search_Window_Half_Size+1):(x_max-(p+Search_Window_Half_Size)+Search_Window_Half_Size*2+1),(y_min-q+Search_Window_Half_Size+1):(y_max-(q+Search_Window_Half_Size)+Search_Window_Half_Size*2+1)).*Dark_frame_Temp(x_min:x_max,y_min:y_max);
                    Dark_frame_Out(p,q)=max(Window_Temp(:));
                end
            end
            Dark_frame_Temp=Dark_frame_Out;
            disp(w);
        end

        subplot(1,2,1)

        imagesc(DF_map);
            colormap('gray');
            %caxis([0 1]);
            axis equal
            xlim([0 size(DF_map,2)]);
            ylim([0 size(DF_map,1)]);

        subplot(1,2,2)

        imagesc(Dark_frame_Out);
            colormap('gray');
            %caxis([0 1]);
            axis equal
            xlim([0 size(Dark_frame_Out,2)]);
            ylim([0 size(Dark_frame_Out,1)]);
        Dark_frame=Dark_frame_Out;


        mkdir(sprintf('%s_correction_file\\',folder_path_without_index));
        dlmwrite([sprintf('%s_correction_file\\',folder_path_without_index),'Dark_frame.txt'],Dark_frame,'delimiter','\t','newline','pc','precision', '%.6f');

    %%
    elseif state == 2

%     imagesc(total_ave_frame);
%     %caxis([0 1]);
%     colormap(gray);
%     axis equal
%     xlim([0 size(total_ave_frame,2)]);
%     ylim([0 size(total_ave_frame,1)]);

    %% try to generate the correction image based on the total_ave_frame (before correction) 這只是用來測試用的, 實際上設定是在後半
    % the gain difference (not sure why it exist, should be corrected by the bobcat
    Gain_Diff_Bound_Pixel=324;
    Gain_ratio=1.005;    %1.02 for cmin=0   1.045 for cmin=5

    Lastline_pixel=488;
    Blur_Size=2;  %size of the matrix, must be odd number
    
    total_ave_frame_gained=total_ave_frame;

    Weight_mask=fspecial('gaussian', Blur_Size,Blur_Size/3);    %gaussian mask


    total_ave_frame_gained(1:Gain_Diff_Bound_Pixel,:)=total_ave_frame(1:Gain_Diff_Bound_Pixel,:)*Gain_ratio;

    total_ave_frame_gained_diff=diff(total_ave_frame_gained,1);
%     imagesc(total_ave_frame_gained);
    %%  這只是用來測試用的, 實際上設定是在後半
    Lastline_ratio=1.14;       %正確值為1.045 for cmin=0 1.075 for cmin=5
    total_ave_frame_gained_lastline=total_ave_frame_gained;

    total_ave_frame_gained_lastline(:,Lastline_pixel)=total_ave_frame_gained_lastline(:,Lastline_pixel)*Lastline_ratio;
%     subplot(1,1,1)
%     imagesc(total_ave_frame_gained_lastline);
%     xlim([480 488]);
    %%
    total_ave_frame_gained_framed=zeros(size(total_ave_frame_gained_lastline,1)+Blur_Size,size(total_ave_frame_gained_lastline,2)+Blur_Size);
    % inpterpol for total_ave_frame_gained_framed
    X_ori=repmat(((1:size(total_ave_frame_gained_lastline,2))-1),[size(total_ave_frame_gained_lastline,1) 1]);
    Y_ori=repmat(((1:size(total_ave_frame_gained_lastline,1))-1)',[1 size(total_ave_frame_gained_lastline,2)]);

    X_new=repmat(((1:size(total_ave_frame_gained_framed,2))-1)/(size(total_ave_frame_gained_framed,2)-1)*(size(total_ave_frame_gained_lastline,2)-1),[size(total_ave_frame_gained_framed,1) 1]);
    Y_new=repmat(((1:size(total_ave_frame_gained_framed,1))-1)'/(size(total_ave_frame_gained_framed,1)-1)*(size(total_ave_frame_gained_lastline,1)-1),[1 size(total_ave_frame_gained_framed,2)]);

    total_ave_frame_gained_framed_bnd=interp2(X_ori,Y_ori,total_ave_frame_gained_lastline,X_new,Y_new,'linear');
    total_ave_frame_gained_framed=total_ave_frame_gained_framed_bnd;

    total_ave_frame_gained_framed((Blur_Size/2+1):(Blur_Size/2)+size(total_ave_frame_gained_lastline,1),(Blur_Size/2+1):(Blur_Size/2)+size(total_ave_frame_gained_lastline,2))=total_ave_frame_gained_lastline;

    %%
    total_ave_frame_blur_framed=filter2(Weight_mask,total_ave_frame_gained_framed,'same');
    total_ave_frame_blur=total_ave_frame_blur_framed((Blur_Size/2+1):(Blur_Size/2)+size(total_ave_frame_gained,1),(Blur_Size/2+1):(Blur_Size/2)+size(total_ave_frame_gained,2));

    %% Apply the gain correction
    total_ave_frame_blur_final=total_ave_frame_blur;
    total_ave_frame_blur_final(1:Gain_Diff_Bound_Pixel,:)=total_ave_frame_blur(1:Gain_Diff_Bound_Pixel,:)/Gain_ratio;

    total_ave_frame_blur_final(:,Lastline_pixel)=total_ave_frame_blur_final(:,Lastline_pixel)/Lastline_ratio;
    

    %%
    mkdir(sprintf('%s_correction_file\\',folder_path_without_index));
    dlmwrite([sprintf('%s_correction_file\\',folder_path_without_index),'correction_file.txt'],total_ave_frame_blur_final,'delimiter','\t','newline','pc','precision', '%.6f');

    end
    if state<3
        state=state+1;
    else
        break;
    end
end


subplot(1,1,1)

imagesc(stiched_image);
    colormap('gray');
    caxis([0 1]);
    axis equal
    xlim([0 size(stiched_image,2)]);
    ylim([0 size(stiched_image,1)]);
toc 















%% try to find the flat field condition (假設平均而言, 在overlapping外的相鄰兩行之intensity得接近)
Disable=1;
if Disable==0
    Band_Width=10;
    stiched_image_color=uint8(repmat(stiched_image/max(stiched_image(:))*255,[1 1 3]));

    Upper=cell(X_mosaic_number,Y_mosaic_number);        %注意:看到的圖轉了180度, 故upper看起來在右邊
    Lower=cell(X_mosaic_number,Y_mosaic_number);
    Left=cell(X_mosaic_number,Y_mosaic_number);
    Right=cell(X_mosaic_number,Y_mosaic_number);

    %lower
    Lower_Edge_Length=frame_width_eff-X_overlapping;
    for p=1:X_mosaic_number
        for q=1:Y_mosaic_number
            Lower_X_Start_Index=X_overlapping+frame_width_eff*(p-1);
            Lower_Y_Start_Index=Y_overlapping+frame_height_eff*(q-1);
            Lower{p,q}=mean(stiched_image(Lower_X_Start_Index:(Lower_X_Start_Index+Lower_Edge_Length-1),Lower_Y_Start_Index:(Lower_Y_Start_Index+Band_Width-1)),2);
            stiched_image_color(Lower_X_Start_Index:(Lower_X_Start_Index+Lower_Edge_Length-1),Lower_Y_Start_Index:(Lower_Y_Start_Index+Band_Width-1),1)=0;
            stiched_image_color(Lower_X_Start_Index:(Lower_X_Start_Index+Lower_Edge_Length-1),Lower_Y_Start_Index:(Lower_Y_Start_Index+Band_Width-1),2)=255;
            stiched_image_color(Lower_X_Start_Index:(Lower_X_Start_Index+Lower_Edge_Length-1),Lower_Y_Start_Index:(Lower_Y_Start_Index+Band_Width-1),3)=0;
        end
    end


    %upper
    Upper_Edge_Length=frame_width_eff-X_overlapping;
    for p=1:X_mosaic_number
        for q=1:Y_mosaic_number
            Upper_X_Start_Index=X_overlapping+frame_width_eff*(p-1);
            Upper_Y_Start_Index=frame_height_eff+frame_height_eff*(q-1);
            Upper{p,q}=mean(stiched_image(Upper_X_Start_Index:(Upper_X_Start_Index+Upper_Edge_Length-1),(Upper_Y_Start_Index-Band_Width+1):Upper_Y_Start_Index),2);
            stiched_image_color(Upper_X_Start_Index:(Upper_X_Start_Index+Upper_Edge_Length-1),(Upper_Y_Start_Index-Band_Width+1):Upper_Y_Start_Index,1)=255;
            stiched_image_color(Upper_X_Start_Index:(Upper_X_Start_Index+Upper_Edge_Length-1),(Upper_Y_Start_Index-Band_Width+1):Upper_Y_Start_Index,2)=0;
            stiched_image_color(Upper_X_Start_Index:(Upper_X_Start_Index+Upper_Edge_Length-1),(Upper_Y_Start_Index-Band_Width+1):Upper_Y_Start_Index,3)=0;
        end
    end


    %left
    Left_Edge_Length=frame_height_eff-Y_overlapping;
    for p=1:X_mosaic_number
        for q=1:Y_mosaic_number
            Left_X_Start_Index=X_overlapping+frame_width_eff*(p-1);
            Left_Y_Start_Index=Y_overlapping+frame_height_eff*(q-1);
            Left{p,q}=mean(stiched_image(Left_X_Start_Index:(Left_X_Start_Index+Band_Width-1),Left_Y_Start_Index:(Left_Y_Start_Index+Left_Edge_Length-1)),1);
            stiched_image_color(Left_X_Start_Index:(Left_X_Start_Index+Band_Width-1),Left_Y_Start_Index:(Left_Y_Start_Index+Left_Edge_Length-1),1)=0;
            stiched_image_color(Left_X_Start_Index:(Left_X_Start_Index+Band_Width-1),Left_Y_Start_Index:(Left_Y_Start_Index+Left_Edge_Length-1),2)=0;
            stiched_image_color(Left_X_Start_Index:(Left_X_Start_Index+Band_Width-1),Left_Y_Start_Index:(Left_Y_Start_Index+Left_Edge_Length-1),3)=255;
        end
    end


    %right
    Right_Edge_Length=frame_height_eff-Y_overlapping;
    for p=1:X_mosaic_number
        for q=1:Y_mosaic_number
            Right_X_Start_Index=frame_width_eff+frame_width_eff*(p-1);
            Right_Y_Start_Index=Y_overlapping+frame_height_eff*(q-1);
            Right{p,q}=mean(stiched_image((Right_X_Start_Index-Band_Width+1):Right_X_Start_Index,Right_Y_Start_Index:(Right_Y_Start_Index+Right_Edge_Length-1)),1);
            stiched_image_color((Right_X_Start_Index-Band_Width+1):Right_X_Start_Index,Right_Y_Start_Index:(Right_Y_Start_Index+Right_Edge_Length-1),1)=255;
            stiched_image_color((Right_X_Start_Index-Band_Width+1):Right_X_Start_Index,Right_Y_Start_Index:(Right_Y_Start_Index+Right_Edge_Length-1),2)=255;
            stiched_image_color((Right_X_Start_Index-Band_Width+1):Right_X_Start_Index,Right_Y_Start_Index:(Right_Y_Start_Index+Right_Edge_Length-1),3)=0;
        end
    end

    imagesc(stiched_image_color);

        axis equal
        xlim([0 size(stiched_image_color,2)]);
        ylim([0 size(stiched_image_color,1)]);

    imwrite(stiched_image_color,[sprintf('%s_stiched_image\\',folder_path_without_index),'stiched_image_color','.png']);
end

%%  Patial Test
disable=1;
if disable ==0
    
    angle=6;
    
    stiched_image_rotated=imrotate(stiched_image,angle);

    imagesc(stiched_image_rotated);
    colormap('gray');
    caxis([0 1]);
    axis equal
    xlim([0 size(stiched_image_rotated,2)]);
    ylim([0 size(stiched_image_rotated,1)]);
    %%
    %% Patial
    X_start_array=2020;%[645 695 845 1695 2470 5070 5070 4560 3875 3050 1940 2020];
    Y_start_array=3460;%[2605 3505 4755 1905 1705 1155 1755 2425 3290 3740 4300 3460] ;

    X_size=768;
    Y_size=round(X_size/3*4);
    stiched_image_rotated_patial=stiched_image_rotated(X_start_array(1):(X_start_array(1)+X_size-1),Y_start_array(1):(Y_start_array(1)+Y_size-1));

    imagesc(stiched_image_rotated_patial);
    colormap('gray');
    caxis([0 1]);
    axis equal
    xlim([0 size(stiched_image_rotated_patial,2)]);
    ylim([0 size(stiched_image_rotated_patial,1)]);

    %%
    imwrite(stiched_image_rotated_patial,[sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('patial_image_Number%d_D%.1f',p,(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2),'.png']);
end

%%
% Test_Dark_Frame=dlmread(['G:\Everday Experiements\WideScan_20160323_175403_15-1258_B3_S25_R60_1x\_correction_file\','correction_file.txt']);
% 
%     imagesc(Test_Dark_Frame);
%     colormap('gray');
%     %caxis([0 10]);
%     axis equal
%     xlim([0 size(Test_Dark_Frame,2)]);
%     ylim([0 size(Test_Dark_Frame,1)]);
%%
fclose all

%%

plot(squeeze(temp_frame_volume(324,244,:)));