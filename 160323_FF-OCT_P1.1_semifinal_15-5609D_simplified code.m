clear all
tic
% Data reading & stacking
folder_path_without_index='G:\Everday Experiements\WideScan_20160316_191803_15-5609D_C2C3_S55_R70\';
TH=5;
correction_set=1;
If_patial=0;
if_record_dark_frame=1;
If_outputraw=0;
If_load_correction_file=1;
If_load_dark_frame=1;
max_frame_index=328;
if correction_set == 1
    map_method='mean';  %mean std kurtosis skewness diff

    if_glass_interface_searching=1;
    start_index_offset_from_max=200-16*(9)+4*12;%+8*8;%300;%+9*4;%16*4;
    
    cmin_2=0;   %B1:0.025   P1.1_vib:0   C3_not blacken: 0
    cmax_2=10;%30;  %without Dark frame 0.25;     %B1:0.55   P1.1_vib:0.4  C3_not blacken: 0.4
    G_Ratio=10;
    Gaussian_X_width=35*G_Ratio;%350*G_Ratio;
    Gaussian_Y_width=25*G_Ratio;%250*G_Ratio;
    Modulation_depth_X=0;%1%1 for max, -1 for min
    Modulation_depth_Y=0;%-1%1 for max, -1 for min
    X_offset=50;%160;
    Y_offset=220;%-100;
    Dark_frame_Ratio=1.4;
    Dark_frame_level=0;%9;%11;

    X_ratio_total_diff=0;%3;
    Y_ratio_total_diff=0;%3;
    
   
    
elseif correction_set == 2
   
    map_method='std';  %mean std kurtosis skewness diff

    if_glass_interface_searching=1;
    
    start_index_offset_from_max=0;%+9*4;%16*4;
    
    cmin_2=0.4;%0.2 for cmin=0
    cmax_2=0.8;%0.85;%16;%60;%120;%20;
    G_Ratio=10;
    Gaussian_X_width=230*G_Ratio;%350*G_Ratio;
    Gaussian_Y_width=100*G_Ratio;%250*G_Ratio;
    Modulation_depth_X=0;%1 for max, -1 for min
    Modulation_depth_Y=0;%1 for max, -1 for min
    X_offset=-100;%160;
    Y_offset=250;%-100;
    Dark_frame_level=7;

    
    X_ratio_total_diff=0;%3;
    Y_ratio_total_diff=0;%3;
    
elseif correction_set == 3  %for 1x
    map_method='diff';  %mean std kurtosis skewness diff

    if_glass_interface_searching=1;
    start_index_offset_from_max=64;%+9*4;%16*4;
    
    cmin_2=0.025;%0.2 for cmin=0
    cmax_2=0.55;%0.85;%16;%60;%120;%20;
    G_Ratio=10;
    Gaussian_X_width=230*G_Ratio;%350*G_Ratio;
    Gaussian_Y_width=100*G_Ratio;%250*G_Ratio;
    Modulation_depth_X=0;%1 for max, -1 for min
    Modulation_depth_Y=0;%1 for max, -1 for min
    
    Dark_frame_level=5.5;
    X_offset=-100;%160;
    Y_offset=250;%-100;
    
    
    X_ratio_total_diff=0;%3;
    Y_ratio_total_diff=0;%3;
    
end



frame_width=648;
frame_height=488;

X_overlapping=25;
Y_overlapping=21;




frame_width_eff=frame_width-X_overlapping;
frame_height_eff=frame_height-Y_overlapping;


num_of_frame_per_division=16;

normalization_factor=50;


X_mosaic_starting_index_min=12;
Y_mosaic_starting_index_min=8;
X_mosaic_number_max=6;%6;%6;%6;%12;%12;%3;%6;%3;%-9;%6;%3;%12;%12;
Y_mosaic_number_max=16;%16;%8;%16;%8;%16;%4;%8;%16;%4;%-12;%8;%4;%16;%16;   這4個是用來定義glass interface matrix size的, 不然相同斜度, 但不同FOV, 影像深度會不同, 很不方便



X_mosaic_starting_index=12;%+3;%+2;%+6;%+6;%+6+3;%+6;%6;%13;%+6;%6;%1;  %start from 0
Y_mosaic_starting_index=8;%+4;%+2;%+8+4;%+6;%+8;%+8+4;%12;%8;%17;%+8;%8;%1;

X_mosaic_number=6;%6;%6;%6;%12;%12;%3;%6;%3;%-9;%6;%3;%12;%12;
Y_mosaic_number=16;%16;%8;%16;%8;%16;%4;%8;%16;%4;%-12;%8;%4;%16;%16;

%% glass interface searching


frame_average_factor=16;    %averaged to stack

stack_sampling_spacing=4;

total_stack_number=1;%20;%40;

%% Load dark frame

if If_load_dark_frame ==1
    Dark_frame_Old=dlmread([sprintf('%s_correction_file\\',folder_path_without_index),'dark_frame.txt']);
else
    Dark_frame=ones(648,488)*Dark_frame_Level;
end
%% Search for supremum of the Dark frame
if If_load_dark_frame ==1

    Times=4;
    Search_Window_Half_Size=10;
    Search_Window_Weight=fspecial('gaussian', Search_Window_Half_Size*2+1,Search_Window_Half_Size/2);    %gaussian mask
    Search_Window_Weight=Search_Window_Weight/max(Search_Window_Weight(:)); %因為是要取max而非平均, 所以filter的最大值應該是1
    Dark_frame_Temp=Dark_frame_Old;
    Dark_frame_Out=Dark_frame_Old;
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

    imagesc(Dark_frame_Old);
        colormap('gray');
        %caxis([0 1]);
        axis equal
        xlim([0 size(Dark_frame_Temp,2)]);
        ylim([0 size(Dark_frame_Temp,1)]);

    subplot(1,2,2)

    imagesc(Dark_frame_Out);
        colormap('gray');
        %caxis([0 1]);
        axis equal
        xlim([0 size(Dark_frame_Out,2)]);
        ylim([0 size(Dark_frame_Out,1)]);
    Dark_frame=Dark_frame_Out;
end
%%

total_FOV_number=X_mosaic_number*Y_mosaic_number;
k=0;

correction_A=ones(frame_width,frame_height);
Overlapping_Mask=zeros(frame_width,frame_height);
%% revised 2016/3/2 A (overlapping)
for tt=1:X_overlapping
    correction_A(tt,:)=correction_A(tt,:)*((tt-1)/((X_overlapping-1)));
    correction_A(frame_width-tt+1,:)=correction_A(frame_width-tt+1,:)*((tt-1)/((X_overlapping-1))); 
    %correction_A(tt,:)=correction_A(tt,:)*0.5;
    %correction_A(frame_width-tt+1,:)=correction_A(frame_width-tt+1,:)*0.5;   
    Overlapping_Mask(tt,:)=Overlapping_Mask(tt,:)+1;
    Overlapping_Mask(frame_width-tt+1,:)=Overlapping_Mask(frame_width-tt+1,:)+1;

end
for tt=1:Y_overlapping
    correction_A(:,tt)=correction_A(:,tt)*((tt-1)/(Y_overlapping-1));
    correction_A(:,frame_height-tt+1)=correction_A(:,frame_height-tt+1)*((tt-1)/((Y_overlapping-1))); 
    %correction_A(:,tt)=correction_A(:,tt)*0.5;
    %correction_A(:,frame_height-tt+1)=correction_A(:,frame_height-tt+1)*0.5; 
    Overlapping_Mask(:,tt)=Overlapping_Mask(:,tt)+1;
    Overlapping_Mask(:,frame_height-tt+1)=Overlapping_Mask(:,frame_height-tt+1)+1;
end

Overlapping_Mask_Array=Overlapping_Mask;

Overlapping_Mask_Array=Overlapping_Mask_Array(:);

imagesc(correction_A);
colormap('gray');
axis equal
xlim([0 size(correction_A,2)]);
ylim([0 size(correction_A,1)]);

%% B (Gaussian)
correction_B_X=ones(frame_width,frame_height);
correction_B_Y=ones(frame_width,frame_height);

for tt=1:frame_height
    if Modulation_depth_X>0
        correction_B_X(:,tt)=(1-Modulation_depth_X)+Modulation_depth_X*gaussmf((1:frame_width),[Gaussian_X_width frame_width/2+X_offset]);
    elseif Modulation_depth_X<=0
        correction_B_X(:,tt)=(1-Modulation_depth_X)+Modulation_depth_X*gaussmf((1:frame_width),[Gaussian_X_width frame_width/2+X_offset]);
    end
end
for tt=1:frame_width
    correction_B_Y(tt,:)=(1-Modulation_depth_Y)+Modulation_depth_Y*gaussmf((1:frame_height),[Gaussian_Y_width frame_height/2+Y_offset]);
end

if If_load_correction_file==0
    correction_B=1./(correction_B_X.*correction_B_Y);
else
    correction_B_Raw_loaded=dlmread([sprintf('%s_correction_file\\',folder_path_without_index),'correction_file.txt']);
    correction_B_Raw_Normalized=correction_B_Raw_loaded/mean(correction_B_Raw_loaded(:));
    correction_B=1./(correction_B_Raw_Normalized.*correction_B_X.*correction_B_Y);
end

%%
correction_image=correction_A.*correction_B;

imagesc(correction_image);
colormap('gray');
axis equal
xlim([0 size(correction_B,2)]);
ylim([0 size(correction_B,1)]);

%% glass interface searching
%Value_1_1=280;
X_incre=-10;% -25;%-14;%3;
Y_incre=8;%   10;%10;%18;
glass_interface_index_map_set_all=zeros(X_mosaic_number_max,Y_mosaic_number_max);
    
glass_interface_index_map_set_all(:)=0;
    
for p=1:size(glass_interface_index_map_set_all,1)
    for q=1:size(glass_interface_index_map_set_all,2)
        glass_interface_index_map_set_all(p,q)=glass_interface_index_map_set_all(p,q)+X_incre*(p-1)+Y_incre*(q-1);
    end
end
glass_interface_index_map_set_all=glass_interface_index_map_set_all-min(glass_interface_index_map_set_all(:));
glass_interface_index_map_set=glass_interface_index_map_set_all((X_mosaic_number_max-X_mosaic_number+1-X_mosaic_starting_index+X_mosaic_starting_index_min):(size(glass_interface_index_map_set_all,1)-X_mosaic_starting_index+X_mosaic_starting_index_min),(Y_mosaic_number_max-Y_mosaic_number+1-Y_mosaic_starting_index+Y_mosaic_starting_index_min):(size(glass_interface_index_map_set_all,2)-Y_mosaic_starting_index+Y_mosaic_starting_index_min));   %這一行很重要, for P1.1

imagesc(glass_interface_index_map_set);
colormap('gray');
xlim([0 size(glass_interface_index_map_set,2)]);
ylim([0 size(glass_interface_index_map_set,1)]);
axis equal
        
   
%%
total_ave_frame=zeros(648,488);
total_ave_frame_after_correction=zeros(648,488);
TH_map=ones(648,488)*TH;
DF_map=ones(648,488)*999;
DF_map_after=ones(648,488)*999;
for NNN=1:total_stack_number
    if if_glass_interface_searching==1
        glass_interface_index_map_set_now=glass_interface_index_map_set;
        glass_interface_index_map_set_now(glass_interface_index_map_set_now>(max_frame_index-frame_average_factor-1-start_index_offset_from_max-(NNN-1)*stack_sampling_spacing))=max_frame_index-frame_average_factor-1-start_index_offset_from_max-(NNN-1)*stack_sampling_spacing;
    end
        stiched_image=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping);
    flag_map=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping);
    for N=0:(total_FOV_number-1)
        X_number=rem(N,X_mosaic_number)+X_mosaic_starting_index; %0~2
        Y_number=floor(N/X_mosaic_number)+Y_mosaic_starting_index; %0~2
        X_FOV_number=X_mosaic_number-1-X_number+X_mosaic_starting_index;
        Y_FOV_number=Y_mosaic_number-1-Y_number+Y_mosaic_starting_index;
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
    
        if (X_number<10)&&(Y_number<10)
            folder_path=sprintf('%s_% d_% d\\',folder_path_without_index,Y_number,X_number);
        elseif (X_number>9)&&(Y_number<10)
            folder_path=sprintf('%s_ %d_%d\\',folder_path_without_index,Y_number,X_number);
        elseif (Y_number>9)&&(X_number<10)
            folder_path=sprintf('%s_%d_ %d\\',folder_path_without_index,Y_number,X_number);
        else
            folder_path=sprintf('%s_%d_%d\\',folder_path_without_index,Y_number,X_number);
        end
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
        end
        if strcmp(map_method,'mean')==1
            Averaged_frame=mean(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),3);
            if if_record_dark_frame == 1
                Averaged_frame_Temp=Averaged_frame;
                Averaged_frame_Temp(Averaged_frame<TH)=999;
                DF_map=max(TH_map,min(DF_map,Averaged_frame_Temp));
            end
            
            Averaged_frame=Averaged_frame-Dark_frame*Dark_frame_Ratio;
            
           if if_record_dark_frame == 1
                Averaged_frame_Temp=Averaged_frame;
                Averaged_frame_Temp(Averaged_frame<TH)=999;
                DF_map_after=max(TH_map,min(DF_map_after,Averaged_frame_Temp));
            end
            
            Averaged_frame(Averaged_frame<0)=0; 
            if (max(Averaged_frame(:))~=0) % && max(isnan(Averaged_frame(:)))~=1
                total_ave_frame=total_ave_frame+Averaged_frame;
                total_ave_frame_after_correction=total_ave_frame_after_correction+Averaged_frame.*correction_B;
            end
            Averaged_frame=Averaged_frame.*correction_image;   %這次試把]normalization放在corre後

        elseif strcmp(map_method,'std')==1
            Averaged_frame=std(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),0,3)./mean(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),3);
            Averaged_frame=(Averaged_frame-cmin_1)/(cmax_1-cmin_1);  %因為無從做max intensity判斷, 只好用固定值
            Averaged_frame(Averaged_frame<0)=0; 
            Averaged_frame(Averaged_frame>1)=1; 
            Averaged_frame=Averaged_frame.*correction_A;   %因為是normalized的std, 故不需要做flat-field correction
        elseif strcmp(map_method,'kurtosis')==1 %內建的好像沒有比較快, 就先不改了
            temp_mean=repmat(mean(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),3),[1 1 length(the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1))]);
            Averaged_frame=mean(((temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1))-temp_mean).^4),3)./(std(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),0,3).^4);
            Averaged_frame=(Averaged_frame-1.5)/3;   %原本應該要減3 (see kurt定義), 但因為不想要它有負值, so)
            Averaged_frame(Averaged_frame<0)=0; 
            Averaged_frame(Averaged_frame>1)=1; 
            Averaged_frame=Averaged_frame.*correction_A;
        elseif strcmp(map_method,'skewness')==1 %內建的好像沒有比較快, 就先不改了
            temp_mean=repmat(mean(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),3),[1 1 length(the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1))]);
            Averaged_frame=mean(((temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1))-temp_mean).^3),3)./(std(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),0,3).^3);
            Averaged_frame=(Averaged_frame)/2;   %因為是normalized的std, 故不需要做flat-field correction
            Averaged_frame(Averaged_frame<0)=0; 
            Averaged_frame(Averaged_frame>1)=1; 
            Averaged_frame=Averaged_frame.*correction_A;   %因為是normalized的std, 故不需要做flat-field correction
        elseif strcmp(map_method,'diff')==1
            temp_mean_all=mean(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),3);
            temp_mean_Q1=mean(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor/4-1)),3);
            temp_mean_Q4=mean(temp_frame_volume(:,:,(the_starting_frame_index_in_the_first_division+frame_average_factor*3/4):(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),3);
            Averaged_frame=(temp_mean_Q1-temp_mean_Q4)./temp_mean_all;
            Averaged_frame(Averaged_frame<0)=0; 
            Averaged_frame(Averaged_frame>1)=1; 
            Averaged_frame=Averaged_frame.*correction_A;   %因為是normalized的std, 故不需要做flat-field correction
        end
        stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))=stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))+flipud(fliplr(Averaged_frame));
        disp(N);

    end
    
    stiched_image=(stiched_image-cmin_2)/(cmax_2-cmin_2);  %因為無從做max intensity判斷, 只好用固定值
    stiched_image(stiched_image<0)=0; 
    stiched_image(stiched_image>1)=1; 
    %%
        
    stiched_image_write=stiched_image;
    stiched_image_write(stiched_image_write>1)=1;
    mkdir(sprintf('%s_stiched_image\\',folder_path_without_index));
    imwrite(stiched_image_write,[sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_offset%.1f micron_X%d_Y%d',(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index),'.png']);
    if If_outputraw==1
        fout=fopen([sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_offset%.1f micron_X%d_Y%d',(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index)],'w+');
        fwrite(fout,stiched_image,'float32','b');
    end
    
    if If_patial==1
      %% Rotate   
        angle=4;
        stiched_image_rotated=imrotate(stiched_image_write,angle);

        imagesc(stiched_image_rotated);
        colormap('gray');
        caxis([0 1]);
        axis equal
        xlim([0 size(stiched_image_rotated,2)]);
        ylim([0 size(stiched_image_rotated,1)]);

        %% Patial
        X_start_array=[225 515 295 1465 1340 1515 1965 1990 2365];   %12: 1435   11: 835    10: 1660   9: 2035 8: 1885    7: 2085    6: 2135 5:2310 1: 750 2: 975    3: 2150     4: 2325
        Y_start_array=[4600 4320 3620 3320 2420 2095 1145 845 220];   %12: 2600   11: 1700   10: 2250   9: 2250 8: 2075    7: 2000    6: 1750 5:1225 1: 550 2: 525     3: 250     4: 700
        X_size=800;
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
%%
subplot(1,1,1)

imagesc(stiched_image_write);
    colormap('gray');
    caxis([0 1]);
    axis equal
    xlim([0 size(stiched_image,2)]);
    ylim([0 size(stiched_image,1)]);
toc 
    
    
%%  Patial Test
disable=1;
if disable ==0
    angle=4;
    stiched_image_rotated=imrotate(stiched_image_write,angle);

    imagesc(stiched_image_rotated);
    colormap('gray');
    caxis([0 1]);
    axis equal
    xlim([0 size(stiched_image_rotated,2)]);
    ylim([0 size(stiched_image_rotated,1)]);
    %%
    %% Patial
    X_start_array=2365;     %11: XX 10: xx 9: 2365    8: 1990   7: 1965 6: 1515    5: 1340    4: 1465   3: 295     2: 515   1:225     [750 975 2150 2325 2310 2135 2085 1885 2035 1660 835 1435];   %12: 1435   11: 835    10: 1660   9: 2035 8: 1885    7: 2085    6: 2135 5:2310 1: 750 2: 975    3: 2150     4: 2325
    Y_start_array=220;      %11: XX 10: XX 9: 220     8: 845   7: 1145 6: 2095    5: 2420    4: 3320 3: 3620    2: 4320 1:4600     [550 525 250 700 1225 1750 2000 2075 2250 2250 1700 2600];   %12: 2600   11: 1700   10: 2250   9: 2250 8: 2075    7: 2000    6: 1750 5:1225 1: 550 2: 525     3: 250     4: 700
    X_size=800;
    Y_size=X_size/3*4;
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

%%
if disable==0

imagesc(total_ave_frame);
%caxis([0 1]);
colormap(gray);
axis equal
xlim([0 size(total_ave_frame,2)]);
ylim([0 size(total_ave_frame,1)]);
total_ave_frame_gained=total_ave_frame;

%% try to generate the correction image based on the total_ave_frame (before correction) 這只是用來測試用的, 實際上設定是在後半
% the gain difference (not sure why it exist, should be corrected by the bobcat

Gain_Diff_Bound_Pixel=324;
Gain_ratio=1.005;    %1.02 for cmin=0   1.045 for cmin=5

Lastline_pixel=488;
Blur_Size=2;  %size of the matrix, must be odd number

Weight_mask=fspecial('gaussian', Blur_Size,Blur_Size/3);    %gaussian mask


total_ave_frame_gained(1:Gain_Diff_Bound_Pixel,:)=total_ave_frame(1:Gain_Diff_Bound_Pixel,:)*Gain_ratio;

total_ave_frame_gained_diff=diff(total_ave_frame_gained,1);
imagesc(total_ave_frame_gained);
%%  這只是用來測試用的, 實際上設定是在後半
Lastline_ratio=1.14;       %正確值為1.045 for cmin=0 1.075 for cmin=5
total_ave_frame_gained_lastline=total_ave_frame_gained;

total_ave_frame_gained_lastline(:,Lastline_pixel)=total_ave_frame_gained_lastline(:,Lastline_pixel)*Lastline_ratio;
subplot(1,1,1)
imagesc(total_ave_frame_gained_lastline);
xlim([480 488]);
%imagesc(total_ave_frame_gained_diff);
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

subplot(3,1,1)
imagesc(total_ave_frame_gained_framed_bnd);

subplot(3,1,2)
imagesc(total_ave_frame_gained_framed);

subplot(3,1,3)
imagesc(total_ave_frame_gained);

%%
total_ave_frame_blur_framed=filter2(Weight_mask,total_ave_frame_gained_framed,'same');
total_ave_frame_blur=total_ave_frame_blur_framed((Blur_Size/2+1):(Blur_Size/2)+size(total_ave_frame_gained,1),(Blur_Size/2+1):(Blur_Size/2)+size(total_ave_frame_gained,2));
subplot(1,1,1)
imagesc(total_ave_frame_blur);
colormap(gray);

%% Apply the gain correction
total_ave_frame_blur_final=total_ave_frame_blur;
total_ave_frame_blur_final(1:Gain_Diff_Bound_Pixel,:)=total_ave_frame_blur(1:Gain_Diff_Bound_Pixel,:)/Gain_ratio;

total_ave_frame_blur_final(:,Lastline_pixel)=total_ave_frame_blur_final(:,Lastline_pixel)/Lastline_ratio;
imagesc(total_ave_frame_blur_final);

colormap(gray);
axis equal
xlim([0 size(total_ave_frame,2)]);
ylim([0 size(total_ave_frame,1)]);


%%
mkdir(sprintf('%s_correction_file\\',folder_path_without_index));
dlmwrite([sprintf('%s_correction_file\\',folder_path_without_index),'correction_file.txt'],total_ave_frame_blur_final,'delimiter','\t','newline','pc','precision', '%.6f');


dlmwrite([sprintf('%s_correction_file\\',folder_path_without_index),'dark_frame.txt'],DF_map,'delimiter','\t','newline','pc','precision', '%.6f');

end
%%
fclose all