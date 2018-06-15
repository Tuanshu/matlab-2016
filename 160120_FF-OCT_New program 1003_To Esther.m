clear all
tic
folder_path_without_index='G:\Everday Experiements\160119\160119_16th scan_S123_R20_4X_B2+B3\';

correction_set=3;
map_method='mean';  %mean std kurtosis skewness diff

if correction_set == 1
    if_glass_interface_searching=1;
    cmin=12;%12;%11
    cmax_set=16;%16;%60;%120;%20;
    start_index_offset_from_max=-30;
    G_Ratio=0.95;
    Gaussian_X_width=400*G_Ratio;
    Gaussian_Y_width=300*G_Ratio;
    X_offset=10;
    Y_offset=30;

elseif correction_set == 2
    if_glass_interface_searching=1;
    cmin=12;%12;%11
    cmax_set=23;%90;%16;%60;%120;%20;
    start_index_offset_from_max=-30;
    G_Ratio=1.5;
    Gaussian_X_width=400*G_Ratio;
    Gaussian_Y_width=200*G_Ratio;
    X_offset=-50;
    Y_offset=0;
elseif correction_set == 3
    if_glass_interface_searching=0;
    start_index_offset_from_max=0;
    cmin=12;%12;%11
    cmax_set=50;%16;%60;%120;%20;
    G_Ratio=1.5;
    Gaussian_X_width=350*G_Ratio;
    Gaussian_Y_width=250*G_Ratio;
    X_offset=95;
    Y_offset=50;
    
end




if_notch=0;

frame_width=648;
frame_height=488;

X_overlapping=25;
Y_overlapping=21;




frame_width_eff=frame_width-X_overlapping;
frame_height_eff=frame_height-Y_overlapping;


num_of_frame_per_division=16;

normalization_factor=50;

X_mosaic_starting_index=6;%6;%13;%+6;%6;%1;  %start from 0
Y_mosaic_starting_index=8;%12;%8;%17;%+8;%8;%1;


X_mosaic_number=6;%3;%6;%3;%-9;%6;%3;%12;%12;
Y_mosaic_number=16;%4;%8;%16;%4;%-12;%8;%4;%16;%16;

%% glass interface searching

search_starting_frame_index=1;   %0.2 micron per frame, start from 1, search for max
search_ending_frame_index=150;

frame_average_factor=16;    %averaged to stack

stack_sampling_spacing=16;

total_stack_number=1;


%%
%

total_FOV_number=X_mosaic_number*Y_mosaic_number;
k=0;


% correction image generation

% correction 1 (edge summation correction)

correction_A=ones(frame_width,frame_height);

% left&right bound

for tt=1:X_overlapping
    correction_A(tt,:)=correction_A(tt,:)*((tt-1)/((X_overlapping-1)));
    correction_A(frame_width-tt+1,:)=correction_A(frame_width-tt+1,:)*((tt-1)/((X_overlapping-1))); %原本這邊多減了0.5, 移掉15/12/30
end
for tt=1:Y_overlapping
    correction_A(:,tt)=correction_A(:,tt)*((tt-1)/(Y_overlapping-1));
    correction_A(:,frame_height-tt+1)=correction_A(:,frame_height-tt+1)*((tt-1)/((Y_overlapping-1)));
end



correction_B_X=ones(frame_width,frame_height);
correction_B_Y=ones(frame_width,frame_height);

for tt=1:frame_height
    correction_B_X(:,tt)=gaussmf((1:frame_width),[Gaussian_X_width frame_width/2+X_offset]);
end
for tt=1:frame_width
    correction_B_Y(tt,:)=gaussmf((1:frame_height),[Gaussian_Y_width frame_height/2+Y_offset]);
end
correction_B=1./(correction_B_X.*correction_B_Y);

correction_image=correction_A.*correction_B;

if if_notch==1
    notch_height=0.85;
    notch_X=222;
    notch_Y=247;
    notch_width=15;
    notch_X_mat=ones(frame_width,frame_height);
    notch_Y_mat=ones(frame_width,frame_height);
    for tt=1:frame_height
        notch_X_mat(:,tt)=notch_height*gaussmf((1:frame_width),[notch_width notch_X]);
    end
    for tt=1:frame_width
        notch_Y_mat(tt,:)=notch_height*gaussmf((1:frame_height),[notch_width notch_Y]);
    end
    notch_image_1=ones(frame_width,frame_height)-(notch_X_mat.*notch_Y_mat);
    notch_height=0.5;
    notch_X=184;
    notch_Y=434;
    notch_width=32;
    notch_X_mat=ones(frame_width,frame_height);
    notch_Y_mat=ones(frame_width,frame_height);
    for tt=1:frame_height
        notch_X_mat(:,tt)=notch_height*gaussmf((1:frame_width),[notch_width notch_X]);
    end
    for tt=1:frame_width
        notch_Y_mat(tt,:)=notch_height*gaussmf((1:frame_height),[notch_width notch_Y]);
    end
    notch_image_2=ones(frame_width,frame_height)-(notch_X_mat.*notch_Y_mat);
    
    correction_image=correction_image./notch_image_1./notch_image_2;
end

%%

for NNN=1:total_stack_number
    stiched_image=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping);

    for N=0:(total_FOV_number-1)
        X_number=rem(N,X_mosaic_number)+X_mosaic_starting_index; %0~2
        Y_number=floor(N/X_mosaic_number)+Y_mosaic_starting_index; %0~2
        X_FOV_number=X_number-X_mosaic_starting_index;
        Y_FOV_number=Y_number-Y_mosaic_starting_index;
        if if_glass_interface_searching==1
            division_starting_index=floor((start_index_offset_from_max+glass_interface_index_map_set(X_FOV_number+1,Y_FOV_number+1)+(NNN-1)*stack_sampling_spacing)/num_of_frame_per_division);   %start from zero
            the_starting_frame_index_in_the_first_division=start_index_offset_from_max+glass_interface_index_map_set(X_FOV_number+1,Y_FOV_number+1)+(NNN-1)*stack_sampling_spacing-division_starting_index*num_of_frame_per_division+1; %corrected, add 1
            division_end_index=     ceil((start_index_offset_from_max+glass_interface_index_map_set(X_FOV_number+1,Y_FOV_number+1)+(NNN-1)*stack_sampling_spacing+frame_average_factor)/num_of_frame_per_division)-1;
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
            Averaged_frame=(Averaged_frame-cmin)/cmax_set;  %因為無從做max intensity判斷, 只好用固定值
            Averaged_frame(Averaged_frame<0)=0; 
            Averaged_frame(Averaged_frame>1)=1; 
            Averaged_frame=Averaged_frame.*correction_image;   %這次試把]normalization放在corre後
        elseif strcmp(map_method,'std')==1
            Averaged_frame=std(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),0,3)./mean(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),3);
            Averaged_frame(Averaged_frame<0)=0; 
            Averaged_frame(Averaged_frame>1)=1; 
            Averaged_frame=Averaged_frame.*correction_A;   %因為是normalized的std, 故不需要做flat-field correction
        elseif strcmp(map_method,'kurtosis')==1 %內建的好像沒有比較快, 就先不改了
            temp_mean=repmat(mean(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),3),[1 1 length(the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1))]);
            Averaged_frame=mean(((temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1))-temp_mean).^4),3)./(std(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),0,3).^4);
            %Averaged_frame=(Averaged_frame-cmin)/cmax_set;  %因為無從做max intensity判斷, 只好用固定值
            Averaged_frame=(Averaged_frame-1.5)/3;   %原本應該要減3 (see kurt定義), 但因為不想要它有負值, so)
            Averaged_frame(Averaged_frame<0)=0; 
            Averaged_frame(Averaged_frame>1)=1; 
            Averaged_frame=Averaged_frame.*correction_A;
        elseif strcmp(map_method,'skewness')==1 %內建的好像沒有比較快, 就先不改了
            temp_mean=repmat(mean(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),3),[1 1 length(the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1))]);
            Averaged_frame=mean(((temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1))-temp_mean).^3),3)./(std(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),0,3).^3);
            %Averaged_frame=(Averaged_frame-cmin)/cmax_set;  %因為無從做max intensity判斷, 只好用固定值
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
     
        stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))=stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))+Averaged_frame;
   
        disp(N);

    end

        
    stiched_image_write=stiched_image;
    stiched_image_write(stiched_image_write>1)=1;
    mkdir(sprintf('%s_stiched_image\\',folder_path_without_index));
    imwrite(stiched_image_write,[sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_offset%.1f micron_X%d_Y%d',(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index),'.png']);
    fout=fopen([sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_offset%.1f micron_X%d_Y%d',(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index)],'w+');
    fwrite(fout,stiched_image,'float32','b');

end
%%
imagesc(stiched_image_write);
    colormap('gray');
    caxis([0 1]);
    axis equal
    xlim([0 size(stiched_image,2)]);
    ylim([0 size(stiched_image,1)]);
    
    
toc
%
 fclose all