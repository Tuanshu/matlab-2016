clear all
tic
% Data reading & stacking
folder_path_without_index='G:\Everday Experiements\Wide_20160226_140253_pork\';

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
    if_glass_interface_searching=1;
    start_index_offset_from_max=0;%16*4;
    
    %cmin=8;%12;%11
    %cmax_set=10;%16;%60;%120;%20;
    %G_Ratio=1.5;
    %Gaussian_X_width=120*G_Ratio;
    %Gaussian_Y_width=140*G_Ratio;
    %X_offset=25;
    %Y_offset=-50;
    
    cmin=10;%12;%11
    cmax_set=11;%16;%60;%120;%20;
    G_Ratio=1.5;
    Gaussian_X_width=500*G_Ratio;%350*G_Ratio;
    Gaussian_Y_width=400*G_Ratio;%250*G_Ratio;
    X_offset=50;%160;
    Y_offset=-70;%-100;
    
end





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


X_mosaic_number=12;%12;%3;%6;%3;%-9;%6;%3;%12;%12;
Y_mosaic_number=16;%16;%4;%8;%16;%4;%-12;%8;%4;%16;%16;

%% glass interface searching


frame_average_factor=16;    %averaged to stack

stack_sampling_spacing=8;

total_stack_number=1;


%%


total_FOV_number=X_mosaic_number*Y_mosaic_number;
k=0;


correction_A=ones(frame_width,frame_height);


for tt=1:X_overlapping
    correction_A(tt,:)=correction_A(tt,:)*((tt-1)/((X_overlapping-1)));
    correction_A(frame_width-tt+1,:)=correction_A(frame_width-tt+1,:)*((tt-1)/((X_overlapping-1))); %原本這邊多減了0.5, 移掉15/12/30
end
for tt=1:Y_overlapping
    correction_A(:,tt)=correction_A(:,tt)*((tt-1)/(Y_overlapping-1));
    correction_A(:,frame_height-tt+1)=correction_A(:,frame_height-tt+1)*((tt-1)/((Y_overlapping-1)));
end


imagesc(correction_A);
colormap('gray');
axis equal
xlim([0 size(correction_A,2)]);
ylim([0 size(correction_A,1)]);


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


imagesc(correction_image);
colormap('gray');
axis equal
xlim([0 size(correction_B,2)]);
ylim([0 size(correction_B,1)]);

%% glass interface def

Value_1_1=0;
X_incre=0;%3;
Y_incre=0;%18;
    
glass_interface_index_map_set=zeros(X_mosaic_number,Y_mosaic_number);
    
glass_interface_index_map_set(:)=Value_1_1;
    
for p=1:size(glass_interface_index_map_set,1)
    for q=1:size(glass_interface_index_map_set,2)
        glass_interface_index_map_set(p,q)=glass_interface_index_map_set(p,q)+X_incre*(q-1)+Y_incre*(p-1);
    end
end
imagesc(glass_interface_index_map_set);
colormap('gray');
xlim([1 size(glass_interface_index_map_set,2)]);
ylim([1 size(glass_interface_index_map_set,1)]);
axis equal
       
    
%%

for NNN=1:total_stack_number
    stiched_image=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping);

    for N=0:(total_FOV_number-1)
        X_number=rem(N,X_mosaic_number)+X_mosaic_starting_index; %0~2
        Y_number=floor(N/X_mosaic_number)+Y_mosaic_starting_index; %0~2
        X_FOV_number=X_number-X_mosaic_starting_index;
        Y_FOV_number=Y_number-Y_mosaic_starting_index;
        %Y_FOV_number=X_mosaic_number-1-X_number+X_mosaic_starting_index;
        %X_FOV_number=Y_mosaic_number-1-Y_number+Y_mosaic_starting_index;
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
        Averaged_frame=mean(temp_frame_volume(:,:,the_starting_frame_index_in_the_first_division:(the_starting_frame_index_in_the_first_division+frame_average_factor-1)),3);

        Averaged_frame=Averaged_frame.*correction_image;   

        %Averaged_frame(Averaged_frame>1)=1;
        %NOTE: 2016/2/26 在下面這行加上flip
        stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))=stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))+flipud(fliplr(Averaged_frame));
   
        disp(N);

    end
    %if max(max(stiched_image))<cmax_set
    %    cmax=max(max(stiched_image));
    %else
        %cmax=cmax_set;
    %end
    %Normailzed_image=(stiched_image-cmin)/cmax;
    %Normailzed_image(Normailzed_image>1)=1;
    %Normailzed_image(Normailzed_image<0)=0;
    
       %% 2016/2/26 try this again
    stiched_image=(stiched_image-cmin)/cmax_set;
    stiched_image(stiched_image<0)=0; 
    stiched_image(stiched_image>1)=1; 
    %%
        
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
    
    
 fclose all