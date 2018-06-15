clear all
tic
% Data reading & stacking
folder_path_without_index='G:\Everday Experiements\160217\WideScan_20160217_142818\';

%parts=regexp(folder_path_without_index, '\','split');
%Parent_folder=folder_path_without_index(1:(length(folder_path_without_index)-length(parts{end})));
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

X_mosaic_starting_index=12;%6;%13;%+6;%6;%1;  %start from 0
Y_mosaic_starting_index=16;%12;%8;%17;%+8;%8;%1;


X_mosaic_number=6;%3;%6;%3;%-9;%6;%3;%12;%12;
Y_mosaic_number=8;%4;%8;%16;%4;%-12;%8;%4;%16;%16;

%% glass interface searching


frame_average_factor=130;    %averaged to stack

stack_sampling_spacing=8;

total_stack_number=1;


%%
%

total_FOV_number=X_mosaic_number*Y_mosaic_number;
k=0;


%stiched_volume=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping,num_of_frame_per_division);
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
    %notch_2
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

%dlmwrite('correction_B.txt',correction_B,'delimiter','\t','newline','pc');

%correction_image(:)=1;
imagesc(correction_image);
colormap('gray');
axis equal
xlim([0 size(correction_B,2)]);
ylim([0 size(correction_B,1)]);

%% glass interface searching

Value_1_1=0;
X_incre=0;
Y_incre=0;
    
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
        
    
    
%glass_interface_index_map_set=glass_interface_index_map;
%%

%stiched_images=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping,total_stack_number);
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
          %temp_frame_volume(:,:,(NN-1)*num_of_frame_per_division+q)=A(:,(frame_height*(q-1)+1):frame_height*q).*correction_image;
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
            %Averaged_frame=(Averaged_frame-cmin)/cmax_set;  %因為無從做max intensity判斷, 只好用固定值
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
            %Averaged_frame=(Averaged_frame-cmin)/cmax_set;  %因為無從做max intensity判斷, 只好用固定值
            %Averaged_frame=(Averaged_frame)/2;   %因為是normalized的std, 故不需要做flat-field correction
            Averaged_frame(Averaged_frame<0)=0; 
            Averaged_frame(Averaged_frame>1)=1; 
            Averaged_frame=Averaged_frame.*correction_A;   %因為是normalized的std, 故不需要做flat-field correction
        end
        %Averaged_frame(Averaged_frame>1)=1;
     
        stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))=stiched_image(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))+Averaged_frame;
   
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
    
       %% 以下為新嘗試, 原本是放在corr前的 2015/12/30    2015/1/20 back to the old method
    %stiched_image=(stiched_image-cmin)/cmax_set;  %因為無從做max intensity判斷, 只好用固定值
    %stiched_image(stiched_image<0)=0; 
    %stiched_image(stiched_image>1)=1; 
    %%
        
    stiched_image_write=stiched_image;
    stiched_image_write(stiched_image_write>1)=1;
    mkdir(sprintf('%s_stiched_image\\',folder_path_without_index));
    imwrite(stiched_image_write,[sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_offset%.1f micron_X%d_Y%d',(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index),'.png']);
    fout=fopen([sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_offset%.1f micron_X%d_Y%d',(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index)],'w+');
    fwrite(fout,stiched_image,'float32','b');
    %dlmwrite([sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_offset%.1f micron_X%d_Y%d',(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index),'.txt'],stiched_image_write,'delimiter','\t','newline','pc','precision', '%.6f');

end
%%
imagesc(stiched_image_write);
    colormap('gray');
    caxis([0 1]);
    axis equal
    xlim([0 size(stiched_image,2)]);
    ylim([0 size(stiched_image,1)]);
    
    
    %fin=fopen([sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_offset%.1f micron_X%d_Y%d',(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index)]);
    %loaded_file=fread(fin,[size(stiched_image,1),size(stiched_image,2)],'float32','b');
    
    
    %imagesc(loaded_file);
    %colormap('gray');
    %caxis([0 1]);
    %axis equal
    %xlim([0 size(stiched_image,2)]);
    %ylim([0 size(stiched_image,1)]);
    %imagesc(histeq(stiched_image,[0 1]));
    %colormap('gray');
    %caxis([0 1]);
    %axis equal
    %xlim([0 size(stiched_image,2)]);
    %ylim([0 size(stiched_image,1)]);
    
    
    %fin2=fopen([cd,'\divide\',sprintf('stiched_image_%.1f micron',(starting_frame_index+(NNN-1)*frame_average_factor)*0.2)]);
    %QQQ=fread(fin2,[size(stiched_image,1),size(stiched_image,2)],'float32','b');

    %imagesc(histeq(stiched_image,[0.2:0.001:0.5]));

    %for QQQ=1:4
    %    Normailzed_image=(mean(stiched_volume(:,:,((QQQ-1)*4+1):((QQQ)*4)),3)-cmin)/cmax;
    %    Normailzed_image(Normailzed_image>1)=1;
    %    %stiched_images(:,:,NN)=mean(stiched_volume,3);
    %    imwrite(Normailzed_image,[cd,'\divide\',sprintf('stiched_image_%d.png',(NN-1)*4+QQQ),'.png']);
    %end
    %for QQQ=1:16
    %    Normailzed_image=(mean(stiched_volume(:,:,QQQ),3)-cmin)/cmax;
    %    Normailzed_image(Normailzed_image>1)=1;
        %stiched_images(:,:,NN)=mean(stiched_volume,3);
    %    imwrite(Normailzed_image,[cd,'\divide3\',sprintf('stiched_image_%d.png',(NN-1)*16+QQQ),'.png']);
    %end
    %stiched_volume=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping,num_of_frame_per_division);
toc
%
disable=1;
if disable==0
frame_starting_index=85;
total_frame_number=20;

%stiched_image=mean(stiched_volume(:,:,frame_starting_index:(frame_starting_index+total_frame_number-1)),3);


%
NNN=1;
stiched_image=stiched_images(:,:,NNN);

imagesc(stiched_image);
colormap('gray');
caxis([cmin cmax]);
axis equal
xlim([0 size(stiched_image,2)]);
ylim([0 size(stiched_image,1)]);

Normailzed_image=(stiched_image-cmin)/cmax;
Normailzed_image(Normailzed_image>1)=1;

imagesc(Normailzed_image);
colormap('gray');
caxis([0 1]);
axis equal
xlim([0 size(stiched_image,2)]);
ylim([0 size(stiched_image,1)]);

%% Stacking
%Frame_spacing=5;
%Stacking_starting_frame=1;
%Total_stacked_frame_number=38;

%for qq=1:Total_stacked_frame_number
%    findex=Stacking_starting_frame+(qq-1)*Frame_spacing;
%    stiched_image=mean(stiched_volume(:,:,findex:(findex+total_frame_number-1)),3);
%    Normailzed_image=(stiched_image-cmin)/cmax;
%    Normailzed_image(Normailzed_image>1)=1;
%    imwrite(Normailzed_image,[cd,'\divide\',sprintf('averaged_frame_%d.png',qq),'.png']);
%end

end
 fclose all