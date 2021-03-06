clear all
tic
folder_path_without_index='D:\AMO\160128_Nikon microscope\160128 20_55_54(test2)\';
cd(folder_path_without_index);
    
    cmin=12;%12;%11
    cmax_set=50;%16;%60;%120;%20;
    G_Ratio=1.5;
    Gaussian_X_width=350*G_Ratio;
    Gaussian_Y_width=250*G_Ratio;
    X_offset=95;
    Y_offset=50;
    
frame_width=2048;
frame_height=1536;

X_overlapping=30;
Y_overlapping=30;


frame_width_eff=frame_width-X_overlapping;
frame_height_eff=frame_height-Y_overlapping;

X_mosaic_number=16;%3;%6;%3;%-9;%6;%3;%12;%12;
Y_mosaic_number=8;%4;%8;%16;%4;%-12;%8;%4;%16;%16;


%%
%

total_FOV_number=X_mosaic_number*Y_mosaic_number;


%stiched_volume=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping,num_of_frame_per_division);
% correction image generation
%% Correction related
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

imagesc(correction_image);
colormap('gray');
axis equal
xlim([0 size(correction_B,2)]);
ylim([0 size(correction_B,1)]);

%%
stiched_image_R=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping);
stiched_image_G=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping);
stiched_image_B=zeros(frame_width_eff*X_mosaic_number+X_overlapping,frame_height_eff*Y_mosaic_number+Y_overlapping);
for N=0:(total_FOV_number-1)
    X_FOV_number=floor(N/Y_mosaic_number); %0~2
    if mod(floor(N/Y_mosaic_number),2)==1
        Y_FOV_number=mod(N,Y_mosaic_number); %0~2
    elseif mod(floor(N/Y_mosaic_number),2)==0
        Y_FOV_number=Y_mosaic_number-mod(N,Y_mosaic_number)-1; %0~2
    end
    file_path=[sprintf('%06d.png',N)];
    Raw=imread(file_path,'png');

    stiched_image_R(((X_FOV_number)*frame_width_eff+1):((X_FOV_number)*frame_width_eff+frame_width),((Y_FOV_number)*frame_height_eff+1):((Y_FOV_number)*frame_height_eff+frame_height))=mean(Raw,3)';
    disp(N);
end 
%%
stiched_image_write=stiched_image_R;
mkdir(sprintf('%s_stiched_image\\',folder_path_without_index));
imwrite(stiched_image_write,[sprintf('%s_stiched_image\\',folder_path_without_index),'stiched_image','.png']);
fout=fopen([sprintf('%s_stiched_image\\',folder_path_without_index),'stiched_image'],'w+');
fwrite(fout,stiched_image_R,'float32','b');
    %dlmwrite([sprintf('%s_stiched_image\\',folder_path_without_index),sprintf('stiched_image_offset%.1f micron_X%d_Y%d',(start_index_offset_from_max+(NNN-1)*stack_sampling_spacing)*0.2,X_mosaic_starting_index,Y_mosaic_starting_index),'.txt'],stiched_image_write,'delimiter','\t','newline','pc','precision', '%.6f');
imagesc(stiched_image_write);
    colormap('gray');
    axis equal
    xlim([0 size(stiched_image_R,2)]);
    ylim([0 size(stiched_image_R,1)]);
    