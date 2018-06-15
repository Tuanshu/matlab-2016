clear all

tic

Scanning_Paramter=1;        %1 for speed, 2 for it, 0 for thickness
Header_Size=256;

Depth_Sampling_Resolution=0.2;  %micron

First_index_range=[11 11];
Second_index_range=[15 15];

First_Length=First_index_range(2)-First_index_range(1)+1;
Second_Length=Second_index_range(2)-Second_index_range(1)+1;

Number_of_FOV=First_Length*Second_Length;

Colomn=488;
Row=648;

Speed=[1 2 4 8 16];
it_factor=[16];
it_step=75;
it=it_factor*it_step;



DF_Value=ones(((Scanning_Paramter==1)*length(Speed)+(Scanning_Paramter==2)*length(it_factor)),1);
for r=1:((Scanning_Paramter==1)*length(Speed)+(Scanning_Paramter==2)*length(it_factor))
    if Scanning_Paramter == 1
        folder_path=sprintf('I:\\iio Data\\160727_To characterize DF\\%dx_it%d\\',Speed(r),it);
    elseif Scanning_Paramter == 2
        folder_path=sprintf('I:\\iio Data\\160727_To characterize DF\\%dx_it%d\\',Speed,it(r));
    end
    file_list=dir(folder_path);


    file_path=[folder_path sprintf('%03d_%03d.bin',First_index_range(1),Second_index_range(1))];

    fin=fopen(file_path);
    Header=fread(fin,Header_Size,'ulong');


    fin=fopen(file_path);
    fseek(fin, 256*4, 'bof');

    Image_Temp=fread(fin,[Row,inf],'single'); %*Frame
    Frame=size(Image_Temp,2)/Colomn;

    Image_Stack=zeros(Row,Colomn,Frame);
    Image_Array=zeros(Frame,1);
    for p=1:Frame
        Image_Stack(:,:,p)=Image_Temp(:,(1+(p-1)*Colomn):(p*Colomn));
        Image_Array(p)=mean(mean(Image_Temp(:,(1+(p-1)*Colomn):(p*Colomn))));
    end

    %plot(Image_Array);

    %
    % Q=67;
    % Image_Read=Image_Stack(:,:,Q);
    % 
    % 
    % imagesc(Image_Read);
    % colormap(gray);
    % 
    % Min=0;
    % Max=10;
    % 
    % Image_Norm=(Image_Read-Min*1E7)/((Max-Min)*1E7);
    % Image_Norm(Image_Norm>1)=1;
    % 
    % Image_Norm(Image_Norm<0)=0;
    % 
    % 
    % imagesc(Image_Norm);
    % axis equal
    % fclose('all');

    % To add previous Map

    Start_Depth=226;

    End_Depth=326;

    Virtual_Section_Thickness=[20];
    DF_Level_Array=zeros(length(Virtual_Section_Thickness),1);
    for q=1:length(Virtual_Section_Thickness)

        Averaging_Factor=Virtual_Section_Thickness(q)/Depth_Sampling_Resolution;

        Number_of_Averaged_Frame=floor((End_Depth-Start_Depth+1)/Averaging_Factor);

        Image_Stack_ROI=Image_Stack(:,:,Start_Depth:End_Depth);

        Averaged_Image_Stack_ROI=zeros(Row,Colomn,Number_of_Averaged_Frame);

        for p=1:Number_of_Averaged_Frame
            Averaged_Image_Stack_ROI(:,:,p)=mean(Image_Stack_ROI(:,:,(1+(p-1)*Averaging_Factor):(p)*Averaging_Factor),3);
        end
        DF_map=ones(Row,Colomn)*999;
        for p=1:Number_of_Averaged_Frame      %only neccessary when there are "signal", this time, all intensity come from noise, so this is actually not needed

            DF_map=min(DF_map,Averaged_Image_Stack_ROI(:,:,p));

        end
    %imagesc(DF_map);
        DF_Normal=(DF_map-min(DF_map(:)))/(max(DF_map(:))--min(DF_map(:)));
        DF_Normal(DF_Normal>1)=1;
        DF_Normal(DF_Normal<0)=0;
        if Scanning_Paramter == 1
            imwrite(DF_map/150,[folder_path,sprintf('%dx_it%d.png',Speed(r),it)]);
        elseif Scanning_Paramter == 2
            imwrite(DF_map/150,[folder_path,sprintf('%d.png',Speed,it(r))]);
        end

        DF_Level_Array(q)=DF_map(Row/2,Colomn/2);
        %disp(q);
    end
    DF_Value(r)=mean(DF_map(:));
    disp(r);
end
imagesc(DF_map);
std(DF_map(:))

if Scanning_Paramter==1
    plot(Speed,DF_Value,'-ro')
    xlabel('Speed (x)');
    ylabel('DF Level');
    saveas(gcf,['D:\AMO\160727_Dark frame analysis\',sprintf('Speed vs. DF_it%d.png',it)]);

elseif Scanning_Paramter==2
    plot(it,DF_Value,'-ro')
    xlabel('Exposure Time (ms)');
    ylabel('DF Level');
    saveas(gcf,['D:\AMO\160727_Dark frame analysis\',sprintf('It vs. DF_%dx.png',Speed)]);
end
%% SD of SD
%https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
 n=16;    %the number used for 4-point calc
 Sds_coef=(1-2/(n-1)*(gamma(n/2)/gamma((n-1)/2))^2)^0.5
% 
 Sds_unbaised_coef=gamma((n-1)/2)/gamma(n/2)*((n-1)/2-(gamma(n/2)/gamma((n-1)/2))^2)^0.5

% 
% plot(Virtual_Section_Thickness,DF_Level_Array);
% xlabel('Virtual Section Thickness (micron)');
% ylabel('DF Level');