clear all

num_of_frame_per_division=16;

Training_Data_Storage_Folder_Path='D:\160421_Teeth_Trianing\\';

folder_path='G:\Everday Experiements\151229_MMH\151229_Temp_28th scan\';


file_list=dir(folder_path);
Division_array=[];
for p=1:length(file_list)
    if ~isempty(str2double(file_list(p).name)) && ~isnan(str2double(file_list(p).name))
        Division_array(length(Division_array)+1)=str2double(file_list(p).name);
    end
end
[value index]=max(Division_array);
Division_array(index)=[];  %刪掉最後一個
max_frame_index=(max(Division_array)+1)*num_of_frame_per_division;  


cmin=12;%12;%11
cmax_set=16;%16;%60;%120;%20;
X_offset=10;
Y_offset=30;

frame_width=648;
frame_height=488;



%%
temp_frame_volume=zeros(frame_width,frame_height,num_of_frame_per_division*length(Division_array));
cd(folder_path);
for NN=1:length(Division_array)
    file_path=[folder_path sprintf('%08d',Division_array(NN))];
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
NNN=401;
subplot(1,1,1)
imagesc(temp_frame_volume(:,:,NNN));
axis equal
axis off