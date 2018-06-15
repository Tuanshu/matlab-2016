clear all
%%
folder_path='I:\Everday Experiements\P1.1 Data Backup\2015_0622_1725_Car_M7\';
cd(folder_path);

% Data format related
num_of_frame_per_division=32;
frame_width=648;
frame_height=488;

% Processing related
ave_factor=19;
micron_per_frame=0.2/ave_factor/4;
N=4;
Offset_1=1;
Offset_2=1; 

%% 計算for each ave frame, 所需讀的division數及編號


    
    
%%
file_list=dir(folder_path);
Division_array=[];
for p=1:length(file_list)
    if ~isempty(str2double(file_list(p).name)) && ~isnan(str2double(file_list(p).name))
        Division_array(length(Division_array)+1)=str2double(file_list(p).name);
    end
end
[value index]=max(Division_array);
Division_array((index-1):index)=[];  %刪掉最後一個, 16/5/11 delete another one

%%
ave_frame_length=floor(num_of_frame_per_division*length(Division_array)/ave_factor);
ave_frame_volume=zeros(frame_width,frame_height,ave_frame_length);

for NN=1:length(Division_array)
    file_path=[folder_path sprintf('%08d',Division_array(NN))];
    fin=fopen(file_path);
    A=fread(fin,[frame_width,frame_height*num_of_frame_per_division],'int16','b');
    if fin ==-1
        k=k+1;
        fclose('all');
    else
        for q=1:num_of_frame_per_division
        ave_frame_volume(:,:,(NN-1)*num_of_frame_per_division+q)=A(:,(frame_height*(q-1)+1):frame_height*q);
        end               
        fclose('all');
    end
    disp(NN);
end

%% N-point
Data_ave=zeros([round(length(Data)/ave_factor)-1 1]);
X_ave=zeros([round(length(Data)/ave_factor)-1 1]);


for p=1:length(Data_ave)
    Data_ave(p)=mean(Data((ave_factor*(p-1)+1+Offset_1):(ave_factor*p+Offset_1)));
    X_ave(p)=mean(X((ave_factor*(p-1)+1+Offset_1):(ave_factor*p+Offset_1)));

end

plot(X,Data,X_ave,Data_ave)

Data_Npoint=zeros([round(length(Data_ave)/N)-1 1]);
X_Npoint=zeros([round(length(Data_ave)/N)-1 1]);



for p=1:length(Data_Npoint)
    subarray=(Data_ave((N*(p-1)+1+Offset_2):(N*p+Offset_2)));
    
    Data_Npoint(p)=((N*sum(subarray.^2)-sum(subarray)^2).^0.5)*(2^0.5)/N;
    X_Npoint(p)=mean(X_ave((N*(p-1)+1+Offset_2):(N*p+Offset_2)));

end

plot(X,Data-mean(Data_ave),X_ave,Data_ave-mean(Data_ave),X_Npoint,Data_Npoint);

plot(X,Data_filtered_car,X,Data_filtered_env,X_Npoint,Data_Npoint);


plot((X-1000)*micron_per_frame,Data_filtered_car/max(Data_filtered_car),(X-1000)*micron_per_frame,Data_filtered_env/max(Data_filtered_car),(X_Npoint-1000)*micron_per_frame,Data_Npoint/max(Data_filtered_car));
plot((X-1000)*micron_per_frame,Data_filtered_car/max(Data_filtered_car),(X-1000)*micron_per_frame,Data_filtered_env/max(Data_filtered_car));
xlim([-2.5 2.5]);
ylim([-1 1]);
xlabel('Position (\mum)');
ylabel('Signal (a.u.)');


Data_filtered_env_norm=Data_filtered_env/max(Data_filtered_env);

Data_Npoint_norm=Data_Npoint/max(Data_Npoint);

Data_Npoint_norm_interp=interp1(X_Npoint,Data_Npoint_norm,X);

FWHM_Hilbert=X(find(Data_filtered_env_norm>0.5,1,'last')-find(Data_filtered_env_norm>0.5,1,'first'))*micron_per_frame;
FWHM_Np=X(find(Data_Npoint_norm_interp>0.5,1,'last')-find(Data_Npoint_norm_interp>0.5,1,'first'))*micron_per_frame;

Ratio=FWHM_Np/FWHM_Hilbert;

plot((X-1000)*micron_per_frame,Data_filtered_env_norm,(X-1000)*micron_per_frame,Data_Npoint_norm_interp);
xlim([-2.5 2.5]);
ylim([0 1]);
xlabel('Position (\mum)');
ylabel('Signal (norm)');