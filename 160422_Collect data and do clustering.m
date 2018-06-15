clear all;

%%

Training_Data_Storage_Folder_Path='D:\160421_Teeth_Trianing\\';

Number_of_Superpixel=200;

%set_array=[1 3 4 8 9 10 11 12];
set_array=[1 3 4 8 9 10]% 11 12];

Class_index_array=[2 3];

Data_Array=cell([length(Class_index_array) 1]);


for p=1:length(set_array)
    Current_Class_index_array=dlmread([Training_Data_Storage_Folder_Path,sprintf('Set_%d_Class_index_array.txt',set_array(p))]);   
    for q=1:length(Current_Class_index_array) 
        Current_Data=dlmread([Training_Data_Storage_Folder_Path,sprintf('Number_%d_Set_%d_Class_%d_Data_Array.txt',Number_of_Superpixel,set_array(p),Current_Class_index_array(q))]);   
        for w=1:length(Class_index_array)
            if Class_index_array(w)==Current_Class_index_array(q)
                Data_Array{w}=[Data_Array{w}; Current_Data];
            end
        end
    end
end 


%% 1D scatter plot
Color_Matrix=[1 0 0; 0 1 0; 0 0 1];
All_Data=[];
for p=1:length(Class_index_array)
    x = Data_Array{p}(:,1);
    y = Data_Array{p}(:,2)./Data_Array{p}(:,1);
    a = 25;
    scatter(x,y,a,Color_Matrix(p,:),'filled');
    hold on
    All_Data=[All_Data;[x y Class_index_array(p)*ones([length(x) 1])]];
end
hold off
xlim([0 5]);
%ylim([15 50]);


%% SVM
SVMStruct=svmtrain(All_Data(:,1:2),All_Data(:,3),'showplot',true);

Result=classify(All_Data(:,1:2),All_Data(:,1:2),All_Data(:,3));

Good=(Result==All_Data(:,3));
Good(Good==0)=[];
Acc=length(Good)/length(Result)