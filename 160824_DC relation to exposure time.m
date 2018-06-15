clear all
close all

Folder_Path='D:\AMO\RD meeting\160824_DC relation to exposure time\';
Current=[0.25 0.4 0.7];


for p=1:length(Current)
    Data=dlmread([Folder_Path sprintf('%1.3gA.txt',Current(p))]);
    ExpTime=Data(:,1);
    MeanValue=Data(:,4);
    Matrix=[MeanValue ones([length(MeanValue) 1])];
    beta=ExpTime\Matrix;
    MeanValue_Fitted=beta(1)*ExpTime+beta(2);
    hold on
    plot(ExpTime,MeanValue,ExpTime,MeanValue_Fitted);
    %plot(ExpTime,MeanValue-MeanValue_Fitted);
    %plot(ExpTime,(MeanValue-beta(2))./ExpTime);
    %plot(ExpTime,MeanValue);
end

hold off
