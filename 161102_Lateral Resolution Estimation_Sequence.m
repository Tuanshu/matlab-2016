clear all
%%
Folder_Path='D:\AMO\161102_SFR test sequence\4\';
cd(Folder_Path);
Number=1:1000;
Frame_Rate=25;
Scanning_Speed=0.8; %micron/sec
Depth_Sampling_Resolution=Scanning_Speed/Frame_Rate;

for p=1:length(Number)
    File_Name=sprintf('%06d.bmp',Number(p));

    Image=double(imread([Folder_Path File_Name],'bmp'));

%     imagesc(Image);
%     colormap(gray);
    %
    Y=276;
    X=[45:80]';
    Value=Image(X,Y);

%     plot(X,Value);

    Value_DF=Value-min(Value);

    % Flat field corr
    Sampling_Resolution=0.89;   %micron


    X_Fit_Range=[70 80];

    X_index_start=find(X>=X_Fit_Range(1),1,'first');
    X_index_end=find(X<=X_Fit_Range(2),1,'last');

    X_index=(X_index_start:X_index_end)';
    X_Range=X(X_index);
    Value_Range=Value_DF(X_index);

    Fit_Result=fit(X_Range,Value_Range,'poly1');
    Value_Fit=Fit_Result.p1*X+Fit_Result.p2;
%     plot(X_Range,Value_Range,X,Value_Fit);
%     plot(X,Value_DF,X,Value_Fit);

    Value_FF=Value_DF./Value_Fit;
    X_micron=X*Sampling_Resolution;
%     plot(X,Value_FF,'-o');
%     plot(X_micron,Value_FF,'-o');

    %
    X_micron_fine=X_micron(1):(X_micron(2)-X_micron(1))/100:X_micron(end);
    Value_FF_fine=interp1(X_micron,Value_FF,X_micron_fine);
%     plot(X_micron_fine,Value_FF_fine,'-o');

    FWHM(p)=abs(X_micron_fine(find(Value_FF_fine>0.9,1,'first'))-X_micron_fine(find(Value_FF_fine<0.1,1,'last')))
end

Depth=Depth_Sampling_Resolution*[1:length(FWHM)];
figure('Position', [100, 100, 800, 500]);
plot(Depth,FWHM,'LineWidth',1.5)
xlabel('Axial Position (\mum)','fontsize',15);
ylabel('Lateral Resolution (\mum)','fontsize',15);
xlim([Depth(1) Depth(end)])
ylim([1.5 4.5])
set(gca,'fontsize',15)
