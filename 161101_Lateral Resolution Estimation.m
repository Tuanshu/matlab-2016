clear all
%%
Folder_Path='D:\AMO\161101_SFR\';
Number=66;
File_Name=sprintf('Linefield%d.bmp',Number);

Image=double(imread([Folder_Path File_Name],'bmp'));

imagesc(Image);
colormap(gray);
%
Y=276;
X=[45:80]';
Value=Image(X,Y);

plot(X,Value);

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
plot(X_Range,Value_Range,X,Value_Fit);
plot(X,Value_DF,X,Value_Fit);

Value_FF=Value_DF./Value_Fit;
X_micron=X*Sampling_Resolution;
plot(X,Value_FF,'-o');
plot(X_micron,Value_FF,'-o');

%
X_micron_fine=X_micron(1):(X_micron(2)-X_micron(1))/100:X_micron(end);
Value_FF_fine=interp1(X_micron,Value_FF,X_micron_fine);
plot(X_micron_fine,Value_FF_fine,'-o');

FWHM=abs(X_micron_fine(find(Value_FF_fine>0.9,1,'first'))-X_micron_fine(find(Value_FF_fine<0.1,1,'last')))
