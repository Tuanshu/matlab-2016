clear all
date='0412';
Point=3;
Number_Before=3;
Number_After=3;

Depth_Before=16;
Depth_After=17.6;

Range_Angle=5;  
Delta_Angle=0.5;

If_Angle_Position_Known=1;
If_Search_for_Best_Depth=0;

Known_or_Predicted_Angle=8.95;

Known_or_Predicted_Position=[35 40];

Depth_After_Search_Range=16:0.8:31.2;

ImA=sum(imread(sprintf('I:\\Everday Experiements\\160413_AP related\\2016_%s_AP_P%d_0min_%d\\_stiched_image\\stiched_image_offset%.1f micron_X0_Y0_with DF_with corr.png',date,Point,Number_Before,Depth_Before)),3);
ImA=ImA/max(ImA(:));
ImB=sum(imread(sprintf('I:\\Everday Experiements\\160413_AP related\\2016_%s_AP_P%d_60min_%d\\_stiched_image\\stiched_image_offset%.1f micron_X0_Y0_with DF_with corr.png',date,Point,Number_After,Depth_After)),3);
ImB=ImB/max(ImB(:));

%% Blur
Blur_Window_Size=1;
filter=fspecial('gaussian',Blur_Window_Size,round(Blur_Window_Size/2));
filter=filter/sum(sum(filter));

ImA_Blur=conv2(ImA,filter,'same');

ImB_Blur=conv2(ImB,filter,'same');
%%
subplot(2,1,1)
imagesc(ImA_Blur);
colormap(gray);
axis equal
xlim([0 size(ImA_Blur,2)]);
ylim([0 size(ImA_Blur,1)]);

subplot(2,1,2)
imagesc(ImB_Blur);
colormap(gray);
axis equal
xlim([0 size(ImB_Blur,2)]);
ylim([0 size(ImB_Blur,1)]);

%%  先旋轉, offset就依賴數學的力量
Center_Angle=Known_or_Predicted_Angle;


Angle_Array=Known_or_Predicted_Angle+(-1*Range_Angle:Delta_Angle:Range_Angle);
Window_Size_X=min(size(ImA_Blur,1),size(ImB_Blur,1));
Window_Size_Y=min(size(ImA_Blur,2),size(ImB_Blur,2));

%%
Temp_A=ImA_Blur(round((size(ImA_Blur,1)-Window_Size_X)/2)+(1:Window_Size_X),round((size(ImA_Blur,2)-Window_Size_Y)/2)+(1:Window_Size_Y));
Temp_B=ImB_Blur(round((size(ImB_Blur,1)-Window_Size_X)/2)+(1:Window_Size_X),round((size(ImB_Blur,2)-Window_Size_Y)/2)+(1:Window_Size_Y));

FFT_A=fft2(Temp_A);
FFT_B=fft2(Temp_B);            
V_Temp=ifftshift(ifftshift(ifft2(conj(FFT_A).*FFT_B),1),2);

imagesc(V_Temp);
%% Grid generation
[value index]=max(V_Temp(:));



X_grid=repmat((1:size(V_Temp,1))',1,size(V_Temp,2));
Y_grid=repmat((1:size(V_Temp,2)),size(V_Temp,1),1);
X_grid_array=X_grid(:);
Y_grid_array=Y_grid(:);

%%
V=0;
Temp_A=ImA_Blur(round((size(ImA_Blur,1)-Window_Size_X)/2)+(1:Window_Size_X),round((size(ImA_Blur,2)-Window_Size_Y)/2)+(1:Window_Size_Y));


Zero_Padding_N=1;

Temp_A_ZP=zeros((2*Zero_Padding_N+1)*size(Temp_A,1),(2*Zero_Padding_N+1)*size(Temp_A,2));
Temp_B_ZP=zeros((2*Zero_Padding_N+1)*size(Temp_B,1),(2*Zero_Padding_N+1)*size(Temp_B,2));
Temp_A_ZP((Zero_Padding_N*size(Temp_A,1)+1):((Zero_Padding_N+1)*size(Temp_A,1)),(Zero_Padding_N*size(Temp_A,2)+1):((Zero_Padding_N+1)*size(Temp_A,2)))=Temp_A;

Value_Temp=0;
Best_MSE=99999999;
if If_Angle_Position_Known == 0
    for p=1:length(Angle_Array)
        ImB_Rotate=imrotate(ImB_Blur,Angle_Array(p));
        Temp_B=ImB_Rotate(round((size(ImB_Rotate,1)-Window_Size_X)/2)+(1:Window_Size_X),round((size(ImB_Rotate,2)-Window_Size_Y)/2)+(1:Window_Size_Y));
        Temp_B_ZP((Zero_Padding_N*size(Temp_B,1)+1):((Zero_Padding_N+1)*size(Temp_B,1)),(Zero_Padding_N*size(Temp_B,2)+1):((Zero_Padding_N+1)*size(Temp_B,2)))=Temp_B/mean(Temp_B(:))*mean(Temp_A_ZP(:));
        FFT_A=fft2(Temp_A_ZP);
        FFT_B=fft2(Temp_B_ZP);            
        V_Temp=ifftshift(ifftshift(ifft2(conj(FFT_A).*FFT_B),1),2);
        V=V_Temp((Zero_Padding_N*size(Temp_B,1)+1):((Zero_Padding_N+1)*size(Temp_B,1)),(Zero_Padding_N*size(Temp_B,2)+1):((Zero_Padding_N+1)*size(Temp_B,2)));
        [value index]=max(V(:));
        %MSE=(sum(sum((Temp_A_ZP-Temp_B_ZP).^2)).^0.5)/(sum(Temp_A_ZP(:))*sum(Temp_B_ZP(:)))^0.5;

        if value>Value_Temp%MSE<Best_MSE%
            Angle_record=Angle_Array(p);
            V_X=X_grid_array(index)-size(V,1)/2-1;
            V_Y=Y_grid_array(index)-size(V,2)/2-1;
            Temp_A_record=Temp_A_ZP;
            Temp_B_record=Temp_B_ZP;
            Value_Temp=value;
        end

        disp(p);
    end
else
    V_X=Known_or_Predicted_Position(1);
    V_Y=Known_or_Predicted_Position(2);
    Angle_record=Known_or_Predicted_Angle;
    ImB_Rotate=imrotate(ImB_Blur,Angle_record);

    Temp_B=ImB_Rotate(round((size(ImB_Rotate,1)-Window_Size_X)/2)+(1:Window_Size_X),round((size(ImB_Rotate,2)-Window_Size_Y)/2)+(1:Window_Size_Y));

    Temp_B_ZP((Zero_Padding_N*size(Temp_B,1)+1):((Zero_Padding_N+1)*size(Temp_B,1)),(Zero_Padding_N*size(Temp_B,2)+1):((Zero_Padding_N+1)*size(Temp_B,2)))=Temp_B;
    Temp_A_record=Temp_A_ZP;
    Temp_B_record=Temp_B_ZP;
end
V_X
V_Y
Angle_record

Temp_A_final=Temp_A_record((Zero_Padding_N*size(Temp_B,1)+1):((Zero_Padding_N+1)*size(Temp_B,1)),(Zero_Padding_N*size(Temp_B,2)+1):((Zero_Padding_N+1)*size(Temp_B,2)));
Temp_B_record_2=circshift(Temp_B_record,[-V_X -V_Y]);
Temp_B_final=Temp_B_record_2((Zero_Padding_N*size(Temp_B,1)+1):((Zero_Padding_N+1)*size(Temp_B,1)),(Zero_Padding_N*size(Temp_B,2)+1):((Zero_Padding_N+1)*size(Temp_B,2)));

Temp_A_final(isnan(Temp_A_final))=0;
Temp_B_final(isnan(Temp_B_final))=0;
Temp_B_final=Temp_B_final/mean(Temp_B_record_2(:))*mean(Temp_A_record(:));
Temp_B_final(Temp_B_final>1)=1;
subplot(1,2,1)
imagesc(Temp_A_final);
colormap(gray);
axis equal
xlim([0 size(Temp_A_final,2)]);
ylim([0 size(Temp_A_final,1)]);

subplot(1,2,2)
imagesc(Temp_B_final);
colormap(gray);
axis equal
xlim([0 size(Temp_B_final,2)]);
ylim([0 size(Temp_B_final,1)]);

%% imshowpair
Comp_result=imfuse(Temp_A_final,Temp_B_final);
% imagesc(Comp_result);
%  
% colormap(gray);
% axis equal
% xlim([0 size(Comp_result,2)]);
% ylim([0 size(Comp_result,1)]);
% axis off

name=sprintf('P%d_B%d_A%d',Point,Number_Before,Number_After);

imwrite(Comp_result(:,:,2),sprintf('D:\\160425_Skin Comp\\%s_B.png',name),'png');
imwrite(Comp_result(:,:,3),sprintf('D:\\160425_Skin Comp\\%s_A.png',name),'png');

imwrite(Comp_result,sprintf('D:\\160425_Skin Comp\\%s.png',name),'png');

% %%
% imagesc(Comp_result(:,:,1));
%  
% colormap(gray);
% axis equal
% xlim([0 size(Comp_result,2)]);
% ylim([0 size(Comp_result,1)]);

%%
if If_Search_for_Best_Depth == 1
    Best_MSE=999999999999;
    for p=1:length(Depth_After_Search_Range)
        ImA=sum(imread(sprintf('I:\\Everday Experiements\\160413_AP related\\2016_%s_AP_P%d_0min_%d\\_stiched_image\\stiched_image_offset%.1f micron_X0_Y0_with DF_with corr.png',date,Point,Number_Before,Depth_Before)),3);
        ImA=ImA/max(ImA(:));
        ImB=sum(imread(sprintf('I:\\Everday Experiements\\160413_AP related\\2016_%s_AP_P%d_60min_%d\\_stiched_image\\stiched_image_offset%.1f micron_X0_Y0_with DF_with corr.png',date,Point,Number_After,Depth_After_Search_Range(p))),3);
        ImB=ImB/max(ImB(:));

        %% Blur
        Blur_Window_Size=1;
        filter=fspecial('gaussian',Blur_Window_Size,round(Blur_Window_Size/2));
        filter=filter/sum(sum(filter));

        ImA_Blur=conv2(ImA,filter,'same');

        ImB_Blur=conv2(ImB,filter,'same');

        %%  先旋轉, offset就依賴數學的力量
        Center_Angle=Known_or_Predicted_Angle;
        Range_Angle=1;  
        Delta_Angle=0.05;

        Angle_Array=Known_or_Predicted_Angle+(-1*Range_Angle:Delta_Angle:Range_Angle);
        Window_Size_X=min(size(ImA_Blur,1),size(ImB_Blur,1));
        Window_Size_Y=min(size(ImA_Blur,2),size(ImB_Blur,2));

        %%
        Temp_A=ImA_Blur(round((size(ImA_Blur,1)-Window_Size_X)/2)+(1:Window_Size_X),round((size(ImA_Blur,2)-Window_Size_Y)/2)+(1:Window_Size_Y));
        Temp_B=ImB_Blur(round((size(ImB_Blur,1)-Window_Size_X)/2)+(1:Window_Size_X),round((size(ImB_Blur,2)-Window_Size_Y)/2)+(1:Window_Size_Y));

        Zero_Padding_N=1;

        Temp_A_ZP=zeros((2*Zero_Padding_N+1)*size(Temp_A,1),(2*Zero_Padding_N+1)*size(Temp_A,2));
        Temp_B_ZP=zeros((2*Zero_Padding_N+1)*size(Temp_B,1),(2*Zero_Padding_N+1)*size(Temp_B,2));
        Temp_A_ZP((Zero_Padding_N*size(Temp_A,1)+1):((Zero_Padding_N+1)*size(Temp_A,1)),(Zero_Padding_N*size(Temp_A,2)+1):((Zero_Padding_N+1)*size(Temp_A,2)))=Temp_A;

        ImB_Rotate=imrotate(ImB_Blur,Known_or_Predicted_Angle);
        Temp_B=ImB_Rotate(round((size(ImB_Rotate,1)-Window_Size_X)/2)+(1:Window_Size_X),round((size(ImB_Rotate,2)-Window_Size_Y)/2)+(1:Window_Size_Y));
        Temp_B_ZP((Zero_Padding_N*size(Temp_B,1)+1):((Zero_Padding_N+1)*size(Temp_B,1)),(Zero_Padding_N*size(Temp_B,2)+1):((Zero_Padding_N+1)*size(Temp_B,2)))=Temp_B;

        Temp_A_final=Temp_A_ZP((Zero_Padding_N*size(Temp_B,1)+1):((Zero_Padding_N+1)*size(Temp_B,1)),(Zero_Padding_N*size(Temp_B,2)+1):((Zero_Padding_N+1)*size(Temp_B,2)));
        Temp_B_record_2=circshift(Temp_B_ZP,[-V_X -V_Y]);
        Temp_B_final=Temp_B_record_2((Zero_Padding_N*size(Temp_B,1)+1):((Zero_Padding_N+1)*size(Temp_B,1)),(Zero_Padding_N*size(Temp_B,2)+1):((Zero_Padding_N+1)*size(Temp_B,2)));

        Temp_A_final(isnan(Temp_A_final))=0;
        Temp_B_final(isnan(Temp_B_final))=0;
        Temp_B_final=Temp_B_final/mean(Temp_B_record_2(:))*mean(Temp_A_ZP(:));
        MSE=sum(sum((Temp_A_final-Temp_B_final).^2)).^0.5;
        
        if MSE<Best_MSE
            Record_Depth=Depth_After_Search_Range(p);
            Best_MSE=MSE;
        end
        
    end
    Record_Depth
end
