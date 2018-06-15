clear all

% Data reading & stacking
Folder_Path='D:\160506_SFR test\';
File_Name='29';


Input_Image=imread([Folder_Path File_Name],'bmp')';

imagesc(Input_Image);
colormap(gray);
axis equal
axis off


ROI_Size=[100 100];

Center_Pixel=[floor(size(Input_Image,1)/2) floor(size(Input_Image,2)/2)];

ROI=[Center_Pixel(1)-floor(ROI_Size(1)/2) Center_Pixel(1)-floor(ROI_Size(1)/2)+ROI_Size(1)-1;Center_Pixel(2)-floor(ROI_Size(2)/2) Center_Pixel(2)-floor(ROI_Size(2)/2)+ROI_Size(2)-1];

ROI_Image=Input_Image(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2));

imagesc(ROI_Image);
colormap(gray);
axis equal
axis off

ROI_Image=double(ROI_Image(size(ROI_Image,1):-1:1,:));

Diff_X=zeros(ROI_Size);
Diff_Y=zeros(ROI_Size);

Diff_X_0=diff(ROI_Image,1,1);
Diff_Y_0=diff(ROI_Image,1,2);

Diff_X(2:end,:)=Diff_X_0;
Diff_Y(:,2:end)=Diff_Y_0;


Diff_All=Diff_X+Diff_Y;
Diff_All=Diff_All/max(Diff_All(:));
Diff_All_BW=im2bw(Diff_All,0.3);

imagesc(Diff_All_BW);
Diff_All_Ske=bwareaopen(bwmorph(Diff_All_BW,'skel',Inf),20);



imagesc(Diff_All_Ske);
colormap(gray);
axis equal
axis off
%%
Offset=1;
Image_Mask=Diff_All_Ske;

Last_One_Array=ROI_Size(2)*ones([ROI_Size(1) 1]);

for p=1:ROI_Size(1)
    if ~isempty(find(Diff_All_Ske(p,:)==1,1,'last'))
        Last_One_Array(p)=find(Diff_All_Ske(p,:)==1,1,'last');
        if (Last_One_Array(p)+Offset)<=ROI_Size(2)
            Last_One_Array(p)=Last_One_Array(p)+Offset;
        end
    end
    Image_Mask(p,1:Last_One_Array(p))=1;
end

imagesc(Image_Mask);


Diff_X_Masked=Image_Mask.*Diff_X;

Diff_X_Masked(1,:)=ROI_Image(1,:);
Recover_ROI_Image=cumsum(Diff_X_Masked,1);


subplot(2,2,1);
imagesc(Diff_X);
colormap(gray);
axis equal
axis off

subplot(2,2,2);
imagesc(Diff_X_Masked);
colormap(gray);
axis equal
axis off

subplot(2,2,3);
imagesc(ROI_Image);
colormap(gray);
axis equal
axis off

subplot(2,2,4);
imagesc(Recover_ROI_Image);
colormap(gray);
axis equal
axis off
%%
imwrite(ROI_Image/max(ROI_Image(:)),[Folder_Path sprintf('%s_ROI.bmp',File_Name)],'bmp');

imwrite(Recover_ROI_Image/max(Recover_ROI_Image(:)),[Folder_Path sprintf('%s_Recover_ROI.bmp',File_Name)],'bmp');

%%
% imagesc(Recover_ROI_Image);
% 
% ROI_Image_with_Boundary = imfuse(ROI_Image,Diff_All_Ske);
% 
% 
% imagesc(ROI_Image_with_Boundary);