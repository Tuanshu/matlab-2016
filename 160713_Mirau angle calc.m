clear all


ImA=sum(imread('D:\AMO\160713_Mirau\W3.jpg'),3);
ImA=ImA/max(ImA(:));
ImB=sum(imread('D:\AMO\160713_Mirau\T1.jpg'),3);
ImB=ImB/max(ImB(:));

FFT_A=fft2(ImA-mean(ImA(:)));
FFT_B=fft2(ImB-mean(ImB(:)));

FFT_A(1:10,:)=0;

FFT_A(:,1:10)=0;

FFT_A((size(FFT_A,1)-9):size(FFT_A,1),:)=0;
FFT_A(:,(size(FFT_A,2)-9):size(FFT_A,2))=0;



FFT_B(1:10,:)=0;

FFT_B(:,1:10)=0;

FFT_B((size(FFT_B,1)-9):size(FFT_B,1),:)=0;
FFT_B(:,(size(FFT_B,2)-9):size(FFT_B,2))=0;


imagesc(abs(FFT_A))
imagesc(abs(FFT_B))
