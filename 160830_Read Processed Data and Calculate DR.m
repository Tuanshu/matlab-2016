clear all
fclose('all')
%%

last_folder_name='160830_PCO_4160fps_0x_200microsec_200x200_forearm_5_no scan';
%Axial_Decimation_Factor=[1 2 4 8 16 32 64];
Axial_Decimation_Factor=[1 2:2:64];
    for p=1:length(Axial_Decimation_Factor)
    Data_Save_Folder='F:\P1.2 Test\Processed Data\';

    Processed_Data_Path=[Data_Save_Folder last_folder_name sprintf('Axial_Decimation_Factor_%d.bin',Axial_Decimation_Factor(p))];


    Row=200;
    Colomn=200;


    fin = fopen(Processed_Data_Path);
    Image_Temp=fread(fin,[Row,inf],'double');
    fclose(fin);
    %
    Colomn_Total=size(Image_Temp,2);

    Frame=Colomn_Total/Colomn;
    Image_Stack=zeros(Row,Colomn,Frame);
    for r=1:Frame
        Image_Stack(:,:,r)=Image_Temp(:,(1+(r-1)*Colomn):(r*Colomn));
    end
    %

    X_Show=size(Image_Stack,1)/2;
    Y_Show=size(Image_Stack,2)/2;

    Array_Show=squeeze(Image_Stack(X_Show,Y_Show,:));

    plot(Array_Show);

    %
    Noise_Start_Index=100;
    Noise_End_Index=400;

    Glass_Interface_Position=214;

    Noise_STD(:,:)=std(Image_Stack(:,:,Noise_Start_Index:Noise_End_Index),0,3);

    imagesc(Noise_STD);

    Mean_Noise_STD(p)=mean(Noise_STD(:))
end

 plot(64./Axial_Decimation_Factor,Mean_Noise_STD*64./Axial_Decimation_Factor);
 xlabel('Averaging Factor');
 ylabel('Noise (Assume Sum-up Instead of Averaged)');
