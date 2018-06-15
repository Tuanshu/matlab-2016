clear all

Folder_path='D:\160503_PZT calibration related\';

X_f=dlmread([Folder_path 'X_f.txt']);
X_bf=dlmread([Folder_path 'X_bf.txt']);
V_f=dlmread([Folder_path 'V_f.txt']);
V_bf=dlmread([Folder_path 'V_bf.txt']);

%% �@�Ǫ�B�B�z

%�u������V

ind=find(V_f>=0,1,'first');

V_f=V_f(ind:end);
V_bf=V_bf(ind:end);
X_f=X_f(ind:end);
X_bf=X_bf(ind:end);

% X_f���̤p�Ȭ�0
X_f=X_f-min(X_f);

% ����X_bf��X_f�MX_bf form a close loop

Delta=X_bf-X_f;

Delta_linear=[Delta(1):(Delta(end)-Delta(1))/(length(Delta)-1):Delta(end)]';

X_bf=X_bf-Delta_linear;

plot(V_f,X_f,V_f,X_bf);

%%
Vbasis=min(V_f):(max(V_f)-min(V_f))/(length(V_f)-1):max(V_f);
X_f_Vbasis=X_f;
X_bf_Vbasis=X_bf;

If_increasing_f=diff(X_f_Vbasis)>=0;
ind_f=find(If_increasing_f==0,1,'last');

if ~isempty(ind_f)
    Temp=X_f_Vbasis(1:(ind_f*10));
    Index_Array=0:(length(Temp)-1);
    Index_Array_Downsample=(0:9)/9*(length(Temp)-1);
    Temp_downsample=interp1(Index_Array,Temp,Index_Array_Downsample);

    Temp_New=interp1(Index_Array_Downsample,Temp_downsample,Index_Array)';
    X_f_Vbasis(1:(ind_f*10))=Temp_New;
end

If_increasing_bf=diff(X_bf_Vbasis)>=0;
ind_bf=find(If_increasing_bf==0,1,'last');

if ~isempty(ind_bf)
    Temp=X_bf_Vbasis(1:(ind_bf*10));
    Index_Array=0:(length(Temp)-1);
    Index_Array_Downsample=(0:9)/9*(length(Temp)-1);
    Temp_downsample=interp1(Index_Array,Temp,Index_Array_Downsample);

    Temp_New=interp1(Index_Array_Downsample,Temp_downsample,Index_Array)';
    X_bf_Vbasis(1:(ind_bf*10))=Temp_New;
end

If_stric_X_f_Vbasis=[[1;diff(X_f_Vbasis)]>0];
Index_All_X_f_Vbasis=[1:length(If_stric_X_f_Vbasis)]';

Stric_X_f_Vbasis=X_f_Vbasis(If_stric_X_f_Vbasis);
Index_Stric_X_f_Vbasis=Index_All_X_f_Vbasis(If_stric_X_f_Vbasis);

Index_Stric_X_f_Vbasis(end)=Index_All_X_f_Vbasis(end);    %�`�N�o��, �ت��O���̫�@�ӥ��x���̫�@�I
    

X_f_Vbasis=interp1(Index_Stric_X_f_Vbasis,Stric_X_f_Vbasis,Index_All_X_f_Vbasis);
plot(Vbasis,X_f_Vbasis,Vbasis,X_bf_Vbasis);


%% 5/3 Change to uniform X basis (keep the same number of data)

Xbasis=[min(X_f_Vbasis):(max(X_f_Vbasis)-min(X_f_Vbasis))/(length(X_f_Vbasis)-1):max(X_f_Vbasis)]';

V_f_Xbasis=interp1(X_f_Vbasis,Vbasis,Xbasis,'linear','extrap');
V_bf_Xbasis=interp1(X_bf_Vbasis,Vbasis,Xbasis,'linear','extrap');

plot(Xbasis,V_f_Xbasis,Xbasis,V_bf_Xbasis);


dlmwrite([Folder_path 'Vbasis.txt'],Vbasis,'delimiter','\t','newline','pc');
dlmwrite([Folder_path 'Xbasis.txt'],Xbasis,'delimiter','\t','newline','pc');
dlmwrite([Folder_path 'X_f_Vbasis.txt'],X_f_Vbasis,'delimiter','\t','newline','pc');
dlmwrite([Folder_path 'X_bf_Vbasis.txt'],X_bf_Vbasis,'delimiter','\t','newline','pc');
dlmwrite([Folder_path 'V_f_Xbasis.txt'],V_f_Xbasis,'delimiter','\t','newline','pc');
dlmwrite([Folder_path 'V_bf_Xbasis.txt'],V_bf_Xbasis,'delimiter','\t','newline','pc');

%%
delta_V=V_f_Xbasis-V_bf_Xbasis;
V_com=0.5*(V_f_Xbasis+V_bf_Xbasis);



plot(Xbasis,V_com);

%% 2016/02/16 try to find the X', �]���O�n�D�Ϩ��, �|�Ψ�����interp1(V1,X1,X2)���Φ�(�S����!) > ���G�n���d���F
%V_com_uniform_min=min([V_com]);
%V_com_uniform_max=max([V_com]);
%V_com_uniform=V_com_uniform_min:(V_com_uniform_max-V_com_uniform_min)/(X_basis_NumOfSam-1):V_com_uniform_max;
%X_basis_uniformV=interp1(V_com,X_basis,V_com_uniform);

%plot(V_com,X_basis,V_com_uniform,X_basis_uniformV);


%Coefficient=1;

%X_n=interp1(V_com/max(V_com)*max(X_basis),X_basis,X_basis);
X_n=V_com/max(V_com)*max(Xbasis);


plot(Xbasis,V_com,X_n,V_com);

plot(Xbasis,X_n);      %�@�}�lx���w�q�n���|��o�ӹϳy���v�T, ���է���offset �n, �ڲש󦳳o��registration map�F


plot(X_n,V_f_Xbasis,X_n,V_bf_Xbasis);

%% 2016/02/16 To get the X(X') (uniform X')

X_n_uniform_min=min([X_n]);
X_n_uniform_max=max([X_n]);
X_n_uniform=X_n_uniform_min:(X_n_uniform_max-X_n_uniform_min)/(X_basis_NumOfSam-1):X_n_uniform_max;
X_uniformXn=interp1(X_n,X_basis,X_n_uniform);

plot(X_n,X_basis,X_n_uniform,X_uniformXn);


%% 2016/02/16 To decomposite  the X(X')

Delta_r=5;
Min_X_n=0;
Max_X_n=400;    %���̤j�u����385


r_array=Min_X_n:Delta_r:Max_X_n-(Delta_r);%0:delta_V:max(V_measured)-delta_V;  %[0 0.1 0.2 0.3 0.4 0.5 1 2 4 6 8 9 9.5 9.6 9.7 9.8 9.9];

S_array=zeros(length(r_array),1);   
%forward
index_r=1;
for p=1:length(r_array)
    if p==1
        index_r_next=find(X_n_uniform>r_array(p+1),1,'first');
        S_array(p)=(X_uniformXn(index_r_next)-X_uniformXn(index_r))/(X_n_uniform(index_r_next)-X_n_uniform(index_r));
    elseif p<length(r_array)
        index_r_next=find(X_n_uniform>r_array(p+1),1,'first');
        S_array(p)=(X_uniformXn(index_r_next)-X_uniformXn(index_r))/(X_n_uniform(index_r_next)-X_n_uniform(index_r))-sum(S_array(1:p-1));
    else
        index_r_next=length(X_n_uniform);
        S_array(p)=(X_uniformXn(index_r_next)-X_uniformXn(index_r))/(X_n_uniform(index_r_next)-X_n_uniform(index_r))-sum(S_array(1:p-1));
    end

    index_r=index_r_next;
end
plot(r_array,S_array);
%% 2016/2/16 Test: Given X_n, calculate X
X_n_Given=X_n_uniform;
X_Calculated=zeros(length(X_n_Given),length(r_array));
for s=1:(length(r_array)) %s=NNNNN:NNNNN%
    if r_array(s)==0
        X_Calculated(:,s)=S_array(s).*X_n_Given;
	elseif r_array(s)>0
        X_Calculated(:,s)=S_array(s).*max(X_n_Given-r_array(s),0);
	elseif r_array(s)<0
        X_Calculated(:,s)=S_array(s).*max(X_n_Given-r_array(s),0);
	end
end

plot(X_n_Given,sum(X_Calculated,2),X_n_Given,X_uniformXn);
%% �n���|�Ψ�: To decomposite the X'(X)

Delta_r_rev=5;
Min_X=0;
Max_X=400;    %���̤j�u����385


r_rev_array=Min_X:Delta_r_rev:Max_X-(Delta_r_rev);%0:delta_V:max(V_measured)-delta_V;  %[0 0.1 0.2 0.3 0.4 0.5 1 2 4 6 8 9 9.5 9.6 9.7 9.8 9.9];

S_rev_array=zeros(length(r_rev_array),1);   
%forward
index_r_rev=1;
for p=1:length(r_rev_array)
    if p==1
        index_r_rev_next=find(X_basis>r_rev_array(p+1),1,'first');
        S_rev_array(p)=(X_n(index_r_rev_next)-X_n(index_r_rev))/(X_basis(index_r_rev_next)-X_basis(index_r_rev));
    elseif p<length(r_rev_array)
        index_r_rev_next=find(X_basis>r_rev_array(p+1),1,'first');
        S_rev_array(p)=(X_n(index_r_rev_next)-X_n(index_r_rev))/(X_basis(index_r_rev_next)-X_basis(index_r_rev))-sum(S_rev_array(1:p-1));
    else
        index_r_next=length(X_n_uniform);
        S_rev_array(p)=(X_n(index_r_rev_next)-X_n(index_r_rev))/(X_basis(index_r_rev_next)-X_basis(index_r_rev))-sum(S_rev_array(1:p-1));
    end

    index_r_rev=index_r_rev_next;
end
plot(r_rev_array,S_rev_array);


%% �n���|�Ψ�: Given X, calculate X_n (�|�ݰ_�ӹ諸������n�O�]����X�Ӫ��_�I�@�w�b0����(�Linitial offset), ��X_n�|��, �]���O���Vcom��X�Ӫ�, �����u�O�toffset�Ӥw)
X_Given=X_basis;
X_n_Calculated=zeros(length(X_Given),length(r_rev_array));
for s=1:(length(r_rev_array)) %s=NNNNN:NNNNN%
    if r_rev_array(s)==0
        X_n_Calculated(:,s)=S_rev_array(s).*X_Given;
	elseif r_rev_array(s)>0
        X_n_Calculated(:,s)=S_rev_array(s).*max(X_Given-r_rev_array(s),0);
	elseif r_rev_array(s)<0
        X_n_Calculated(:,s)=S_rev_array(s).*max(X_Given-r_rev_array(s),0);
	end
end
%plot(X_n_Given,sum(X_Calculated,2),X_n_Given,X_uniformXn);

plot(X_Given,sum(X_n_Calculated,2),X_Given,X_n);

%% 2016/2/17 �ΤW�zfunction generate X_for_n�MX_bak_n                    (�ݷ|���K
X_f_n_Calculated=zeros(length(X_f),length(r_rev_array));
X_bf_n_Calculated=zeros(length(X_bf),length(r_rev_array));
for s=1:(length(r_rev_array)) %s=NNNNN:NNNNN%
    if r_rev_array(s)==0
        X_f_n_Calculated(:,s)=S_rev_array(s).*X_f;
        X_bf_n_Calculated(:,s)=S_rev_array(s).*X_bf;
	elseif r_rev_array(s)>0
        X_f_n_Calculated(:,s)=S_rev_array(s).*max(X_f-r_rev_array(s),0);
        X_bf_n_Calculated(:,s)=S_rev_array(s).*max(X_bf-r_rev_array(s),0);
	elseif r_rev_array(s)<0
        X_f_n_Calculated(:,s)=S_rev_array(s).*max(X_f-r_rev_array(s),0);
        X_bf_n_Calculated(:,s)=S_rev_array(s).*max(X_bf-r_rev_array(s),0);
	end
end

plot(V_f,X_f,V_f,X_bf,V_f,sum(X_f_n_Calculated,2),V_f,sum(X_bf_n_Calculated,2));
plot(V_f,sum(X_f_n_Calculated,2),V_f,sum(X_bf_n_Calculated,2),V_f,0.5*(sum(X_f_n_Calculated,2)+sum(X_bf_n_Calculated,2)));  %%�_��, �������ᱵ����I�B���I����?

%% 2016/2/16 ���U�ӨD���X'(V)��w�Mvth
Delta_vth=1;
Min_V=0;
Max_V=10;


vth_array=Min_V:Delta_vth:Max_V-(Delta_vth);%0:delta_V:max(V_measured)-delta_V;  %[0 0.1 0.2 0.3 0.4 0.5 1 2 4 6 8 9 9.5 9.6 9.7 9.8 9.9];


V_measured_forward=V_f;
X_measured_forward=X_f_Ratioed_adjusted;%Position_Patial_Forward;
V_measured_backward=V_f;
X_measured_backward=X_bf_Ratioed_adjusted;%Position_Patial_Backward;


%plot(V_measured_backward,X_measured_backward,V_measured_forward,X_measured_forward,V_measured_forward,X_mean_ideal);
%xlabel('V');
%ylabel('X');


w_measured_forward=zeros(length(vth_array),1);   
w_measured_backward=zeros(length(vth_array),1);   
%forward
index_vth=1;
for p=1:length(vth_array)
    if p==1
        index_vth_next=find(V_measured_forward>vth_array(p+1),1,'first');
        w_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth));
    elseif p<length(vth_array)
        index_vth_next=find(V_measured_forward>vth_array(p+1),1,'first');
        w_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth))-sum(w_measured_forward(1:p-1));
    else
        index_vth_next=length(V_measured_forward);
        w_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth))-sum(w_measured_forward(1:p-1));
    end

    index_vth=index_vth_next;
end
%backward 16/1/27 NEW!!!!!!!!!!!!!! FORCE�Ĥ@����w�Pforward�ۦP
index_vth=length(X_measured_backward);
for p=1:length(vth_array)
    if p==1
        index_vth_next=find((V_measured_backward)<(Max_V-vth_array(p+1)),1,'last');
        w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth));%w_measured_forward(p);%
    elseif p<length(vth_array)
        index_vth_next=find((V_measured_backward)<(Max_V-vth_array(p+1)),1,'last');
        w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth))-sum(w_measured_backward(1:p-1));
    else
        index_vth_next=1;
        w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth))-sum(w_measured_backward(1:p-1));
    end

    index_vth=index_vth_next;
end

%w_measured_forward(w_measured_forward<0)=0;
%w_measured_backward(w_measured_backward<0)=0;



%%



Center_X=195;

Xc=abs(X_basis-Center_X);

Xc_no_abs=X_basis-Center_X;

index_Xc_went_positive=find(Xc_no_abs>0,1,'first');

Xc_1=Xc(1:index_Xc_went_positive);
Xc_2=Xc((index_Xc_went_positive+1):end);



plot(Xc,delta_V);

%% �]����narray���׬�����, �ڴN�����Ҽ{�_�ƪ����p�F
%
Amp=1;
Offset=2245;        %Check: �ORatio�j�P���k���; ��������Ratio_Check���n�X�{0
clear Ratio
for p=1:((length(X_diff))-1-Offset)
    Ratio(p)=X_diff(length(X_diff)-p-Offset)/X_diff(p);
end
Ratio_root=Amp.*(Ratio).^0.5;
Ratio_root(length(Ratio_root):length(X_diff))=0;        %�|���@��0, �����O�]��X_diff�e�X���]��0
X_diff_Better=X_diff.*Ratio_root';
plot(1:length(Ratio_root),Ratio_root);
%%
plot(X_diff_Better);

% ���
Offset=2250;

for p=1:((length(X_diff_Better))-1-Offset)
    Ratio_Check(p)=X_diff_Better(length(X_diff_Better)-p-Offset)/X_diff_Better(p);
end
plot(1:length(Ratio_Check),Ratio_Check);
%%
%plot(1:length(Ratio),Ratio,1:length(X_diff),X_diff)
% x'=a*x+b��b, offset correction

a=Ratio_root';

X_f_Ratioed=a.*X_f;
X_bf_Ratioed=a.*X_bf;

XX=X_bf_Ratioed-X_f_Ratioed;
plot(XX);
X_mean_ratioed=a.*(X_f+X_bf)/2;
[maxvalue maxindex]=max(X_mean_ratioed);

X_mean_ratioed_ideal=[X_mean_ratioed(1):(maxvalue-X_mean_ratioed(1))/(maxindex-1):X_mean_ratioed(1)+(maxvalue-X_mean_ratioed(1))/(maxindex-1)*(length(X_mean_ratioed)-1)]';
plot(X_mean_ratioed_ideal);
X_mean_ratioed_adjust_to_ideal=X_mean_ratioed_ideal-X_mean_ratioed;
b=X_mean_ratioed_adjust_to_ideal;
X_f_Ratioed_adjusted=X_f_Ratioed+b;
X_bf_Ratioed_adjusted=X_bf_Ratioed+b;


plot(1:length(X_f_Ratioed),X_f_Ratioed_adjusted,1:length(X_bf_Ratioed),X_bf_Ratioed_adjusted);
plot(X_f,X_f_Ratioed_adjusted,X_bf,X_bf_Ratioed_adjusted);
%% (16/02/02 NEW!) To Decomposing the V-X relation    w, forward, measured based on X_f_Ratioed_adjusted and X_bf_Ratioed_adjusted)
Delta_vth=1;
Min_V=0;
Max_V=10;


vth_array=Min_V:Delta_vth:Max_V-(Delta_vth);%0:delta_V:max(V_measured)-delta_V;  %[0 0.1 0.2 0.3 0.4 0.5 1 2 4 6 8 9 9.5 9.6 9.7 9.8 9.9];


V_measured_forward=V_f;
X_measured_forward=X_f_Ratioed_adjusted;%Position_Patial_Forward;
V_measured_backward=V_f;
X_measured_backward=X_bf_Ratioed_adjusted;%Position_Patial_Backward;


%plot(V_measured_backward,X_measured_backward,V_measured_forward,X_measured_forward,V_measured_forward,X_mean_ideal);
%xlabel('V');
%ylabel('X');


w_measured_forward=zeros(length(vth_array),1);   
w_measured_backward=zeros(length(vth_array),1);   
%forward
index_vth=1;
for p=1:length(vth_array)
    if p==1
        index_vth_next=find(V_measured_forward>vth_array(p+1),1,'first');
        w_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth));
    elseif p<length(vth_array)
        index_vth_next=find(V_measured_forward>vth_array(p+1),1,'first');
        w_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth))-sum(w_measured_forward(1:p-1));
    else
        index_vth_next=length(V_measured_forward);
        w_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth))-sum(w_measured_forward(1:p-1));
    end

    index_vth=index_vth_next;
end
%backward 16/1/27 NEW!!!!!!!!!!!!!! FORCE�Ĥ@����w�Pforward�ۦP
index_vth=length(X_measured_backward);
for p=1:length(vth_array)
    if p==1
        index_vth_next=find((V_measured_backward)<(Max_V-vth_array(p+1)),1,'last');
        w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth));%w_measured_forward(p);%
    elseif p<length(vth_array)
        index_vth_next=find((V_measured_backward)<(Max_V-vth_array(p+1)),1,'last');
        w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth))-sum(w_measured_backward(1:p-1));
    else
        index_vth_next=1;
        w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth))-sum(w_measured_backward(1:p-1));
    end

    index_vth=index_vth_next;
end

%w_measured_forward(w_measured_forward<0)=0;
%w_measured_backward(w_measured_backward<0)=0;

plot(1:length(w_measured_forward),w_measured_forward,1:length(w_measured_forward),w_measured_backward);

%%

% plot(Voltage_Patial_Forward,Position_Patial_Forward,flipud(Voltage_Patial_Backward),flipud(Position_Patial_Backward));
% 
% 
% V_f=Voltage_Patial_Forward;
% X_f=Position_Patial_Forward;
% V_bf=flipud(Voltage_Patial_Backward);
% X_bf=flipud(Position_Patial_Backward);
% 
% if length(X_f)>length(X_bf)
%     X_bf((length(X_bf)+1):length(X_f))=X_bf(length(X_bf));
% else
%     X_f((length(X_f)+1):length(X_bf))=X_f(length(X_f));
% end
% plot(V_f,X_f,V_f,X_bf);
% 
% %V_diff=V_bf-V_f;
% X_diff=X_bf-X_f;
% 
% X_mean=(X_bf+X_f)/2;
% 
% X_mean_ideal=[X_mean(1):(X_mean(end)-X_mean(1))/(length(X_mean)-1):X_mean(end)]';
% 
% plot(V_f,X_diff);
% plot(V_f,X_mean,V_f,X_mean_ideal);
% X_f_better=X_mean_ideal-X_diff/2;
% X_bf_better=X_mean_ideal+X_diff/2;
% 
% 
% X_b_better=flipud(X_bf_better);
% 
% if length(Voltage_Patial_Backward)>length(X_b_better)
%     X_b_better((length(X_b_better)+1):length(Voltage_Patial_Backward))=X_b_better(length(X_b_better));
% else
%     X_b_better=X_b_better(1:length(Voltage_Patial_Backward));
% end
% 
% 
% plot(V_f,X_f_better,V_f,X_bf_better,V_f,X_f,V_f,X_bf);
% legend('X_f_better','X_bf_better','X_f','X_bf');
%% (16/02/01 NEW!) To Decomposing the V-X relation    w, forward, measured based on V_f, X_f_better, X_bf_better (note! not X_b_better)
% 
% Delta_vth=1;
% Min_V=0;
% Max_V=10;
% 
% 
% vth_array=Min_V:Delta_vth:Max_V-(Delta_vth);%0:delta_V:max(V_measured)-delta_V;  %[0 0.1 0.2 0.3 0.4 0.5 1 2 4 6 8 9 9.5 9.6 9.7 9.8 9.9];
% 
% 
% V_measured_forward=V_f;
% X_measured_forward=X_f_better;%Position_Patial_Forward;
% V_measured_backward=V_f;
% X_measured_backward=X_bf_better;%Position_Patial_Backward;
% 
% 
% plot(V_measured_backward,X_measured_backward,V_measured_forward,X_measured_forward,V_measured_forward,X_mean_ideal);
% xlabel('V');
% ylabel('X');
% 
% 
% w_measured_forward=zeros(length(vth_array),1);   
% w_measured_backward=zeros(length(vth_array),1);   
% %forward
% index_vth=1;
% for p=1:length(vth_array)
%     if p==1
%         index_vth_next=find(V_measured_forward>vth_array(p+1),1,'first');
%         w_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth));
%     elseif p<length(vth_array)
%         index_vth_next=find(V_measured_forward>vth_array(p+1),1,'first');
%         w_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth))-sum(w_measured_forward(1:p-1));
%     else
%         index_vth_next=length(V_measured_forward);
%         w_measured_forward(p)=(X_measured_forward(index_vth_next)-X_measured_forward(index_vth))/(V_measured_forward(index_vth_next)-V_measured_forward(index_vth))-sum(w_measured_forward(1:p-1));
%     end
% 
%     index_vth=index_vth_next;
% end
% %backward 16/1/27 NEW!!!!!!!!!!!!!! FORCE�Ĥ@����w�Pforward�ۦP
% index_vth=length(X_measured_backward);
% for p=1:length(vth_array)
%     if p==1
%         index_vth_next=find((V_measured_backward)<(Max_V-vth_array(p+1)),1,'last');
%         w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth));%w_measured_forward(p);%
%     elseif p<length(vth_array)
%         index_vth_next=find((V_measured_backward)<(Max_V-vth_array(p+1)),1,'last');
%         w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth))-sum(w_measured_backward(1:p-1));
%     else
%         index_vth_next=1;
%         w_measured_backward(p)=(X_measured_backward(index_vth_next)-X_measured_backward(index_vth))/(V_measured_backward(index_vth_next)-V_measured_backward(index_vth))-sum(w_measured_backward(1:p-1));
%     end
% 
%     index_vth=index_vth_next;
% end
% 
% %w_measured_forward(w_measured_forward<0)=0;
% %w_measured_backward(w_measured_backward<0)=0;
% 
% plot(1:length(w_measured_forward),w_measured_forward,1:length(w_measured_forward),w_measured_backward);
%% Calibrate with 2nd loop of unipolar poling (�]�N�O�Dinitial loading curve, �D�o���Ow�M2*vth)
V_start=0.01;
V_max_1=10;
V_max_2=10;
V_max_3=10;
V_min_1=0;
V_min_2=0;
V_min_3=0;
d_V=0.01;
V_forward_1=V_start:d_V:V_max_1;
V_backward_1=(V_max_1-d_V):-d_V:V_min_1;
V_forward_2=(V_min_1+d_V):d_V:V_max_2;
V_backward_2=(V_max_2-d_V):-d_V:V_min_2;
V_forward_3=(V_min_2+d_V):d_V:V_max_3;
V_backward_3=(V_max_3-d_V):-d_V:V_min_3;

V=[V_forward_1 V_backward_1];% V_forward_2 V_backward_2 V_forward_3 V_backward_3];
plot(V)

Rate=0.01; %V/sec
Time=Rate:Rate:Rate*length(V);
vth_measured_another_assumption=vth_array/2;
X_another_assumption=zeros(length(V),length(w_measured_forward));
NNNNN=10;
for s=1:(length(w_measured_forward)) %s=NNNNN:NNNNN%
    for p=1:length(V)
        if p==1
            X_another_assumption(p,s)=w_measured_backward(s)*(V(p)+vth_measured_another_assumption(s));
            %X_another_assumption(p,s)=w_measured_forward(s)*(V(p)+vth_measured_another_assumption(s));
        else
            X_another_assumption(p,s)=max(w_measured_forward(s)*(V(p)-vth_measured_another_assumption(s)),min(w_measured_backward(s)*(V(p)+vth_measured_another_assumption(s)),X_another_assumption(p-1,s)));
            %X_another_assumption(p,s)=max(w_measured_forward(s)*(V(p)-vth_measured_another_assumption(s)),min(w_measured_forward(s)*(V(p)+vth_measured_another_assumption(s)),X_another_assumption(p-1,s)));
        end
    end
end
X_total_another_assumption=sum(X_another_assumption,2);

plot(V,X_another_assumption,V,X_total_another_assumption);
xlabel('V (volt)');
ylabel('X (micron)');
legend('X_1','X_2','X_3','X_4','X_5','X_6','X_7','X_8','X_9','X_1_0','X_1_1','X_1_2','X_t_o_t_a_l');
plot(V,X_total_another_assumption,V_f,X_f_Ratioed_adjusted,V_f,X_bf_Ratioed_adjusted);
xlabel('V (volt)');
ylabel('X (micron)');
%%
plot(V,X_total_another_assumption-(X_total_another_assumption(1))+Position_Patial_Roundtrip(1),V_measured_forward,X_measured_forward,V_measured_backward,X_measured_backward);
plot(V,X_total_another_assumption-(X_total_another_assumption(1))+Position_Patial_Roundtrip(1),Voltage_Patial_Roundtrip,Position_Patial_Roundtrip);

        xlabel('V (volt)');
    ylabel('X (micron)');
    legend('V-X (calculated)','V-X (measured)');
    
    plot(V,X_total_another_assumption-(X_total_another_assumption(1))+X_measured_forward(1),V_measured_forward,X_measured_forward,V_measured_backward,X_measured_backward);
%% Trial: try to transfrom the X_total_another_assumption to the original one

%1st step: to break down the X_total_another_assumption
X_total_another_assumption_f=X_total_another_assumption(1:length(V_forward_1));
X_total_another_assumption_b=X_total_another_assumption((length(V_forward_1)+1):end);
X_total_another_assumption_bf=flipud(X_total_another_assumption_b);
V_f_Trial=V_forward_1';
V_bf_Trial=flipud(V_backward_1');
plot(X_total_another_assumption_bf);
[maxvalue maxindex]=max(V_f);
V_f_for_interpolation=[V_f(1):(maxvalue-V_f(1))/(maxindex-1):V_f(1)+(maxvalue-V_f(1))/(maxindex-1)*(length(V_f)-1)]';
a_Trial=interp1(V_f_for_interpolation,a,V_f_Trial);
b_Trial=interp1(V_f_for_interpolation,b,V_f_Trial);
plot(V_f_for_interpolation,a,V_f_Trial,a_Trial);
X_f_Derived_Original=(X_total_another_assumption_f-(X_total_another_assumption_f(1))+X_f(1)-b_Trial)./a_Trial;
X_bf_Derived_Original=(X_total_another_assumption_bf-(X_total_another_assumption_bf(1))+X_bf(1)-b_Trial)./a_Trial;

%plot(V_f_Trial,X_total_another_assumption_f-X_total_another_assumption_f(1)+X_f_Ratioed_adjusted(1),V_f_Trial,X_total_another_assumption_bf-X_total_another_assumption_bf(1)+X_bf_Ratioed_adjusted(1),V_f,X_f_Ratioed_adjusted,V_f,X_bf_Ratioed_adjusted);%,V_f,X_f,V_f,X_bf);
plot(V_f_Trial,X_total_another_assumption_f,V_f_Trial,X_total_another_assumption_bf,V_f,X_f_Ratioed_adjusted,V_f,X_bf_Ratioed_adjusted);%,V_f,X_f,V_f,X_bf);
X_Test_f=(X_f_Ratioed_adjusted-b)./a;
X_Test_bf=(X_bf_Ratioed_adjusted-b)./a;
X_Trial_f=(X_total_another_assumption_f-X_total_another_assumption_f(1)+X_f_Ratioed_adjusted(1)-b_Trial)./a_Trial;
X_Trial_bf=(X_total_another_assumption_bf-X_total_another_assumption_bf(1)+X_bf_Ratioed_adjusted(1)-b_Trial)./a_Trial;
plot(V_f,X_Test_f,V_f,X_Test_bf);

plot(V_f_Trial,X_Trial_f,V_f_Trial,X_Trial_bf,V_f,X_Test_f,V_f,X_Test_bf);
% �ҥH���׬O: OK, �������`�N�@�}�l��offset, �]���Ounipolar
%% �M��nbreakdown a and b array, to match the form of vth etc.
a_array=zeros([length(vth_array) 1]);
b_array=zeros([length(vth_array) 1]);
index_vth=1;
for p=1:length(vth_array)
    if p<length(vth_array)
        index_vth_next=find(V_f>vth_array(p+1),1,'first');
        a_array(p)=mean(a(index_vth:index_vth_next));
        b_array(p)=min(b(index_vth:index_vth_next));
    else
        index_vth_next=length(V_f);
        a_array(p)=mean(a(index_vth:index_vth_next));
        b_array(p)=min(b(index_vth:index_vth_next));
        
    end
    index_vth=index_vth_next;
end
%% Note: �ѩ�x'=ax+b, �o��transform�Oapply�bx�W, �G�|���@�ix versus x'���@��, �b�ӹϤW, �i�H��������w�Mvth���覡, ��XS function���Ѽ�
%% (2/5 ��M�o�{, �ڤ@�}�l��a�Mb���w�q�i��N���Ӧn) > �o�N��ڥ����Otransform��?
plot(X_f_Derived_Original-X_f_Derived_Original(1),X_total_another_assumption_f-X_total_another_assumption_f(1),X_bf_Derived_Original-X_bf_Derived_Original(1),X_total_another_assumption_bf-X_total_another_assumption_bf(1));
for s=1:(length(a_array)) %s=NNNNN:NNNNN%   %�C�@��
    for p=1:length(V)
        if p==1
            X_another_assumption_back_to_Original(p,s)=w_measured_backward(s)*(V(p)+vth_measured_another_assumption(s));
            %X_another_assumption(p,s)=w_measured_forward(s)*(V(p)+vth_measured_another_assumption(s));
        else
            X_another_assumption_back_to_Original(p,s)=max(w_measured_forward(s)*(V(p)-vth_measured_another_assumption(s)),min(w_measured_backward(s)*(V(p)+vth_measured_another_assumption(s)),X_another_assumption(p-1,s)));
            %X_another_assumption(p,s)=max(w_measured_forward(s)*(V(p)-vth_measured_another_assumption(s)),min(w_measured_forward(s)*(V(p)+vth_measured_another_assumption(s)),X_another_assumption(p-1,s)));
        end
    end
end





        xlabel('V (volt)');
    ylabel('X (micron)');
    legend('V-X (calculated)','V-X (measured)');

    %%
    NNN=2;
    plot(X_another_assumption(:,NNN));