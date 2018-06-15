clear all


%t1= @ (n) (2*n1_Considered)./(n1_Considered+n); 

V_start=0;
V_max_1=10;
V_max_2=8;
V_max_3=6;
V_min_1=-10;
V_min_2=0;
V_min_3=0;
d_V=0.1;
V_forward_1=V_start:d_V:V_max_1;
V_backward_1=(V_max_1-d_V):-d_V:V_min_1;
V_forward_2=(V_min_1+d_V):d_V:V_max_2;
V_backward_2=(V_max_2-d_V):-d_V:V_min_2;
V_forward_3=(V_min_2+d_V):d_V:V_max_3;
V_backward_3=(V_max_3-d_V):-d_V:V_min_3;

Rate=0.01; %V/sec

V=[V_forward_1 V_backward_1 V_forward_2]'% V_backward_2 V_forward_3 V_backward_3];
plot(V)
w=[5];% 5];% 5 5 5 5 5];%slope
Time=Rate:Rate:Rate*length(V);
vth=[2];% 4]%; 5 6 7 8 9];%threshold 

X=zeros(length(V),length(w));
for s=1:length(w)
    for p=1:length(V)
        if p==1
            X(p,s)=max(w(s)*(V(p)-vth(s)),min(w(s)*(V(p)+vth(s)),0));
        else
            X(p,s)=max(w(s)*(V(p)-vth(s)),min(w(s)*(V(p)+vth(s)),X(p-1,s)));
        end
    end
end
X_total=sum(X,2);

plot(V,X,V,X_total);
plot(V,X_total);

    xlabel('V (volt)');
    ylabel('X (micron)');
    %xlim([-0.5 11]);
    %ylim([-5 150]);
%%
    plot(Time,V);
    xlabel('Time (second)');
    ylabel('V (volt)');
    
    plot(Time,X_total);
    xlabel('Time (second)');
    ylabel('Displacement (micron)');
%legend('X_1','X_2','X_3','X_4','X_t_o_t_a_l');
%% ´«¦¨displacement