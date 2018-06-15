clear all

qb_initial=0;
q_start=0;
q_max=0.4;
q_min=-0.4;
d_q=0.01;
q_forward_1=q_start:d_q:q_max;
q_backward=(q_max-d_q):-d_q:q_min;
q_forward=(q_min+d_q):d_q:q_max;

q=[q_forward_1 q_backward q_forward];

C=[1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2];%[1/2 1/0.6 1/0.3 1/0.26 1/0.06 1/0.1 1/0.05 1/0.03 1/0.1 1/0.5];        %capacitance

v=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];%[0.2 0.3 0.3 2.6 0.9 2 1.5 1.2 7 80];        %frictional voltage

V=ones(length(q),length(C));
qb=ones(1,length(C))*qb_initial;
for p=1:length(q)
    for s=1:length(C)
        if abs((q(p)-qb(s))/C(s))<v(s)
            V(p,s)=(q(p)-qb(s))/C(s);
        else
            if p>2
                V(p,s)=v(s)*sign(q(p)-q(p-1));
                qb(s)=q(p)-C(s)*v(s)*sign(q(p)-q(p-1));
            else
                V(p,s)=v(s)*sign(q(p+1)-q(p));
                qb(s)=q(p)-C(s)*v(s)*sign(q(p+1)-q(p));
            end   
        end
    end
end
V_total=sum(V,2);

plot(q,V_total);
    xlabel('q');
    ylabel('V');
    %xlim([-12 12]);
    %ylim([-30 30]);

%% ´«¦¨displacement
