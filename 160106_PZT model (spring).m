clear all

q_forward=0:0.01:1;
q_backward=1:-0.01:0;

q=[q_forward q_backward];

C=[1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2 1/2];%[1/2 1/0.6 1/0.3 1/0.26 1/0.06 1/0.1 1/0.05 1/0.03 1/0.1 1/0.5];        %capacitance

v=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];%[0.2 0.3 0.3 2.6 0.9 2 1.5 1.2 7 80];        %frictional voltage

V_forward=ones(length(q_forward),length(C));
V_backward=ones(length(q_backward),length(C));
qb=zeros(1,length(C));
for p=1:length(q_forward)
    for s=1:length(C)
        if abs((q_forward(p)-qb(s))/C(s))<v(s)
            V_forward(p,s)=(q_forward(p)-qb(s))/C(s);
        else
            if p>2
                V_forward(p,s)=v(s)*sign(q_forward(p)-q_forward(p-1));
                qb(s)=q_forward(p)-C(s)*v(s)*sign(q_forward(p)-q_forward(p-1));
            else
                V_forward(p,s)=v(s)*sign(q_forward(p+1)-q_forward(p));
                qb(s)=q_forward(p)-C(s)*v(s)*sign(q_forward(p+1)-q_forward(p));
            end   
        end
    end
end

for p=1:length(q_backward)
    for s=1:length(C)
        if abs((q_backward(p)-qb(s))/C(s))<v(s)
            V_backward(p,s)=(q_backward(p)-qb(s))/C(s);
        else
            if p>2
                V_backward(p,s)=v(s)*sign(q_backward(p)-q_backward(p-1));
                qb(s)=q_backward(p)-C(s)*v(s)*sign(q_backward(p)-q_backward(p-1));
            else
                V_backward(p,s)=v(s)*sign(q_backward(p+1)-q_backward(p));
                qb(s)=q_backward(p)-C(s)*v(s)*sign(q_backward(p+1)-q_backward(p));
            end

        end
    end
end
V_forward_total=sum(V_forward,2);
V_backward_total=sum(V_backward,2);

plot(q_forward,V_forward_total,q_backward,V_backward_total);
    xlabel('q');
    ylabel('V');
    %xlim([-12 12]);
    %ylim([-30 30]);

%% ´«¦¨displacement

plot(q_forward,V_forward_total);
    xlabel('q');
    ylabel('V');
    %xlim([-12 12]);