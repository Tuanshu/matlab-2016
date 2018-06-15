clear all

q_forward=-10:0.1:10;
q_backward=10:-0.1:-10;

C=1;        %capacitance

v=10;        %frictional voltage

V_forward=ones(length(q_forward),1);
V_backward=ones(length(q_backward),1);

qb=0;
for p=1:length(q_forward)
    if (q_forward(p)-qb)/C<v
        V_forward(p)=(q_forward(p)-qb)/C;
    else
        if p>2
            V_forward(p)=v*sign(q_forward(p)-q_forward(p-1));
            qb=q_forward(p)-C*v*sign(q_forward(p)-q_forward(p-1));
        else
            V_forward(p)=v*sign(q_forward(p+1)-q_forward(p));
            qb=q_forward(p)-C*v*sign(q_forward(p+1)-q_forward(p));
        end   
    end
end
    
for p=1:length(q_backward)
    if (q_backward(p)-qb)/C<v
        V_backward(p)=(q_backward(p)-qb)/C;
    else
        if p>2
            V_backward(p)=v*sign(q_backward(p)-q_backward(p-1));
            qb=q_backward(p)-C*v*sign(q_backward(p)-q_backward(p-1));
        else
            V_backward(p)=v*sign(q_backward(p+1)-q_backward(p));
            qb=q_backward(p)-C*v*sign(q_backward(p+1)-q_backward(p));
        end

    end
end

plot(q_forward,V_forward,q_backward,V_backward);
    xlabel('q');
    ylabel('V');
    xlim([-12 12]);
    ylim([-20 20]);
