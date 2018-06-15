clear all
close all
%Xinc=[-45.7 -25.95 29.92 9 -0.006 2.24 1.42];
Xinc=[14.62 4.08 -0.39];% -0.59 -0.18];
%Yinc=[-6.64 -14.17 -1.22 -9.3 -4.01 -2.83 -0.52];
Yinc=[-3.28 -4.72 0.85];% 0.7 0.29];
Xinc_micron=Xinc*24/5;
Yinc_micron=Yinc*32/5;
figure('Units','normalized','Position',[.1 .1 .5 .5])
plot(Xinc_micron,Yinc_micron,'.');
axis equal
grid on
grid minor
xlim([-100 100]);
ylim([-100 100]);
xlabel('Z-offset along X');
ylabel('Z-offset along Y');
%set(gca, 'XColor', 'gray');
%set(gca, 'YColor', 'g');

factor=1.1;

hold on



for p=1:(length(Xinc)-1)
    X_center=Xinc_micron(p);%(Xinc(p)+Xinc(p+1))/2;
    Y_center=Yinc_micron(p);%(Yinc(p)+Yinc(p+1))/2;
    X_length=(Xinc_micron(p+1)-Xinc_micron(p))*factor;
    Y_length=(Yinc_micron(p+1)-Yinc_micron(p))*factor;
    quiver(X_center,Y_center,X_length,Y_length);
end

for p=1:length(Xinc)
    txt1 = sprintf('Total=%.2f',abs(Xinc_micron(p))+abs(Yinc_micron(p)));
    text(Xinc_micron(p),Yinc_micron(p),txt1,'FontSize',12,'Color','b')
    
end

hold off