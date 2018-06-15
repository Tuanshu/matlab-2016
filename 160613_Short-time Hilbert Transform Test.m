clear all

X=pi/2:pi/2:4*pi;

Y=sin(X)+0.5;

plot(X,Y)

S=abs(fft(Y));

plot(S);


XX=[1 4 6 4 1 4 6 4];   %bi nomial
XX=[1 4 6 4 1 4 6 4];

XX=[0 1 1 1 0 1 1 1];


XX=[0 1 2 1 0 1 2 1];


SS=fft(XX);


plot(abs(SS));