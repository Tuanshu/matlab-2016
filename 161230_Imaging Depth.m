clear all

n=1.36;

DOFobj=4;


n0=1.3:0.001:1.4;

DOF=abs(DOFobj./((n-n0)+(1./n-1./n0)));

plot(n0,DOF,'linewidth',3);
xlabel('Oil/Gel Refractive Index','fontsize',15);
ylabel('Depth of Field (micron)','fontsize',15);
set(gca,'fontsize',15)
xlim([1.3 1.41]);
ylim([0 1000]);