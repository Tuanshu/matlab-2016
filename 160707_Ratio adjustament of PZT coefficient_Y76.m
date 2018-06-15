clear all

Supposed_N=64;

Measured_N=79.857;


Old_a0=-1.7256*10^(-1);

Old_a1=4.6*10^(-3);

Old_a2=6.039*10^(-5);


New_a0=Old_a0

New_a1=Old_a1*(Measured_N/Supposed_N)^0.5

New_a2=Old_a2*Measured_N/Supposed_N
