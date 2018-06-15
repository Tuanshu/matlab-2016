clear all
cd('D:\AMO\160926_TiSa LS Development\161006_Data Figures');
Current=[0.3
0.4
0.5
0.6
0.7
0.8
0.9
1
]';

Green=[88.9
179.1
265
350
435
513
587
660
];

Power_60x=[0.93
1.723
2.526
3.3
4.04
4.744
5.4
6
];

Power_40x=[0.737
1.41
2.047
2.658
3.251
3.785
4.29
4.76
];


Res_60x=[23.5
45.4
67.3
89
110.3
131
150.5
169
];

Res_40x=[20.21343
38.5299
56.60433
74.72062
92.96589
109.50615
124.8931
138.3364
];


figure('Position', [100, 100, 800, 500]);
plot(Green,Power_60x,'-o',Green,Power_40x,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('L3: 60x, B Wavelength Range','L3: 40x, B Wavelength Range','Location','NorthWest')
title('Output Ti:sapphire Power (mW @780nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)


figure('Position', [100, 100, 800, 500]);
plot(Green,Res_60x,'-o',Green,Res_40x,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('L3: 60x, B Wavelength Range','L3: 40x, B Wavelength Range','Location','NorthWest')
title('Residual Pump Power (mW @520nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)