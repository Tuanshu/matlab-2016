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

Power_2micron=[0.93
1.723
2.526
3.3
4.04
4.744
5.4
6
];

Power_1micron=[1.246
2.502
3.724
4.933
6.143
7.279
8.371
9.4
];


Power_withAl=[0.78
1.7
];


Green_withAl=[78.8
234
];


figure('Position', [100, 100, 800, 500]);
plot(Green,Power_2micron,'-o',Green,Power_1micron,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('Side Coating: 2 micron (New Recipe)','Side Coating: 1 micron (Old Recipe)','Location','NorthWest')
title('Output Ti:sapphire Power (mW @780nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)


figure('Position', [100, 100, 800, 500]);
plot(Green,Res_B40,'-o',Green,Res_A60,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('L3: 40x, B Wavelength Range','L3: 60x, A Wavelength Range','Location','NorthWest')
title('Residual Pump Power (mW @520nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)