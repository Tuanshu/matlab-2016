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

Power_BackOnly=[0.43
0.851
1.248
1.635
2.023
2.378
2.712
3.041
];

Power_withMirror=[0.505
1
1.47
1.93
2.402
2.83
3.24
3.63
];


Power_withAl=[0.78
1.7
];


Green_withAl=[78.8
234
];


figure('Position', [100, 100, 800, 500]);
plot(Green,Power_BackOnly,'-o',Green,Power_withMirror,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
hold on
plot(Green_withAl,Power_withAl,'o','Color','r');
legend('Backward Only','Backward + Forward (Mirror)','Backward + Forward (Aluminum Paper)','Location','NorthWest')
title('Output Ti:sapphire Power (mW @780nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)


figure('Position', [100, 100, 800, 500]);
plot(Green,Res_B40,'-o',Green,Res_A60,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('L3: 40x, B Wavelength Range','L3: 60x, A Wavelength Range','Location','NorthWest')
title('Residual Pump Power (mW @520nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)