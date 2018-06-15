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

Power_Forward=[0.93
1.723
2.526
3.3
4.04
4.744
5.4
6
];

Power_Backward=[0.43
0.851
1.248
1.635
2.023
2.378
2.712
3.041
];


Res_B40=[11.09
22.02
32.34
42.59
53.11
62.51
71.62
79.79
];

Res_A60=[14.68
29.2
43.05
56.54
70.69
83.24
95.4
106.81
];


figure('Position', [100, 100, 800, 500]);
plot(Green,Power_Forward,'-o',Green,Power_Backward,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('Optimized Forward Power','Optimized Backward Power','Location','NorthWest')
title('Output Ti:sapphire Power (mW @780nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)


figure('Position', [100, 100, 800, 500]);
plot(Green,Res_B40,'-o',Green,Res_A60,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('L3: 40x, B Wavelength Range','L3: 60x, A Wavelength Range','Location','NorthWest')
title('Residual Pump Power (mW @520nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)