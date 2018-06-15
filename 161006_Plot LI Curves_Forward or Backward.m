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


Res_Forward=[
];

Res_Backward=[
];


figure('Position', [100, 100, 800, 500]);
plot(Green,Power_Forward,'-o',Green,Power_Backward,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('Optimized Forward Power','Optimized Backward Power','Location','NorthWest')
title('Output Ti:sapphire Power (mW @780nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)


figure('Position', [100, 100, 800, 500]);
plot(Green,Res_Before,'-o',Green,Res_After,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('Backward Beam Not Collimated','Backward Beam  Collimated','Location','NorthWest')
title('Residual Pump Power (mW @520nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)