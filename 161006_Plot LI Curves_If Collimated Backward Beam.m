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

Power_Before=[0.49
0.956
1.4
1.82
2.24
2.61
2.96
3.3
];

Power_After=[0.44
0.873
1.28
1.679
2.084
2.445
2.8
3.138
];


Res_Before=[16.905763
32.6659172
47.67518
62.203834
77.397488
91.634207
104.450752
116.62371
];

Res_After=[15.642928
31.3806451
46.203136
61.5747573
76.9335308
91.4555215
105.32536
118.8326006
];


figure('Position', [100, 100, 800, 500]);
plot(Green,Power_Before,'-o',Green,Power_After,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('Backward Beam Not Collimated','Backward Beam Collimated','Location','NorthWest')
title('Output Ti:sapphire Power (mW @780nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)


figure('Position', [100, 100, 800, 500]);
plot(Green,Res_Before,'-o',Green,Res_After,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('Backward Beam Not Collimated','Backward Beam  Collimated','Location','NorthWest')
title('Residual Pump Power (mW @520nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)