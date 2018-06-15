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

Power_25Deg=[0.585
1.12
1.63
2.12
2.604
3.04
3.465
3.875
];

Power_65Deg=[0.777
1.495
2.17
2.8
3.455
4.01
4.54
5
];

Power_90Deg=[0.93
1.723
2.526
3.3
4.04
4.744
5.4
6
];



Res_25Deg=[30.55815
58.4968
85.7757
112.8868
140.10756
165.1056
189.42135
212.76125
];
Res_65Deg=[21.34903
41.19305
59.9063
78.492
97.03745
114.9439
130.6906
145.95
];
Res_90Deg=[22.0027
42.62597
63.23314
83.687
103.7956
123.36216
141.806
159.34
];

figure('Position', [100, 100, 800, 500]);
plot(Green,Power_25Deg,'-o',Green,Power_65Deg,'-o',Green,Power_90Deg,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('Theta = 25 Degree','Theta = 65 Degree','Theta = 90 Degree','Location','NorthWest')
title('Output Ti:sapphire Power (mW @780nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)



figure('Position', [100, 100, 800, 500]);
plot(Green,Res_25Deg,'-o',Green,Res_65Deg,'-o',Green,Res_90Deg,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('Theta = 25 Degree','Theta = 65 Degree','Theta = 90 Degree','Location','NorthWest')
title('Residual Pump Power (mW @520nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)