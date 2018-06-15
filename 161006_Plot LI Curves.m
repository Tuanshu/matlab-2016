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

Power_1micron_1=[0.55
1
1.5
1.93
2.36
2.75
3.1
3.46
];
Power_1micron_2=[0.545
1.028
1.495
1.94
2.38
2.77
3.12
3.47
];
Power_1micron_3=[0.545
1.04
1.505
1.95
2.38
2.76
3.13
3.48
];

Power_2micron_1=[0.515
0.973
1.42
1.825
2.21
2.56
2.9
3.22
];
Power_2micron_2=[0.545
1.02
1.478
1.91
2.33
2.715
3.05
3.37
];
Power_2micron_3=[0.525
1
1.44
1.86
2.28
2.65
2.97
3.31
];






Res_1micron_1=[24.843735
35.0475
51.335
67.3052
83.4904
98.9875
113.259
126.8694
];
Res_1micron_2=[17.02255
31.94492
45.79305
58.8766
73.1682
86.3403
99.4768
112.4133
];
Res_1micron_3=[18.62255
35.1256
51.27695
67.5605
83.9682
99.4564
113.4607
126.3972
];

Res_2micron_1=[25.07835
47.45347
68.2338
89.83425
113.6169
138.8534
157.656
186.1658
];
Res_2micron_2=[26.32255
50.8578
74.72042
102.4249
121.7487
144.12885
163.0895
181.5743
];
Res_2micron_3=[26.5975
50.89
74.3816
98.0054
121.6292
143.9335
165.2183
184.6709
];


Power_1micron_Mean=mean([Power_1micron_1 Power_1micron_2 Power_1micron_3],2);
Power_1micron_Error=std([Power_1micron_1 Power_1micron_2 Power_1micron_3],0,2);

Power_2micron_Mean=mean([Power_2micron_1 Power_2micron_2 Power_2micron_3],2);
Power_2micron_Error=std([Power_2micron_1 Power_2micron_2 Power_2micron_3],0,2);



Res_1micron_Mean=mean([Res_1micron_1 Res_1micron_2 Res_1micron_3],2);
Res_1micron_Error=std([Res_1micron_1 Res_1micron_2 Res_1micron_3],0,2);

Res_2micron_Mean=mean([Res_2micron_1 Res_2micron_2 Res_2micron_3],2);
Res_2micron_Error=std([Res_2micron_1 Res_2micron_2 Res_2micron_3],0,2);

figure('Position', [100, 100, 800, 500]);
errorbar(Green,Power_1micron_Mean,Power_1micron_Error,'-b','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
hold on 
errorbar(Green,Power_2micron_Mean,Power_2micron_Error,'-g','LineWidth',1.5);
legend('Side Coating = 1 micron','Side Coating = 2 micron','Location','NorthWest')
title('Output Ti:sapphire Power (mW @780nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)



figure('Position', [100, 100, 800, 500]);
errorbar(Green,Res_1micron_Mean,Res_1micron_Error,'-b','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Residual Pump Power (mW, 520nm)','fontsize',15);
hold on 
errorbar(Green,Res_2micron_Mean,Res_2micron_Error,'-g','LineWidth',1.5);
legend('Side Coating = 1 micron','Side Coating = 2 micron','Location','NorthWest')
title('Residual Pump Power (mW @520nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)