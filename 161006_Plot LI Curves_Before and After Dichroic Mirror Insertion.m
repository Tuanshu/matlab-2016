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

Power_Before=[0.548
1.051
1.536
2.029
2.509
2.95
3.365
3.776
];

Power_After=[0.49
0.956
1.4
1.82
2.24
2.61
2.96
3.3
];


Res_Before=[19.3085676
37.5164937
54.8319632
73.8593023
92.3824783
109.393165
125.7945255
141.8054512
];

Res_After=[16.905763
32.6659172
47.67518
62.203834
77.397488
91.634207
104.450752
116.62371
];


figure('Position', [100, 100, 800, 500]);
plot(Green,Power_Before,'-o',Green,Power_After,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('Before Dichroic Mirror Insertion','After Dichroic Mirror Insertion','Location','NorthWest')
title('Output Ti:sapphire Power (mW @780nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)


figure('Position', [100, 100, 800, 500]);
plot(Green,Res_Before,'-o',Green,Res_After,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('Before Dichroic Mirror Insertion','After Dichroic Mirror Insertion','Location','NorthWest')
title('Residual Pump Power (mW @520nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)