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

Power_A=[0.737
1.41
2.047
2.658
3.251
3.785
4.29
4.76
];

Power_B=[0.532
1.028
1.5
1.97
2.42
2.83
3.22
3.57
];


Res_A=[26.8838619
51.244767
75.2837589
99.3784246
123.6446337
145.6431795
166.107823
183.987412
];

Res_B=[19.2098284
37.4327436
55.44105
74.251639
92.174054
109.650121
126.238014
141.714559
];


figure('Position', [100, 100, 800, 500]);
plot(Green,Power_A,'-o',Green,Power_B,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('L2: 60x, A Wavelength Range','L2: 60x, B Wavelength Range','Location','NorthWest')
title('Output Ti:sapphire Power (mW @780nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)


figure('Position', [100, 100, 800, 500]);
plot(Green,Res_A,'-o',Green,Res_B,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('L2: 60x, A Wavelength Range','L2: 60x, B Wavelength Range','Location','NorthWest')
title('Residual Pump Power (mW @520nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)