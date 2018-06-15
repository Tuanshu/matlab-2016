clear all
cd('D:\AMO\160926_TiSa LS Development\161006_Data Figures');
Current=[0.1
0.2
0.25
0.3
0.35
0.4
0.45
0.5
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

Power_Old=[0.34
0.644
30.6
75.1
117
158
200
240
];

Power_New=[0.34
0.726
44.2
92.5
140
186
233
277
];


Res_Old=[26.8838619
51.244767
75.2837589
99.3784246
123.6446337
145.6431795
166.107823
183.987412
];

Res_New=[19.2098284
37.4327436
55.44105
74.251639
92.174054
109.650121
126.238014
141.714559
];


figure('Position', [100, 100, 800, 500]);
plot(Green,Power_Old,'-o',Green,Power_New,'-o','LineWidth',1.5);
xlabel('Pump Current (A)','fontsize',15);
ylabel('Green Laser Power (mW, 520nm)','fontsize',15);
legend('Homemade Heat Sink','Newport Heat Sink','Location','NorthWest')
title('Output Ti:sapphire Power (mW @780nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)


figure('Position', [100, 100, 800, 500]);
plot(Green,Res_A,'-o',Green,Res_B,'-o','LineWidth',1.5);
xlabel('Pump Power After L1 (mW, 520nm)','fontsize',15);
ylabel('Ti:sapphire ASE (mW, 780nm)','fontsize',15);
legend('L2: 60x, A Wavelength Range','L2: 60x, B Wavelength Range','Location','NorthWest')
title('Residual Pump Power (mW @520nm)','fontsize',15,'color','b')
set(gca,'fontsize',15)