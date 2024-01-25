%% Vm

figure(1);clf;

<<<<<<< HEAD
VmprefSlope = [0.358117441, 0.146685786;
0.341561492, 0.360022674;
0.092926175, -0.071015986;
0.704461468, 0.59405219;
0.671399216, 0.607275522;
0.389983334, 0.335456921;
0.217166885, 0.287753273;
0.349583424, 0.337883086;
0.732236645, 0.604441746;
0.780270134, 0.659184765;
1.054735802, 1.054094046;
0.95119543, 0.694086179;
0.07364525, -0.090828203;
-0.216655419, -0.419806067
]; 


VmprefSlope(:,1) = [-54.71964164
-37.97041722
-41.28883179
-45.59063601
-53.18154375
-61.31466889
-63.12366006
-65.21967153
-54.54449197
-56.096458
-53.78887705
-53.69407771
=======
VmprefSlope(:,1) = [-54.64029677
-38.57471691
-41.32621205
-45.19658243
-53.13956433
-61.05042303
-62.80031104
-65.17395361
-53.80294634
-55.7364531
-53.79788534
-53.914117
-48.72048818
-53.12903688

>>>>>>> 9b03edbff2b8ba294a89866cb656fbc3d452e101
];

VmprefSlope(:,2) = [-55.65812664
-42.09737336
-41.20272321
-45.35980804
-53.68109267
-62.65394451
-65.18871776
-65.38559355
-54.24380683
-56.08222067
-53.38284746
-59.68738423
-49.83257737
-54.10130719

];

[hVm, pVm] = ttest(VmprefSlope(:,1),VmprefSlope(:,2)); 
[pVm_wil, hVm_wil] = signrank(VmprefSlope(:,1),VmprefSlope(:,2));

x = [0.5 1.5];
y1 = VmprefSlope(1,:);
y2 = VmprefSlope(2,:);
y3 = VmprefSlope(3,:);
y4 = VmprefSlope(4,:);
y5 = VmprefSlope(5,:);
y6 = VmprefSlope(6,:);
y7 = VmprefSlope(7,:);
y8 = VmprefSlope(8,:);
y9 = VmprefSlope(9,:);
y10 = VmprefSlope(10,:);
y11 = VmprefSlope(11,:);
y12 = VmprefSlope(12,:);
<<<<<<< HEAD
y13 = VmprefSlope(13,:);
y14 = VmprefSlope(14,:);
=======
y13 = VmprefSlope(11,:);
y14 = VmprefSlope(12,:);
>>>>>>> 9b03edbff2b8ba294a89866cb656fbc3d452e101

meanPref = mean(VmprefSlope(:,1));
SEM_prefHead = std(VmprefSlope(:,1))./sqrt(length(VmprefSlope(:,1)));                               % Calculate Standard Error Of The Mean
CI95 = bsxfun(@plus, mean(VmprefSlope(:,1)), bsxfun(@times, [-1  1]*1.96, SEM_prefHead));   % 95% Confidence Intervals

meanOpp = mean(VmprefSlope(:,2));

figure(1); clf; 
plot(x,y1, '-ok','MarkerFaceColor','w')
hold on
plot(x,y2, '-ok','MarkerFaceColor','w')
plot(x,y3, '-ok','MarkerFaceColor','w')
plot(x,y4, '-ok','MarkerFaceColor','w')
plot(x,y5, '-ok','MarkerFaceColor','w')
plot(x,y6, '-ok','MarkerFaceColor','w')
plot(x,y7, '-ok','MarkerFaceColor','w')
plot(x,y8, '-ok','MarkerFaceColor','w')
plot(x,y9, '-ok','MarkerFaceColor','w')
plot(x,y10, '-ok','MarkerFaceColor','w')
plot(x,y11, '-ok','MarkerFaceColor','w')
plot(x,y12, '-ok','MarkerFaceColor','w')
plot(x,y13, '-ok','MarkerFaceColor','w')
plot(x,y14, '-ok','MarkerFaceColor','w')
plot(x, [meanPref, meanOpp],'-or','MarkerFaceColor','r');
ylabel('Vm (mV)')
xlim([0 2])
xNames = {'Pref head'; 'Anti-pref head'};
set(gca, 'xtick',[0.5 1.5], 'xticklabel', xNames)
<<<<<<< HEAD
ylabel('linear relationship b/w Vm & translational velocity')
%% FR

FRSlope(:,1) = [2.495649908
2.754142131
2.986221791
1.500160872
3.142462723
2.280675616
0.280685178
1.910691306
1.610940606
0.465339899
1.258078659
0.22459875
-1.01559856
];

FRSlope(:,2) = [1.088584279
2.156430777
1.707329914
1.188555525
2.868764062
1.286131721
0.101794243
1.367595558
0.96573278
0.125824435
0.949221508
-0.254836422
-0.881048197]

FRSlope(:,1) = [10.09048953
10.79203142
18.39034289
3.023171857
6.698409824
4.533820603
0.432528196
3.018536965
1.9920643
0.186513267
1.590878389
=======
set(gcf, 'color', 'w')

%% FR

FRSlope(:,1) = [9.444853925
10.66509511
18.4721368
2.496853587
5.960294282
4.273948645
0.318502339
2.736193798
1.162285185
0.276042678
1.506017063
11.53755068
8.35701775
>>>>>>> 9b03edbff2b8ba294a89866cb656fbc3d452e101

];

FRSlope(:,2) = [4.763466248
6.318999513
18.52583211
1.432827756
5.584488112
2.794225648
0.159824357
2.662847384
0.96764833
0.315971937
0.190896981
8.78627719
7.325139825
];

[hFR, pFR] = ttest(FRSlope(:,1),FRSlope(:,2)); 
[pFR_wil, hFR_wil] = signrank(FRSlope(:,1),FRSlope(:,2));

x = [0.5 1.5];
y1 = FRSlope(1,:);
y2 = FRSlope(2,:);
y3 = FRSlope(3,:);
y4 = FRSlope(4,:);
y5 = FRSlope(5,:);
y6 = FRSlope(6,:);
y7 = FRSlope(7,:);
y8 = FRSlope(8,:);
y9 = FRSlope(9,:);
y10 = FRSlope(10,:);
y11 = FRSlope(11,:);
y12 = FRSlope(12,:);
y13 = FRSlope(13,:);

meanPref = mean(FRSlope(:,1));
% SEM_prefHead = std(VmprefSlope(:,1))./sqrt(length(VmprefSlope(:,1)));                               % Calculate Standard Error Of The Mean
% CI95 = bsxfun(@plus, mean(VmprefSlope(:,1)), bsxfun(@times, [-1  1]*1.96, SEM_prefHead));   % 95% Confidence Intervals

meanOpp = mean(FRSlope(:,2));

figure(2); clf; 
plot(x,y1, '-ok','MarkerFaceColor','w')
hold on
plot(x,y2, '-ok','MarkerFaceColor','w')
plot(x,y3, '-ok','MarkerFaceColor','w')
plot(x,y4, '-ok','MarkerFaceColor','w')
plot(x,y5, '-ok','MarkerFaceColor','w')
plot(x,y6, '-ok','MarkerFaceColor','w')
plot(x,y7, '-ok','MarkerFaceColor','w')
plot(x,y8, '-ok','MarkerFaceColor','w')
plot(x,y9, '-ok','MarkerFaceColor','w')
plot(x,y10, '-ok','MarkerFaceColor','w')
plot(x,y11, '-ok','MarkerFaceColor','w')
plot(x,y12, '-ok','MarkerFaceColor','w')
plot(x,y13, '-ok','MarkerFaceColor','w')
plot(x, [meanPref, meanOpp],'-or','MarkerFaceColor','r');
ylabel('spikes/sec')
xlim([0 2])
xNames = {'Pref head'; 'Anti-pref head'};
set(gca, 'xtick',[0.5 1.5], 'xticklabel', xNames)
<<<<<<< HEAD
ylabel('linear relationship b/w firing rate & translational velocity')
=======
set(gcf, 'color', 'w')
>>>>>>> 9b03edbff2b8ba294a89866cb656fbc3d452e101
