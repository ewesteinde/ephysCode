%folders = get_folders_ephys(rootDir);

summaryArray = summaryArray(~ismembertol(summaryArray.rho,1,10^-10),:); % gets rid of trials w/ no heading change --> indicates problem
summaryArray = summaryArray(summaryArray.timeMov > 5,:); % only look at trials where fly's vel was above threshold for at least 5 seconds

%     if meno
%         summaryArray = summaryArray(summaryArray.rho > 0.5,:);
%     else
%         summaryArray = summaryArray(summaryArray.rho < 0.5,:);
%     end

savePlots = 1; 
count = 1; 
for f = 1:size(summaryArray,1)
   % try
        folder = summaryArray.Folder(f); 
            if strcmp(folder(end),'.')
                folder = folder(1:end-2); 
            end

            processedDir = fullfile(folder,'processedData');
            load(fullfile(processedDir,'pro_behaviourData.mat'))
            load(fullfile(processedDir,'pro_trialData.mat'))

        t = summaryArray.numTrial(f);
        bData = pro_behaviourData{t};
        tData = pro_trialData{t};
        bData = bData(summaryArray.Indices{f},:);
        tData = tData(summaryArray.Indices{f},:);
        [vecStr, goal] = CalculateAverageHeading_ephys(bData,1.5, 'all');
        prefHead = calcPrefHead(bData, tData);
        %[vf_coeff, vy_coeff] = PFL2_3_behaviourVSactivity_lineplots_new(prefHead, 0.5, 120, tData, bData) ;
        [vfCoeff,vyCoeff] = PFL2_3_LinePlotsFitLine_ephys(folder, bData, tData, 0);

        %saveas(gcf,fullfile(folder,'figures','activityVsbehaviour_sepHeadings.fig'))
        
        isRight = regexp(folder,'RPFL3');
        isPFL2 = regexp(folder,'PFL2');
        if isempty(isRight)
            Summary.RPFL3(count,1) = 0; 
        else
            Summary.RPFL3(count,1) = 1; 
        end
        
        if isempty(isPFL2)
            Summary.PFL2(count,1) = 0;
        else
            Summary.PFL2(count,1) = 1;
        end
        
        Summary.folder(count,1) = string(folder);
        Summary.trial(count,1) = t;
        Summary.prefHead(count,1) = prefHead;
        Summary.goal(count,1) = goal;
        Summary.vecStr(count,1) = vecStr;
        Summary.vfPrefSlope(count,1) = vfCoeff(1);
        Summary.vfPrefyInt(count,1) = vfCoeff(2);
        %Summary.vfantiPrefSlope(count,1) = vf_coeff(2,1);
        %Summary.vfantiPrefyInt(count,1) = vf_coeff(2,2);
        Summary.vyPrefSlope(count,1) = vyCoeff(1);
        Summary.vyPrefyInt(count,1) = vyCoeff(2);
        %Summary.vyantiPrefSlope(count,1) = vy_coeff(2,1);
        %Summary.vyantiPrefyInt(count,1) = vy_coeff(2,2);

        count = count + 1;                  
end
close all 
summaryTable = struct2table(Summary');
%%
goaldiff = rad2deg(angdiff(-deg2rad(summaryTable.prefHead),-summaryTable.goal)); 

goaldiff = goaldiff;%(summaryTable.vecStr > 0.5);
tempTable = summaryTable;%(summaryTable.vecStr > 0.5,:);

figure()
set(gcf,'color','w')
set(gcf,'renderer','painters')
hold on
scatter(goaldiff(tempTable.RPFL3 == 1),tempTable.vyPrefSlope(tempTable.RPFL3 == 1),[],'r')
scatter(goaldiff(tempTable.RPFL3 == 0 & tempTable.PFL2 == 0),tempTable.vyPrefSlope(tempTable.RPFL3 == 0 & tempTable.PFL2 == 0),[],'b')
scatter(goaldiff(tempTable.PFL2 == 1),tempTable.vyPrefSlope(tempTable.PFL2 == 1),[],'k')
[vycoeff, vyparams] = polyfit(goaldiff, tempTable.vyPrefSlope, 3);
xFit = linspace(min(goaldiff), max(goaldiff), 1000);
yFit = polyval(vycoeff, xFit);
plot(xFit,yFit)
ylabel('slope of yaw relationship')
xlabel('goal - prefHead')
xlim([-180,180])
legend({'RPFL3','LPFL3','PFL2'})

figure()
set(gcf,'color','w')
set(gcf,'renderer','painters')
hold on
scatter(goaldiff(tempTable.RPFL3 == 1),tempTable.vyPrefyInt(tempTable.RPFL3 == 1),[],'r')
scatter(goaldiff(tempTable.RPFL3 == 0 & tempTable.PFL2 == 0),tempTable.vyPrefyInt(tempTable.RPFL3 == 0 & tempTable.PFL2 == 0),[],'b')
%scatter(goaldiff(tempTable.PFL2 == 1),tempTable.vyPrefyInt(tempTable.PFL2 == 1),[],'k')

[vycoeff, vyparams] = polyfit(goaldiff(tempTable.RPFL3 == 1),tempTable.vyPrefyInt(tempTable.RPFL3 == 1), 3);
xFit = linspace(min(goaldiff(tempTable.RPFL3 == 1)), max(goaldiff(tempTable.RPFL3 == 1)), 1000);
yFit = polyval(vycoeff, xFit);
plot(xFit,yFit,'r')

[vycoeff, vyparams] = polyfit(goaldiff(tempTable.RPFL3 == 0 & tempTable.PFL2 == 0), tempTable.vyPrefyInt(tempTable.RPFL3 == 0 & tempTable.PFL2 == 0), 3);
xFit = linspace(min(goaldiff(tempTable.RPFL3 == 0 & tempTable.PFL2 == 0)), max(goaldiff(tempTable.RPFL3 == 0 & tempTable.PFL2 == 0)), 1000);
yFit = polyval(vycoeff, xFit);
plot(xFit,yFit,'b')

% [vycoeff, vyparams] = polyfit(goaldiff(tempTable.PFL2 == 1), tempTable.vyPrefyInt(tempTable.PFL2 == 1), 4);
% xFit = linspace(min(goaldiff(tempTable.PFL2 == 1)), max(goaldiff(tempTable.PFL2 == 1)), 1000);
% yFit = polyval(vycoeff, xFit);
% plot(xFit,yFit,'k')

ylabel('yInt of yaw relationship')
xlabel('goal - prefHead')
xlim([-180,180])
legend({'RPFL3','LPFL3'})

figure()
set(gcf,'color','w')
set(gcf,'renderer','painters')
hold on
scatter(goaldiff,tempTable.vfPrefSlope,[],tempTable.vecStr)
vfcoeff = polyfit(goaldiff, tempTable.vfPrefSlope, 3);
xFit = linspace(min(goaldiff), max(goaldiff), 1000);
yFit = polyval(vfcoeff, xFit);
plot(xFit,yFit)
ylabel('slope of + yaw rel')
xlabel('goal - prefHead')
xlim([-180,180])
%%
figure()
for idx = 1:size(tempTable,1)
    plot([1,1.5],[tempTable.vfPrefyInt(idx),tempTable.vfantiPrefyInt(idx)],'o-')
    hold on
end
xlim([0.5,2])

figure()
for idx = 1:size(tempTable,1)
    plot([1,1.5],[tempTable.vyPrefSlope(idx),tempTable.vyantiPrefSlope(idx)],'o-')
    hold on
end
xlim([0.5,2])

figure()
for idx = 1:size(tempTable,1)
    plot([1,1.5],[tempTable.vyPrefyInt(idx),tempTable.vyantiPrefyInt(idx)],'o-')
    hold on
end
xlim([0.5,2])


