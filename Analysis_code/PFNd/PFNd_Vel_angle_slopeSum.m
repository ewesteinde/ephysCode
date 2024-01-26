function [] = PFNd_Vel_angle_slopeSum(date, cell, trial, p)
    varNames = {'Date', 'Cell', 'Trial', 'VmprefSlope','VmprefInt' 'VmoppSlope','VmoppInt','FRprefSlope','FRprefInt' 'FRoppSlope','FRoppInt'};
    trialTable = table(cellstr(date), cell, trial, p(1,1),p(1,2),p(2,1),p(2,2),p(1,3),p(1,4),p(2,3),p(2,4),'VariableNames',varNames);
    writetable(trialTable,'C:\Users\ewest\OneDrive\Documents\Data\PFNd_headingVelTrendSum.xlsx','WriteMode','Append','WriteVariableNames',false);
end