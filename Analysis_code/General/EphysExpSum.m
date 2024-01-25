function [] = EphysExpSum(trialMeta, date, cell, trial)
    varNames = {'Date','Experiment', 'Cell', 'Trial', 'Trial Type', 'Input R', 'Access R', 'Seal R', 'Pipette R', 'Notes'};
    trialTable = table(cellstr(date), cellstr(trialMeta.fly.flyExp), cell, trial, cellstr(trialMeta.trialType), trialMeta.inputR, trialMeta.accessR, trialMeta.sealR, trialMeta.pipetteR, cellstr(trialMeta.notes),'VariableNames',varNames);
    writetable(trialTable,'C:\Users\ewest\OneDrive\Documents\Data\ExpSum.xlsx','WriteMode','Append','WriteVariableNames',false);
end