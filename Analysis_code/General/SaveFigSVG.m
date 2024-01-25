rootDir = 'Z:\Dropbox (HMS)\Wilson_Lab_Data\Analysis\Westeinde2022_ms\revision_plots\ephys_plots\vm_notAbs';

folders = dir(rootDir);

for f = 1:length(folders)
    if regexp(folders(f).name,'.fig')
        figName = fullfile(folders(f).folder,folders(f).name); 
        openfig(figName)
        folderBase = folders(f).name(1:end-4); 
        
        saveas(gcf,fullfile(folders(f).folder,[folderBase,'.svg']));
    end
end