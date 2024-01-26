function newTable = convertStruct2Table(struct, save, rootPath)

    f = fieldnames(struct);
    for k=1:numel(f)
        if size(struct.(f{k}), 2) > size(struct.(f{k}), 1)
            struct.(f{k}) = struct.(f{k})'; 
        end
    end
    
    newTable = struct2table(struct);
    
    if strcmp(save, 'y')
        fileName = input('File name? include extension: ','s'); 
        pathName = fullfile(rootPath, fileName); 
        writetable(newTable, pathName)
    end
end
        
        