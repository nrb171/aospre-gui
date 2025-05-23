app.FileList = dir('/glade/u/home/nbarron/scratch/wrf/wrfv4.6.1/run/stitched/');

app.FileList = app.FileList(~ismember({app.FileList.name},{'.','..'}));
app.FileList = app.FileList(contains({app.FileList.name},'wrfout'));

variablesToKeep = ["W", "QRAIN"];

requiredVariables = ["XLAT", "XLONG", "P", "PB", "PH", "PHB"];

%% check existence of folder containing abbreviated files
abbreviatedFolder = [app.FileList(1).folder, '/abbreviated/'];
if ~exist(abbreviatedFolder, 'dir')
    mkdir(abbreviatedFolder)
end

%% create variable slice at 2 km
for iFile = 1:12:numel(app.FileList)
    nci = ncinfo([app.FileList(iFile).folder, '/', app.FileList(iFile).name]);
    allVarNames = string({nci.Variables.Name});
    for varName = [requiredVariables, variablesToKeep]
        varData = ncread([app.FileList(iFile).folder, '/', app.FileList(iFile).name], varName);
        varIndx = find(allVarNames == varName);
        varDimNames = string({nci.Variables(varIndx).Dimensions.Name});
        varDimLength = [nci.Variables(varIndx).Dimensions.Length];
        varDimCell = {};
        for iDim = 1:length(varDimLength)
            if ~contains(varDimNames(iDim), 'bottom_top')
                varDimCell{end+1} = varDimNames(iDim);
                varDimCell{end+1} = varDimLength(iDim);
            else
                varDimCell{end+1} = varDimNames(iDim);
                varDimCell{end+1} = 2;
            end
        end
        nccreate([abbreviatedFolder,app.FileList(iFile).name], varName, "Dimensions",varDimCell)
        if numel(varDimCell) == 8
            ncwrite([abbreviatedFolder,app.FileList(iFile).name], varName, varData(:,:,28:29,:))
        else
            ncwrite([abbreviatedFolder,app.FileList(iFile).name], varName, varData)
        end 
    end
end