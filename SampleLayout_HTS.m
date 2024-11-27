clearvars;
close all;

fcsPath = ['/home/ariana/Documents/CompExpRobot/EvolandComp_Bs210_Bs224/Competition_FCS/20240321_210BC_Bs211_4h']
layoutPath = '/home/ariana/Documents/CompExpRobot/EvolandComp_Bs210_Bs224/Layout/Layout_20240321_Bs210BC.txt'

%%import layout 

fid = fopen(layoutPath);
f1 = textscan(fid, '%s %s %*c', 'delimiter', ',', 'HeaderLines', 1);
well = cellfun(@(x) string(x), f1{1});
sample= cellfun(@(x) string(x), f1{2});
fclose(fid);

% open dirs

fcsFiles=struct2cell(dir(fcsPath + "/*.fcs"));
fcsWellName = cellfun(@(x) string(x), fcsFiles(1, :), 'UniformOutput', false);
fcsWell = cellfun(@(x) x{3}, cellfun(@(x) strsplit(x, '_'), fcsWellName, 'UniformOutput', false), 'UniformOutput', false);
fcsWell = cellfun(@(x) string(x), fcsWell);





fileNameArr(1:length(fcsFiles), 1:3)= "zero";

[idx, where] = ismember(fcsWell, well);

fileNameArr(:, 1) = fcsWellName(:) ;
fileNameArr(:, 2) = well(where);
fileNameArr(:, 3) = sample(where);


%%

for i=1:length(fcsFiles)
    
    movefile(fcsFiles{2, i} + "/" + string(fcsFiles{1, i}),  fcsFiles{2, i} + "/" + fileNameArr(i, 3) + ".fcs");

end




