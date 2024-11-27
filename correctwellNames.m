clearvars;
close all;

fcsPath = ['/home/ariana/Documents/CompExpRobot/EvolandComp_Bs210_Bs224/Competition_FCS/20240104_Bs210ctl_Bs211_4h']



% open dirs

fcsFiles=struct2cell(dir(fcsPath + "/*.fcs"));
fcsWellName = cellfun(@(x) string(x), fcsFiles(1, :), 'UniformOutput', false);
[fcsWellNum, sdx] = sort(cellfun(@(x) x{4}, cellfun(@(x) strsplit(x, '_'), fcsWellName, 'UniformOutput', false), 'UniformOutput', false));
fcsWellName = [fcsWellName{sdx} "Specimen_019_H12_H12_095.fcs"];




fcsWell = reshape(fcsWellName,[12, 8])';

fcsRot=rot90(rot90(fcsWell));
fcsRot2=reshape(fcsRot', [96, 1]);

fileNameArr(96, 3)= "zero";
 

[idx, where] = ismember([fcsWellName]', fcsRot2);
fileNameArr(:, 1) = fcsRot2;
fileNameArr(:, 2) = "new" + fcsWellName(where);
fileNameArr(:, 3) = [fcsWellName];


%%

for i=1:95

   
   movefile(fcsFiles{2, i} + "/" + fileNameArr(i, 3),  fcsFiles{2, i} + "/" + fileNameArr(i, 2));

end




