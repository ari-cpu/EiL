% % %% %% %% %% FACS_Plot.m %% %% %% %% % %

% This scripts reads in your measured competition experiment fractions,
% calculates the selection coeffients.

% You can choose to correct for
% 1. the fraction of non-fluorescent events in the pure GFP signal,
% ---------------- corrNonFluorescentGFP
% 2. the fitness difference between ancestor and GFP reporter strain
% ---------------- corrANCvsGFP
%
% When reading in the results of more than one experiment, the different
% experiments will be compared.
%
%Functions needed:
% ----> scoeff -- gives the selection coefficient and a error based on the
%       uncertainty of the cytometer measurement of the fractionsuuu(203)
% ----> getMeanError -- gives you the error through error propagation of a
%       mean
%
% ------ Outputs ------
%
% **outSIndividual**
% contains the start/end fraction and the calculated selection coefficient
% for very single measured sample
%
% **outS**
% contains mean(selCoeff) and std(selCoeff) of a sample for each run
% --> the mean equals the points you get from a standard competition exp.
%
% **outCollect**
% collects the mean of the means per run
%
% ----------------------
% % %% %% %% %% %% %% %% %% %% %% %% %% % %


clearvars
close all


%% define paths and inputs

basePath = "/home/ariana/Documents/CompExpRobot/EvolandComp_Bs210_Bs224/";


% new experiments June 2023
%facsdata(1) = basePath + "out/" + "20240124_Bs224hyb_A.txt";
%facsdata(2) = basePath + "out/" + "20240124_Bs224hyb_B.txt";
%facsdata(3) = basePath + "out/" + "20240125_Bs224hyb_A.txt";
%facsdata(4) = basePath + "out/" + "20240125_Bs224hyb_B.txt";

%facsdata(1) = basePath + "out/Filaments_corr/" + "20240124_Bs224hyb_A_Corr4.txt";
%facsdata(2) = basePath + "out/Filaments_corr/" + "20240124_Bs224hyb_B_Corr4.txt";
%facsdata(3) = basePath + "out/Filaments_corr/" + "20240125_Bs224hyb_A_Corr4.txt";
%facsdata(4) = basePath + "out/Filaments_corr/" + "20240125_Bs224hyb_B_Corr4.txt";

%facsdata(1) = basePath + "out/" + "20241601_Bs224ctl_A.txt";
%facsdata(2) = basePath + "out/" + "20241601_Bs224ctl_B.txt";
%facsdata(3) = basePath + "out/" + "20241701_Bs224ctl_A.txt";
%facsdata(4) = basePath + "out/" + "20241701_Bs224ctl_B.txt";


%facsdata(1) = basePath + "out/Filaments_corr/" + "20241601_Bs224ctl_A_Corr4.txt";
%facsdata(2) = basePath + "out/Filaments_corr/" + "20241601_Bs224ctl_B_Corr4.txt";
%facsdata(3) = basePath + "out/Filaments_corr/" + "20241701_Bs224ctl_A_Corr4.txt";
%facsdata(4) = basePath + "out/Filaments_corr/" + "20241701_Bs224ctl_B_Corr4.txt";

%facsdata(1) = basePath + "out/" + "20241001_Bs210Hyb_A.txt";
%facsdata(2) = basePath + "out/" + "20241001_Bs210Hyb_B.txt";
%facsdata(3) = basePath + "out/" + "20241101_Bs210Hyb_A.txt";
%facsdata(4) = basePath + "out/" + "20241101_Bs210Hyb_B.txt";

%facsdata(1) = basePath + "out/Filaments_corr/corr4/" + "20241001_Bs210hyb_A_Corr4.txt";
%facsdata(2) = basePath + "out/Filaments_corr/corr4/" + "20241001_Bs210hyb_B_Corr4.txt";
%facsdata(3) = basePath + "out/Filaments_corr/corr4/" + "20241101_Bs210hyb_A_Corr4.txt";
%facsdata(4) = basePath + "out/Filaments_corr/corr4/" + "20241101_Bs210hyb_B_Corr4.txt";

%facsdata(1) = basePath + "out/Filaments_corr/MonaScript_out/" + "20241001A_Bs210hyb_woGaps_filamentCorr_v4.txt";
%facsdata(2) = basePath + "out/Filaments_corr/MonaScript_out/" + "20241001B_Bs210hyb_woGaps_filamentCorr_v4.txt";
%facsdata(3) = basePath + "out/Filaments_corr/MonaScript_out/" + "20241101A_Bs210hyb_woGaps_filamentCorr_v4.txt";
%facsdata(4) = basePath + "out/Filaments_corr/MonaScript_out/" + "20241101B_Bs210hyb_woGaps_filamentCorr_v4.txt";



%facsdata(1) = basePath + "out/" + "20240401_Bs210Ctl.txt";
%facsdata(2) = basePath + "out/" + "20240501_Bs210Ctl_A.txt";
%facsdata(3) = basePath + "out/" + "20240501_Bs210Ctl_B.txt";
%facsdata(4) = basePath + "out/" + "20240131_Bs210ctl_A.txt";
%facsdata(5) = basePath + "out/" + "20240131_Bs210ctl_B.txt";
%facsdata(6) = basePath + "out/" + "20231220_Bs210Ctl_001.txt";

%facsdata(1) = basePath + "out/Filaments_corr/corr4/" + "20240104_Bs210ctl__Corr4.txt";
%facsdata(2) = basePath + "out/Filaments_corr/corr4/" + "20240105_Bs210ctl_A_Corr4.txt";
%facsdata(3) = basePath + "out/Filaments_corr/corr4/" + "20240105_Bs210ctl_B_Corr4.txt";
%facsdata(4) = basePath + "out/Filaments_corr/corr4/" + "20240131_Bs210ctl_A_Corr4.txt";
%facsdata(5) = basePath + "out/Filaments_corr/corr4/" + "20240131_Bs210ctl_B_Corr4.txt";
%facsdata(6) = basePath + "out/Filaments_corr/corr4/" + "20231220_Bs210ctl_001_Corr4.txt";

%facsdata(1) = basePath + "out/Filaments_corr/MonaScript_out/" + "20240104_Bs210ctl_woGaps_filamentCorr_v4.txt";
%facsdata(2) = basePath + "out/Filaments_corr/MonaScript_out/" + "20240105A_Bs210ctl_woGaps_filamentCorr_v4.txt";
%facsdata(3) = basePath + "out/Filaments_corr/MonaScript_out/" + "20240105B_Bs210ctl_woGaps_filamentCorr_v4.txt";
%facsdata(4) = basePath + "out/Filaments_corr/MonaScript_out/" + "20240131A_Bs210ctl_woGaps_filamentCorr_v4.txt";
%facsdata(5) = basePath + "out/Filaments_corr/MonaScript_out/" + "20240131B_Bs210ctl_woGaps_filamentCorr_v4.txt";
%facsdata(6) = basePath + "out/Filaments_corr/MonaScript_out/" + "20231220001_Bs210ctl_woGaps_filamentCorr_v4.txt";


%facsdata(1) = basePath + "out/" + "20241201_210ANC_A.txt";
%facsdata(2) = basePath + "out/" + "20241201_210ANC_B.txt";
%facsdata(3) = basePath + "out/" + "20241501_210ANC_A.txt";
%facsdata(4) = basePath + "out/" + "20241501_210ANC_B.txt";

facsdata(1) = basePath + "out/Filaments_corr/corr4/" + "20241201_210ANC_A_Corr4.txt";
facsdata(2) = basePath + "out/Filaments_corr/corr4/" + "20241201_210ANC_B_Corr4.txt";
facsdata(3) = basePath + "out/Filaments_corr/corr4/" + "20241501_210ANC_A_Corr4.txt";
facsdata(4) = basePath + "out/Filaments_corr/corr4/" + "20241501_210ANC_B_Corr4.txt";

%facsdata(1) = basePath + "out/Filaments_corr/corr4/" + "20240319_Bs210BC_A_Corr4.txt";
%facsdata(2) = basePath + "out/Filaments_corr/corr4/" + "20240319_Bs210BC_B_Corr4.txt";
%facsdata(3) = basePath + "out/Filaments_corr/corr4/" + "20240321_210BC__Corr4.txt";

ExpOutName="Bs210ANC"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Do you want to compare with the ancestor?
compareAnc = "OFF";
% .. then I need to know which outCollect.mat of the ancestor to use ..
%AncDist = "/home/isabel/Dokumente/P2_FitnessExp/1_CompExp_Robot/2021_AncBs166/20210423_AncRuns_1to3_woLag_outCollect_allPlates.mat";
AncDist = "/media/isabel/IsabelFestplatte/Doktor/P2_FitnessExp/1_CompExp/20210628_outCollect_Bs166_woLag";

%%%%%%% what do you want to plot?
plotSingleRuns          = "OFF";  % default: "ON" 
plotComparingSingleRuns = "OFF";
plotSystematics         = "OFF";
plotErrorPlots          = "OFF";
BCcond = "ON";
printFigs      = "OFF";

corrNonFluorescentGFP = "ON";   % You can switch "OFF" and "ON"
corrANCvsGFP          = "ON";   % You can switch "OFF" and "ON"

samplePrefix = "BVAL";% LibGRBval "LibEvolExp"; % or LibSCW23 or Bs166 or LibSCBval or LibSCBmoj
%ExpOutName   = "Bs210Ctl_noCorr";

% outPath for plots
savePlotPath = "/home/ariana/Documents/CompExpRobot/EvolandComp_Bs210_Bs224/plots/FilamentCorrection/";


% % % Here you have to specify some options


% Here you have the option to exclude samples from the merged plots (100,
% 101, ...) that contain a certain string, e.g. "Bs166"
excludeType = "ismember"; %"contains" ; % "ismember"
% excludeSample = ["Bs19-" + string([12 24 36 48 60])];
% excludeSample = ["Bs19-" + string([12 24 36 48 60]), "Bval" + ["0" + string(1:9), string([10:19,21:30])], "Bval-"];
excludeSample = [];% "Bs166" + string([12 24 36 48 60]); 
%["LibBvalwS3ctld5"] + string([2:6,8:26,28:34,36,37,39:44,46,48,50:56,58:73,75,76,78:82,84,86,87,88,90:96]);%["Bs166" + string([12,24,36,48,60])];

%["Bans" + ["0" + string([1:9]), string(10:20)]] + string(20),...


% OTHER OPTION - ANC DISTR: If you want to correct each day with another set of anc
% (for ANC Distr!!!) and exclude other samples as the correction ones ...
% excludeSampleopt(1,:) = ["Bs166" + string([12 24 36 48 60] - 2)  "Bs166" + string(85:96)]; 
% excludeSampleopt(2,:) = ["Bs166" + string([12 24 36 48 60] - 5)  "Bs166" + string(85:96)];
% excludeSampleopt(3,:) = ["Bs166" + string([12 24 36 48 60] - 11) "Bs166" + string(85:96)];
% excludeSampleopt(4,:) = ["Bs166" + string([12 24 36 48 60] - 4)  "Bs166" + string(85:96)];

% how are your ancestor and gfp samples called -- needed for plot legends
% ancestor = "Bs19";
% gfp      = "Bs205";
ancestor = "Bs210";
gfp      = "Bs211"
%BCcorr = -0.017;  %Bs210
%BCcorr = -0.014;   %Bs224
BCcorr = 0; %ANC

% with which samples do you want to correct in corrANCvsGFP?
% these will also be excluded in all distributions
% corrANCsamples = "Bs19-" + string([12 24 36 48 60]);
corrANCsamples = "Bs210" + ["V" "W" "X" "Y" "Z"];

% OTHER OPTION - ANC DISTR: If you want to correct each day with another set of anc
% (for ANC Distr!!!)
% corrANCopt(1,:) = "Bs166" + string([12 24 36 48 60] -2); 
% corrANCopt(2,:) = "Bs166" + string([12 24 36 48 60] - 5);
% corrANCopt(3,:) = "Bs166" + string([12 24 36 48 60] - 11);
% corrANCopt(4,:) = "Bs166" + string([12 24 36 48 60] - 4);

minFracAnc = 30;
maxFracAnc = 70;

%Bs210 16.1
%Bs224 16.2

genTime = 16.1/60; % for CM medium it is 1 h; for MM it is 175/60

% what is the uncertainty of the fraction measurement?
uncertain = 5;

% Do you want to exclude samples with a start fraction smaller than
minFrac = 40;
% or greater than:
maxFrac = 60;
% from the selection coefficient calculations?

% All selection coefficients that come from less than dataPointExcl data
% points (runs), will be excluded in the outCollect
dataPointExcl = 2;

% plan the histograms
hists = []; histLegends = []; histmax = []; histBinWidth = 0.005; 
histBinEdges = -40*histBinWidth-0.5*histBinWidth:histBinWidth:40*histBinWidth+0.5*histBinWidth;
cmp = repmat([0 76 255; 0 0 0; 189 38 38; 81 189 81; 235 235 84; 130 218 250; 247 134 247]/255, 70,1);
cmp_runs = [3 136 252; 88 166 111; 254 179 67; 250 69 10; 150 150 150; 120 240 100]/255; 

% % which samples did you expect to measure? the missing ones will be saved
% to samplesLost
upperLim = 96;
notIn = 12*[1:8];
nam = repmat(samplePrefix,1,upperLim) + string(1:upperLim);
names_all = nam(~ismember(1:upperLim,notIn));

%% Run through all sets of data seperately (merging begins in line ~126)

for c = 1 : length(facsdata)

    clear data outc meanS stdS time_* t
     
    % set variables 
    if exist("corrANCopt", "var")
        % Check correct length of corrANCopt
        if size(corrANCopt,1)~= size(facsdata,2)
            error("Specify your corrANCopt samples for each input or delete the variable and use corrANCsamples instead!"); 
        end
        corrANCsamples = corrANCopt(c,:);
    end
    
     if exist("excludeSampleopt", "var")
        % Check correct length of corrANCopt
        if size(excludeSampleopt,1)~= size(facsdata,2)
            error("Specify your excludeSampleopt samples for each input or delete the variable and use corrANCsamples instead!"); 
        end
        excludeSample = excludeSampleopt(c,:);
     end
    
    % read data
    fid = fopen(facsdata(c));
    fgetl(fid); % Skip first headerline
    % Read in second headerline and extract t(1) and t(2):
    columnNames = fgetl(fid); columnNamesSep = strsplit(columnNames, {'\tnonGfpStart_', 'h\tnonGfpEnd_', 'h\ttimestampStart\t'});
    t(1) = str2num(columnNamesSep{2}); t(2) = str2num(columnNamesSep{3}); 
    
    imp = textscan(fid,'%s %f %f %s %s','delimiter','\t');
    fclose(fid);
    data.sample=imp{1}; data.start=imp{2}; data.ende=imp{3}; data.startTime = imp{4}; data.endeTime = imp{5};
    clear imp
    
    % get names and times
    nameSplit = cellfun(@(x) strsplit(x, '_'), [data.sample], 'UniformOutput', false);
    name_first = cellfun(@(x) x{1}, nameSplit, 'UniformOutput', false);
    name_last = cellfun(@(x) x{2}, nameSplit, 'UniformOutput', false);
    %plate = cellfun(@(x) x{end}, nameSplit, 'UniformOutput', false);
    
    %plateSplit = cellfun(@(x) strsplit(x,'plate'), plate, 'UniformOutput', false);
    %plateNum = cellfun(@(x) str2num(x{2}), plateSplit);
    
    timeSplit_start = cellfun(@(x) str2double(strsplit(x,':')),data.startTime,'UniformOutput',false);
    timeSplit_ende  = cellfun(@(x) str2double(strsplit(x,':')),data.endeTime,'UniformOutput',false);
    
    for i=1:length(data.sample)
        time_start(i) = timeSplit_start{i}(1)*3600 + timeSplit_start{i}(2)*60 + timeSplit_start{i}(3);
        time_ende(i)  = timeSplit_ende{i}(1)*3600  + timeSplit_ende{i}(2)*60  + timeSplit_ende{i}(3);
    end
    
    % Find gfps
    pureGFP{c} = find(cellfun('isempty', name_first) & cellfun(@(x) contains(x,gfp), name_last))';
    s_woGFP{c} = find(~cellfun('isempty', name_first));
    
     %% Correct raw data - 1 - non-fluorenscent GFP correction:
    
    switch corrNonFluorescentGFP
        
        case "ON"
            % percentage of counts that are nongfp if we measure a 100% gfp sample:
            if ~isempty(pureGFP{c})
                nongfp_start = mean(data.start(pureGFP{c})./(100-data.start(pureGFP{c})));
                nongfp_ende  = mean( data.ende(pureGFP{c})./(100- data.ende(pureGFP{c})));
            else
                nongfp_start = 0; nongfp_ende = 0 ;
            end
            
            % outcommented you see the old correction
            start_corr = (data.start(:) - (100-data.start(:))* nongfp_start);
            data.start_corr = start_corr;
            ende_corr =  (data.ende(:)  - (100-data.ende(:)) * nongfp_ende);
            data.ende_corr = ende_corr;
            
        case "OFF"
            data.start_corr = data.start;
            data.ende_corr = data.ende;
            
        otherwise
            error("Your correction 1 - non-fluorenscent GFP correction - needs to be turned ON or OFF.");
            
    end
    
    %% Calculate selection coefficients
    % in outc, only the samples appear, that are not pure gfp ..
    fprintf("Calculating selection coefficients (t = [" + num2str(t(1))+ " "+ num2str(t(2))+ "])... \n");
    
    for i = 1 : length(s_woGFP{c})
        idx_woGFP = s_woGFP{c}(i);
        outc(i).sample = string(data.sample(idx_woGFP));
        outc(i).start_s = data.start_corr(idx_woGFP);
        outc(i).ende_s = data.ende_corr(idx_woGFP);
        
        [outc(i).s,outc(i).s_Err] = scoeff(outc(i).start_s,outc(i).ende_s,uncertain,genTime,t(2)-t(1));
        outc(i).s_ErrPerc = 100 * outc(i).s_Err/outc(i).s;
        
        outc(i).sampleID = string(name_first{idx_woGFP});
        outc(i).plate = 1; %plateNum(idx_woGFP);
        outc(i).tendency = "neutral";
        
        outc(i).startTime = time_start(idx_woGFP);
        outc(i).endeTime  = time_ende(idx_woGFP);
    end
    
    % TENDENCY from plate to plate 
    % we want to test if for all experiment days the tendency is correct,
    % meaning, if the amount of gfp goes up from plate 1 to 4
%     [namesHere,~]     = unique([outc(:).sampleID]);
%     
%     for i=1:numel(namesHere)
%         idx = find(ismember([outc(:).sampleID],namesHere(i)));
%         starts_test = [outc(idx).start_s];
%         plates_test = [outc(idx).plate];
%         [a,~,where] = unique(plates_test);
%         if numel(a)==1
%             outc(idx).tendency = deal("neutral");
%         elseif numel(a)==2 && mean(starts_test(where==2)) - mean(starts_test(where==1)) < 4
%             [outc(idx).tendency] = deal("good");
%             
%         elseif numel(a)==3 && mean(starts_test(where==2)) - mean(starts_test(where==1)) < 4 ...
%                 && mean(starts_test(where==3)) - mean(starts_test(where==2)) < 4
%             
%             [outc(idx).tendency] = deal("good");
%         elseif numel(a)==4 && mean(starts_test(where==2)) - mean(starts_test(where==1)) < 4 ...
%                 && mean(starts_test(where==3)) - mean(starts_test(where==2)) < 4 ...
%                 && mean(starts_test(where==4)) - mean(starts_test(where==3)) < 4
%             
%             [outc(idx).tendency] = deal("good");
%         else
%             [outc(idx).tendency] = deal("bad");
%         end
%         
%     end
%     
%     countTendency{c} = sum([outc.tendency]=="bad")/numel([outc.tendency]);
    
    clear numbers sortOrder
    %% Correct raw data - 2 - Correction of difference btw reporter and anc PER PLATE:
    
    % Create mask for excluding all ancestor entries with start fractions smaller than
    % minFrac or greater than maxFrac
    fracMaskcorrAnc = ~([outc.start_s]<minFracAnc | [outc.start_s]>maxFracAnc);
    
    switch corrANCvsGFP
        case "ON"
            
            for pl=1:5
                outc_plateMask = [];
                outc_plateMask = [outc.plate]==pl;
                outc_plateIdx = find([outc.plate]==pl);
                if sum(outc_plateMask) > 0
                    
                    % Find ancestors with which to correct
                    ancs = ismember([outc.sampleID], corrANCsamples) & outc_plateMask;
                    fprintf("Correcting plate " + num2str(pl) + " with " + sum(ancs & fracMaskcorrAnc) + ".\n");
                    
                    if sum(ancs & fracMaskcorrAnc) == 0
                        fprintf("There are no ancestors to correct plate " + num2str(pl) + ", \nso all data points from plate " + num2str(pl) + " will be excluded in the following anaylsis!\n")
                        for i = 1 : length(outc_plateIdx)
                            tmp = outc_plateIdx(i);
                            outc(tmp).s_corr = NaN;
                            outc(tmp).s_corrErr = NaN;
                        end
                        continue
                    end
                    
                    [anc_Mean,anc_MeanErr] = getMeanError([outc(ancs & fracMaskcorrAnc).s],[outc(ancs&fracMaskcorrAnc).s_Err]);
                    
                    
                    for i = 1 : length(outc_plateIdx)
                        tmp = outc_plateIdx(i);
                        outc(tmp).s_corr = outc(tmp).s - anc_Mean;
                        
                        outc(tmp).s_corrErr = sqrt( outc(tmp).s_Err ^2 + anc_MeanErr ^2) / 2;
                        
                    end
                    
                end
            end
            
            
        case "OFF"
            for i = 1 : length(s_woGFP{c})
                outc(i).s_corr = outc(i).s;
                outc(i).s_corrErr = outc(i).s_Err;
            end
            
        otherwise
            error("Your ANC vs GFP - correction is not working properly!");
    end
    
    
    %% Exclude samples:
    % % -- from excludeSample
    % % -- the ones that cannot be corrected (no anc on the plate)
    
    % % % % % % % % % % % % % % % % % 
    % % % % A mask for exclude is made
     exclMask = zeros(1,length([outc.sample]));
    if ~isempty(excludeSample) && excludeType == "ismember"
        exclMask = ismember([outc.sampleID],excludeSample);
    elseif ~isempty(excludeSample) && excludeType == "contains"
        exclMask = contains([outc.sampleID],excludeSample);
    end
    
    % find all samples that cannot be corrected
    exclMask = exclMask | isnan([outc.s_corr]);
    
    % % % % % % % % % % % % % % % % % 
    
    % indices of samples of interest:
    s_idx = find(~exclMask)';
    
    % apply mask
    outc = outc(s_idx);
    fracMaskcorrAnc = fracMaskcorrAnc(s_idx);
    
    if isempty(outc)
        error("Remove " +facsdata(c)+ " ! All datapoints are excluded from the analysis!")
    end
    
    % sort the samples numerically according to their number
    numbers = regexp([outc.sampleID],'\d*','Match');  % extract the numbers
    
    % if there are more numbers in some names, then dont sort
    if sum(cellfun(@(x) numel(x)>1,numbers))>0
        sortOrder = [1:numel(numbers)];  
    else
        [~,sortOrder] = sort(double([numbers{:}]));       % sort
    end
    
    % % % %
    % this is the preliminary data 
    outc = outc(sortOrder);
    fracMaskcorrAnc = fracMaskcorrAnc(sortOrder);
    
    samples{c} = unique([outc.sampleID],'stable');
    
    % with the data at this point, create a frac mask
    fracMask = ~([outc.start_s]<minFrac | [outc.start_s]>maxFrac);
    
    %% Plotting results of the single runs
    
    if plotSingleRuns == "ON"
        
        fprintf("Plotting results ... \n");
        
        % Plot with UNCORRECTED s ... figure(1),(7),(13)...
        figure(c*6-5);
        title("Sample" + c + "-- the UNCORRECTED s")
        hold on;
        set(gcf, 'Position', [25 570 560 420], 'Renderer', 'painters')
        plot([0.5 length(samples{c})+0.5],[0 0], '--', 'Color', [0.2 0.2 0.2])
        
        tmpS = [outc.s];
        tmpID = [outc.sampleID];
        boxplot(tmpS(fracMask),tmpID(fracMask)); ylabel('selection coeff uncorrected')
        
        if length(samples{c}) > 4
            set(gca, 'XTickLabelRotation', 15)
        end
        
        % Plot with CORRECTED s ... figure(2),(8),(14)...
        figure(c*6-4); hold on;
        title("Sample" + c + "-- the CORRECTED s")
        set(gcf, 'Position', [600 570 560 420], 'Renderer', 'painters')
        plot([0.5 length(samples{c})+0.5],[0 0], '--', 'Color', [0.2 0.2 0.2])
        tmpSCorr = [outc.s_corr];
        tmpID = [outc.sampleID];
        boxplot(tmpSCorr(fracMask),tmpID(fracMask));
        ylabel('selection coeff corrected')
        
        if length(samples{c}) > 4
            set(gca, 'XTickLabelRotation', 15)
        end
        
    end
    
    if plotSystematics == "ON"
        % Plot CORRECTED s against start frac for each plate. .. figure(310)..
        % Plot CORRECTED s against meas. order for each plate. .. figure(410)..
        %%% Is there a systematic? %%%
        
        outc_plateMask = zeros(4,size(outc,2));
        
        for pl=0:5
            outc_plateMask = [];
            outc_plateMask = [outc.plate]==pl;
            
            if sum(outc_plateMask) > 0
                corrMask = ismember([outc.sampleID],corrANCsamples);
                sanity   = [outc.start_s] > 0 & ~isnan([outc.s_corr]);
                
                collect_s_corr    = [outc(outc_plateMask & sanity).s_corr];
                collect_start_s   = [outc(outc_plateMask & sanity).start_s];
                collect_ende_s    = [outc(outc_plateMask & sanity).ende_s];
                collect_timeStart = [outc(outc_plateMask & sanity).startTime];
                collect_timeEnde = [outc(outc_plateMask & sanity).endeTime];
                lengthStartMeas = max(collect_timeStart) - min(collect_timeStart);
                lengthEndeMeas = max(collect_timeEnde) - min(collect_timeEnde);
                
            figure(300 + 10*c); hold on; set(gcf, 'Renderer', 'painters'); box on;
                m_fD     = scatter(collect_start_s,collect_s_corr,'*','MarkerEdgeColor',cmp(pl,:), 'LineWidth', 1.2); hold on; %cmp_runs(pl,:)); hold on
                m_fDcorr = scatter([outc(outc_plateMask & sanity & corrMask).start_s],[outc(outc_plateMask & sanity & corrMask).s_corr],'d','MarkerEdgeColor','m'); hold on
                
            figure(400 + 10*c + pl); set(gcf, 'Position', [1200 70 560 850], 'Renderer', 'painters')
                subplot(3, 1, 1); box on; hold on;
                m_fD2     = scatter(collect_timeStart,collect_s_corr,'*','MarkerEdgeColor','k'); hold on%cmp_runs(pl,:)); hold on
                m_fD2corr = scatter([outc(outc_plateMask & sanity & corrMask).startTime],[outc(outc_plateMask & sanity & corrMask).s_corr],'d','MarkerEdgeColor','m');
                local_corrMean = mean([outc(outc_plateMask & sanity & corrMask).s_corr]);
                plot([0 1e7],[local_corrMean local_corrMean],'m--')
                
                
            figure(300 + 10*c);
                xlim([minFrac-10 maxFrac+10]);
                ax = gca;
                ylim([min([collect_s_corr ax.YLim+0.1])-0.1 max([collect_s_corr ax.YLim-0.1])+0.1]);
                
                if numel(collect_start_s) >= 5
                    [corrCov,~] = corrcoef(collect_start_s,collect_s_corr);
                    text(max(collect_start_s)-5,max(collect_s_corr)+0.05,"Correlation = " + num2str(corrCov(2,1), "%0.3f"), "Color", cmp(pl,:));
                end
                
                title("Dependency start fraction and s" + "-- sample" + string(c))
                xlabel('Start Frac Sample'); ylabel('Selection Coeff corrected')
                
            figure(400 + 10*c + pl); subplot(3, 1, 1);
                xlim([min(collect_timeStart)-50 max(collect_timeStart)+50]);
                ylim([min(collect_s_corr)-0.1 max(collect_s_corr)+0.1]);
                
                if numel(collect_timeStart) >= 5
                    [corrCov,~] = corrcoef(collect_timeStart,collect_s_corr);
                    text(max(collect_timeStart)-0.35*lengthStartMeas,max(collect_s_corr+0.05),"Correlation = " + num2str(corrCov(2,1), "%0.3f"))
                end
                
                title("Dependency on order of measurement -- sample " + string(c) + " -- plate" + string(pl))
                xlabel('Time [s]'); ylabel('Selection Coeff Corrected')
                
            figure(400 + 10*c + pl); subplot(3, 1, 2); box on;
                hold on; ylabel("Start Fraction [%]")
                plot(collect_timeStart, collect_start_s, "kx", "LineWidth", 1.4)
                xlim([min(collect_timeStart)-50 max(collect_timeStart)+50]);
                ylim([20 80])
                xlabel("Duration: " + num2str(lengthStartMeas/60, "%0.1f") + " min")
                if numel(collect_timeStart) >= 5
                    [corrCov,~] = corrcoef(collect_timeStart, collect_start_s);
                    text(max(collect_timeStart)-0.35*lengthStartMeas,max(collect_start_s)+10,"Correlation = " + num2str(corrCov(2,1), "%0.3f"))
                end
                
            figure(400 + 10*c + pl); subplot(3, 1, 3); box on;
                 hold on; ylabel("End Fraction [%]")
                plot(collect_timeEnde, collect_ende_s, "kx", "LineWidth", 1.4)
                xlim([min(collect_timeEnde)-50 max(collect_timeEnde)+50]);
                ylim([20 80]);
                xlabel("Duration: " + num2str(lengthEndeMeas/60, "%0.1f") + " min")
                if numel(collect_timeEnde) >= 5
                    [corrCov,~] = corrcoef(collect_timeEnde, collect_ende_s);
                    text(max(collect_timeEnde)-0.35*lengthEndeMeas,max(collect_ende_s)+10,"Correlation = " + num2str(corrCov(2,1), "%0.3f"))
                end
            end
        end
        
    end
    
    %% Exclude samples
    
    % % -- all samples with too high/low start fractions 
    
    outc = outc(fracMask & ~contains([outc.sample], corrANCsamples) | fracMaskcorrAnc & contains([outc.sample], corrANCsamples) );
    samples{c} = unique([outc.sampleID],'stable');
    
    samplesLost{c} = names_all(~ismember(names_all,samples{c}));
    
    
    
    %% Calculating means and stds, 
    % % -- saving them
    % % -- plotting histograms   
    % (here, the values with too high/low start fraction are excluded!!
    
    meanS = []; errorS = []; % -> the arrays will fit samples{c}
    for i = 1 : length(samples{c})
        
        idx = find(strcmp(samples{c}(i), [outc.sampleID]));
        
        if BCcond == "ON"
            tmpSCorr = [outc(idx).s_corr- BCcorr];
        else
             tmpSCorr = [outc(idx).s_corr];
        end     
        tmpSError= [outc(idx).s_corrErr];
        
        % Here, instead of just calculating the std, we take the errors into
        % account and evaluate error of mean with error propagation
        [meanS(i),errorS(i)] = getMeanError(tmpSCorr,tmpSError);
        
    end
    
    samplesC = cellfun(@(x) {string(x)},cellstr(samples{c})); % convert samples to cell
    
    if ~exist("outSIndividual","var")
        % Saving all information about the individual samples
        outSIndividual = struct("Sample", {outc.sampleID},"Exp", num2cell(c*ones(1,length(outc))), "startFrac", num2cell([outc.start_s]), "endFrac", num2cell([outc.ende_s]),...
            "s", num2cell([outc.s_corr]),"sErr", num2cell([outc.s_corrErr]));
        % Saving the mean and std per sample PER RUN
        outS = struct("Sample", samplesC, "Run", num2cell(c*ones(1,length(samples{c}))), "meanS", num2cell(meanS), "errorS", num2cell(errorS));
    else
        outSaveAdd = struct("Sample", {outc.sampleID},"Exp", num2cell(c*ones(1,length(outc))), "startFrac", num2cell([outc.start_s]), "endFrac", num2cell([outc.ende_s]),...
            "s", num2cell([outc.s_corr]),"sErr", num2cell([outc.s_corrErr]));
        outSIndividual = [outSIndividual outSaveAdd];
        outSAdd = struct("Sample", samplesC, "Run", num2cell(c*ones(1,length(samples{c}))), "meanS", num2cell(meanS), "errorS", num2cell(errorS));
        outS = [outS outSAdd];
        clear *Add
    end
    
    
    % overlay the histograms
    f10 = figure(10); set(gcf, 'Position', [550 350 960/2 420/2], 'Renderer', 'painters');
    hold on; box on;
    set(gca, 'FontSize', 9);
    xlabel("selection coefficient"); ylabel("Counts");
    
    h_fDSingle = histogram(meanS, "BinEdges", histBinEdges, "Normalization", "probability");
    h_fDSingle.FaceColor = "none";
    h_fDSingle.EdgeColor = cmp(c,:);
    h_fDSingle.LineWidth = 1.2;
    
    m_fD = plot([mean(meanS) mean(meanS)], [0 4], "--", "LineWidth", 1.2,'Color',cmp(c,:));
    
    histmax = [histmax max(h_fDSingle.Values)];
    hists = [hists,h_fDSingle,m_fD];
    dateTmp = regexp(facsdata(c), "\d*", "Match"); dateMeasured = dateTmp(cellfun("length", dateTmp)==8);
    histLegends = [histLegends,"Run "+c + ": " + dateMeasured,"Mean Run " + c];
    legend(hists,histLegends,'Location','NorthWest')
    
    ylim([0 max(histmax) + 0.05]);
    
end

%% Plotting the results merged together

% ----------------------------------------------------------------------
% ------------------PLOT ------- outSIndividual ------------------------
% -----------------------------------------------------------------------

% if plotComparingSingleRuns == "ON"
    
%     % Give every sampleSet a unique name to not mix up the same replicates
%     % measured at different days
%     for i = 1 : length(outSIndividual)
%         outSIndvNames(i) = outSIndividual(i).Sample + "_Run" + num2str(outSIndividual(i).Exp);
%     end
%     
%     % Sort the samples for a better over view
%     [sortSIndv, sortIdxSIndv] = sort(outSIndvNames);
%     unsortSIndv = [outSIndividual(:).s];
%     unsortSIndvErr = [outSIndividual(:).sErr];
%     sortedSIndv = unsortSIndv(sortIdxSIndv);
%     sortedSIndvErr = unsortSIndvErr(sortIdxSIndv);
    
%     % Plot the selection coefficent seperately for each run
%     f100 = figure(100); hold on; grid on;
%     set(gcf, 'Position', [0 50 760/1.9 420/2], 'Renderer', 'painters') % for only one set of data
%     set(gca, 'FontSize', 9);
%     ylabel("Selection Coefficient");
%     % set(gcf, 'Position', [0 50 960 420], 'Renderer', 'painters')
%     set(gca, 'XTickLabelRotation', 15)
%     boxplot(sortedSIndv,cellstr(sortSIndv))
%     plot([0.5 length(unique(sortSIndv))+0.5],[0 0], '--', 'Color', [0.2 0.2 0.2])
%     ylim([min(sortedSIndv)-0.1 max(sortedSIndv)+0.1])
%     
    % plot each run with error
%     f101 = figure(101); hold on; grid on;
%     set(gcf, 'Position', [0 50 760/1.9 420/2], 'Renderer', 'painters') % for only one set of data
%     set(gca, 'XTick',[1:numel(outSIndvNames)],'XTickLabel',outSIndvNames,'XTickLabelRotation', 15)
%     ylabel("Selection Coefficient");
%     % set(gcf, 'Position', [0 50 960 420], 'Renderer', 'painters')
%     set(gca, 'XTickLabelRotation', 15)
%     errorbar(sortedSIndv,sortedSIndvErr,'o');
%     ylim([min(sortedSIndv)-0.1 max(sortedSIndv)+0.1])
% end


%% Collecting all the runs for each sample - create outCollect
clear sortOrder numbers numbers_temp

% Sort the sample summary outS for a better over view
numbers = regexp([outS.Sample],'\d*','Match');

% if there are more numbers in some names, then take first
if sum(cellfun(@(x) numel(x)>1,numbers))>0
    sortOrder = [1:numel(numbers)];
else
    [~,sortOrder] = sort(double([numbers{:}]));       % sort
end

MeanSamples_sort = [outS(sortOrder).Sample];
MeanS_sort = [outS(sortOrder).meanS];
Run_sort = [outS(sortOrder).Run];

% Calculating means and std again!: (over all experiment days!)
samplesOutput = unique(MeanSamples_sort,'stable');

for i = 1 : length(samplesOutput)
    idx = find(strcmp(samplesOutput{i}, [outS(:).Sample]));
    dataPoints(i) = length(idx);
    
    tmpS      = [outS(idx).meanS];
    tmpSError = [outS(idx).errorS];
    
    [meanMeanS(i),ErrorMeanS(i)] = getMeanError(tmpS,tmpSError);
    stdErrorMeanS(i)   = std(tmpS) / sqrt(length(tmpS));
end

samplesOut = cellfun(@(x) {string(x)},cellstr(samplesOutput));


% Saving the mean and std per sample per run
outCollect = struct("Sample", samplesOut,"dataPoints", num2cell(dataPoints),...
    "meanS", num2cell(meanMeanS), "stdErrorS", num2cell(stdErrorMeanS), "errorS", num2cell(ErrorMeanS));

% % Exclude all samples from outCollect
% % % --with less than dataPointExcl data points
dataPoExclMask = [outCollect.dataPoints]<dataPointExcl;

% This is just printing while running the script, to know which samples are excluded
if sum(dataPoExclMask)>0
    fprintf("\nRemoving samples with less than "+ num2str(dataPointExcl)+ " data points........\n") ;
    exclSampl = [outCollect(dataPoExclMask).Sample];
    for i = 1 : sum(dataPoExclMask)
        fprintf("--> " + exclSampl(i) + "\n");
    end
    fprintf("\n");
end
outCollect_all = outCollect;
outCollect     = outCollect(~dataPoExclMask);

%% Plot collected data 
% ----------------------------------------------------------------------
% ------------------PLOT ------- outCollect----- ------------------------
% -----------------------------------------------------------------------

ymax = max([outSIndividual.s]) + 0.1;
ymin = min([outSIndividual.s]) - 0.05;

if plotErrorPlots == "ON"
    % for each sample the boxplot
    f200 = figure(200); hold on;
    ylabel("Selection Coefficient");
    title('s for each sample (each run contributes with max. 1 data point)')

    set(gcf, 'Position', [850 350 1000/1.5 1000/3], 'Renderer', 'painters')
    set(gca, 'XTick',1:length(MeanS_sort),'XTickLabelRotation', 15,'FontSize', 9)

    plot([0 length(MeanS_sort)+1],[0 0], '--', 'Color', [0.2 0.2 0.2])
    plot([1:length([outCollect_all.meanS])],[outCollect_all.meanS],'m*','HandleVisibility','off')
    boxplot(MeanS_sort,cellstr(MeanSamples_sort))
    ylim([ymin ymax])
end
    % for each sample the boxplot with supermuch infos extra
    f201 = figure(201); hold on; 
    ylabel("Selection Coefficient");
    title('s for each sample (each run contributes with max. 1 data point)')

    set(gcf, 'Position', [850 350 1000/1.5 1000/3], 'Renderer', 'painters')
    set(gca, 'XTick',1:length(MeanS_sort),'XTickLabelRotation', 50,'FontSize', 9)

    plot([0 length(MeanS_sort)+1],[0 0], '--', 'Color', [0.2 0.2 0.2],'HandleVisibility','off')
    boxplot(MeanS_sort,cellstr(MeanSamples_sort),'Colors',[0.8 0.8 0.8])
    plot([1:length([outCollect_all.meanS])],[outCollect_all.meanS],'kd','HandleVisibility','off','MarkerFaceColor','k', 'MarkerSize', 5)

    shift = linspace(-0.01,0.01,numel(unique(Run_sort)));
    % run over each run so that runs are plotted individually
    for i = unique(Run_sort)
        whereRun  = (Run_sort==i);
        
        [~,whereNames] = ismember(MeanSamples_sort(whereRun),[outCollect_all.Sample]);
        
        scatter(whereNames + shift(i),MeanS_sort(whereRun),'x','MarkerEdgeColor',cmp_runs(i,:),...
            'Displayname',"Run " + i, 'LineWidth', 1.2);
    end
    legend('NumColumns',4)
    %ylim([ymin ymax])
    ylim([-0.35 0.35])
    
if plotErrorPlots == "ON"
    
    % for teach sample the boxplot with supermuch infos extra
    f202 = figure(202); hold on;
    ylabel("Selection Coefficient");
    title('s for each sample -- each measurement is shown from diff. days')

    set(gcf, 'Position', [850 350 1000/1.5 1000/3], 'Renderer', 'painters')
    set(gca, 'XTick',1:length(MeanS_sort),'XTickLabelRotation', 15,'FontSize', 9)

    plot([0 length(MeanS_sort)+1],[0 0], '--', 'Color', [0.2 0.2 0.2],'HandleVisibility','off')
    boxplot(MeanS_sort,cellstr(MeanSamples_sort),'Colors',[0.8 0.8 0.8])
    plot([1:length([outCollect_all.meanS])],[outCollect_all.meanS],'dm','HandleVisibility','off','MarkerFaceColor','m')

    shift = linspace(-0.3,0.3,numel(unique([outSIndividual.Exp])));
    for i = unique([outSIndividual.Exp])
        whereS  = ([outSIndividual.Exp]==i);
        [~,whereNames] = ismember([outSIndividual(whereS).Sample],[outCollect_all.Sample]);

        scatter(whereNames,[outSIndividual(whereS).s],'x','MarkerEdgeColor',cmp_runs(i,:),...
            'Displayname',"Run " + i);
    end
    legend('NumColumns',4)
    ylim([ymin ymax])

    % for each sample the error bars from measurement inaccuracy
    f210 = figure(210); hold on;
    ylabel("Selection Coefficient");

    set(gcf, 'Position', [850 350 1000/1.5 1000/3], 'Renderer', 'painters')
    set(gca, 'XTick',1:length([outCollect.Sample]),'XTickLabel',[outCollect.Sample],'XTickLabelRotation', 15,'FontSize', 9)
    errorbar([outCollect.meanS],[outCollect.errorS],'k.','HandleVisibility','off');
    plot(1:length([outCollect.Sample]),[outCollect.meanS],'m*','Displayname','s mean with error');

    plot([0 length([outCollect.Sample])+1],[0 0],'--','Color',[0.2 0.2 0.2],'HandleVisibility','off')
    legend('NumColumns',4)
    ylim([ymin ymax])

    % for each sample the standard mean error
    f220 = figure(220); hold on;
    ylabel("Selection Coefficient");

    set(gcf, 'Position', [850 350 1000/1.5 1000/3], 'Renderer', 'painters')
    set(gca, 'XTick',1:length([outCollect.Sample]),'XTickLabel',[outCollect.Sample],'XTickLabelRotation', 15,'FontSize', 9)
    errorbar([outCollect.meanS],[outCollect.stdErrorS],'k.','HandleVisibility','off');
    plot(1:length([outCollect.Sample]),[outCollect.meanS],'m*','Displayname','s mean with error');

    plot([0 length([outCollect.Sample])+1],[0 0],'--','Color',[0.2 0.2 0.2],'HandleVisibility','off')
    legend('NumColumns',4)
    ylim([ymin ymax])
end




%% Fitness distribution of all runs taken together
% but if there were ancestors, then these are plotted ontop

f11 = figure(11); set(gcf, 'Position', [1150 350 990/2 440/2], 'Renderer', 'painters');
title("All samples against the ancestor samples used for correction")
hold on; box on;
set(gca, 'FontSize', 9);

h_fD = histogram([outCollect(~ismember([outCollect.Sample], corrANCsamples)).meanS], 'BinEdges', histBinEdges,...
    'Normalization', 'count','FaceColor',cmp(1,:),'EdgeColor','none','LineWidth',1.5);
h_fDAnc = histogram([outCollect(ismember([outCollect.Sample], corrANCsamples)).meanS], 'BinEdges', histBinEdges,...
    'Normalization', 'count','FaceColor',[0.2 0.2 0.2],'FaceAlpha', 1,'EdgeColor','none','LineWidth',1.5);

ylabel("Counts");

ylim([0 1.15*max([h_fD.Values h_fDAnc.Values])]);
xlim([-max(abs(h_fD.BinLimits))-h_fD.BinWidth max(abs(h_fD.BinLimits))+h_fD.BinWidth]);

mean_fD    = mean([outCollect(~ismember([outCollect.Sample], corrANCsamples)).meanS]);
std_fD     = std([outCollect(~ismember([outCollect.Sample], corrANCsamples)).meanS]);
mean_fDAnc = mean([outCollect(ismember([outCollect.Sample], corrANCsamples)).meanS]);
std_fDAnc  = std([outCollect(ismember([outCollect.Sample], corrANCsamples)).meanS]);

m_fD = plot([mean_fD mean_fD], [0 1.15*max([h_fD.Values h_fDAnc.Values])], "-",...
    "LineWidth", 0.8, "Color", cmp(1,:));


%%% Here check if there are any ancestor data points -- if not, then skip %
%%% plotting the data in the legend

if all(isnan(h_fDAnc.Values))
    legend([h_fD m_fD ], "DFE library",...
        "Mean = " + num2str(mean_fD, "%0.4f") + "\newlineStd = " + num2str(std_fD, "%0.4f"),...
        'Location','NorthWest')
else
    m_fDAnc = plot([mean_fDAnc mean_fDAnc], [0 1.15*max([h_fD.Values h_fDAnc.Values])], "-",...
        "LineWidth", 0.8, "Color", [0.2 0.2 0.2 0.8]);
    legend([h_fD h_fDAnc m_fD m_fDAnc], "DFE library", "corr "+ ancestor + " samples",...
        "Mean = " + num2str(mean_fD, "%.1d") + "\newlineStd = " + num2str(std_fD, "%.1d"),...
        "Mean = " + num2str(mean_fDAnc, "%.1d") + "\newlineStd = " + num2str(std_fDAnc, "%.1d"),...
        'Location','NorthWest')
end



%% Now compare the ancestor distribution to the evaluated data.
% Plot the distributions with mean and std
% do the ks test --  are they the same?

if compareAnc == "ON"
    
    % take the locally analysed data from outCollect
    % BUT --> avoid that ancestors get mixed up in the Lib, write output
    noAncs = ~ismember([outCollect.Sample], corrANCsamples);
    LibmeanS = [outCollect(noAncs).meanS];
    LibMeanSBC = LibMeanS - BCcorr; 
    
    fprintf("--> %i ancestor sample*s in your local data excluded \n",sum(noAncs==0))
    
    %%%
    collectAnc = load(AncDist);
    AncmeanS = [collectAnc.outCollect(:).meanS];
    AncMeanSBC = AncMeanS - BCcorr;
    
    f30 = figure(30); set(gcf, 'Position', [550 350 700 420], 'Renderer', 'painters');

    hold on; box on;
    set(gca, 'FontSize', 15);
    
    h_fD = histogram(LibmeanS, 'BinEdges', histBinEdges,...
        'Normalization', 'probability','FaceColor','none','EdgeColor',cmp(1,:),'LineWidth',1.5);
    h_fDAnc = histogram(AncmeanS, 'BinEdges', histBinEdges,...
        'Normalization', 'probability','FaceColor',[0.2 0.2 0.2],'FaceAlpha', 0.4, 'EdgeColor','none','LineWidth',1.5);
    
    ylabel("Probability");
    
    ylim([0 1.15*max([h_fD.Values h_fDAnc.Values])]);
    %xlim([-max(abs(h_fD.BinLimits))-h_fD.BinWidth max(abs(h_fD.BinLimits))+h_fD.BinWidth]);
    xlim([-0.4 0.4])
    mean_fD    = mean(LibmeanS);
    std_fD     = std(LibmeanS);
    mean_fDAnc = mean(AncmeanS);
    std_fDAnc  = std(AncmeanS);
    
    m_fD = plot([mean_fD mean_fD], [0 1.15*max([h_fD.Values h_fDAnc.Values])], "--",...
        "LineWidth", 1.2, "Color", [cmp(1,:) 0.5]);
    m_fDAnc = plot([mean_fDAnc mean_fDAnc], [0 1.15*max([h_fD.Values h_fDAnc.Values])], "-.",...
        "LineWidth", 1.2, "Color", [0.2 0.2 0.2 0.5]);
    legend([h_fD m_fD m_fDAnc], ExpOutName,...
        "Mean = " + num2str(mean_fD, "%.4f") + "\newlineStd = " + num2str(std_fD, "%.4f"),...
        "Mean = " + num2str(mean_fDAnc, "%.4f") + "\newlineStd = " + num2str(std_fDAnc, "%.4f"),...
        'Location','NorthWest')
    
    [h,p] = kstest2(AncmeanS,LibmeanS);
    if h==1
        testResult = "The null hypothesis was rejected, the distributions are significantly different! \n";
    elseif h==0
        testResult = "The null hypothesis was NOT rejected, the distributions are the same!";
    end
    sprintf(testResult + "The p-value is %0d",p)
    
else
    noAncs = ~ismember([outCollect.Sample], corrANCsamples);
    LibMeanS = [outCollect(noAncs).meanS];
 %   LibMeanSBC = LibmeanS - BCcorr; 

     mean_fD    = mean(LibMeanS);
    std_fD     = std(LibMeanS);
    
    fprintf("--> %i ancestor sample*s in your local data excluded \n",sum(noAncs==0))
    
    f30 = figure(30); set(gcf, 'Position', [600 350 600 260], 'Renderer', 'painters');
    hold on; box on;
    set(gca, 'FontSize', 12);
    
    h_fD = histogram(LibMeanS, 'BinEdges', histBinEdges,...
        'Normalization', 'probability','FaceColor','none','EdgeColor',cmp(1,:),'LineWidth',1.5);
    m_fD = plot([mean_fD mean_fD], [0 1.15*max([h_fD.Values h_fDAnc.Values])], "--",...
        "LineWidth", 1.2, "Color", [cmp(2,:)]);
    
    ylabel("Probability"); title("Distribution of fitness effects");
    
    ylim([0 round(1.15*max(h_fD.Values),2)]);
    xlim([-0.1 0.1])
    xticks([-0.1 -0.05 0 0.05 0.1])
   
    legend([h_fD m_fD], ExpOutName,...
        "Mean = " + num2str(mean_fD, "%.4f") + "\newlineStd = " + num2str(std_fD, "%.4f"),...
        'Location','NorthEast')
end

%outCollect = outCollect(noAncs);



%% print
outCollect = outCollect(noAncs)
save(basePath + "out/BCcorr/me/" + datestr(now, "yyyymmdd") + '_' + ExpOutName+ '_' + "outCollect", "outCollect")

if printFigs == "ON"
    print(f10,  '-painters', '-dpng', savePlotPath + datestr(now, "yyyymmdd") + "_" + ExpOutName + "_fitnessDist_separat");
    print(f11,  '-painters', '-dpng', savePlotPath + datestr(now, "yyyymmdd") + "_"+ ExpOutName + "_plate2_20210728_fitnessDist_wCorrSamples");
    
    if plotErrorPlots == "ON"
        print(f200,  '-painters', '-dpng', savePlotPath + datestr(now, "yyyymmdd") + "_"  + ExpOutName + "_BoxPlots" );
        print(f210,  '-painters', '-dpng', savePlotPath + datestr(now, "yyyymmdd") + "_"  + ExpOutName + "_MeansSError" );
        print(f220,  '-painters', '-dpng', savePlotPath + datestr(now, "yyyymmdd") + "_"  + ExpOutName + "_MeansSstdError" );
        print(f201,  '-painters', '-dpng', savePlotPath + datestr(now, "yyyymmdd") + "_"  + ExpOutName + "_BoxPlots+" );
    end 
    
    if compareAnc == "ON"
        print(figure(30),  '-painters', '-dpng', savePlotPath + datestr(now, "yyyymmdd") + "_" + ExpOutName  + "_fitnessDist" );
    end
end

%print(figure(30),  '-painters', '-dpng', savePlotPath + datestr(now, "yyyymmdd") + "_" + ExpOutName  + "_fitnessDist" );
