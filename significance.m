clearvars
close all

Path = "/home/ariana/Documents/CompExpRobot/EvolandComp_Bs210_Bs224/out/BCcorr/me/" ;
%Path =  "/home/ariana/Documents/MATLAB/Mona/";

ctlData= load(Path+"20241118_Bs224ctl_outCollect.mat") ;
%ctlData= load(Path+"20240325_Bs210CtlEop5_outCollect_6runs_filCorr_v4.mat") ;

hybData= load(Path+"20241118_Bs224hyb_outCollect.mat") ;
%hybData= load(Path+"20240325_Bs210HybEop5_outCollect_4runs_filCorr_v4.mat") ;

ancData = load(Path+"20241125_Bs210ANC_outCollect.mat") ;

%comp = [ctlData.LibMeanSBC; hybData.LibMeanSBC]




test = kstest2([ctlData.outCollect.meanS], [hybData.outCollect.meanS]);
        
ctlCon = bootci(10000, @mean,[ctlData.outCollect.meanS]);
hybCon = bootci(10000, @mean, [hybData.outCollect.meanS]);




%% %% Outlier bestimmen %%%%%%%%%%%%%

out_clt_mask=zeros(2, length([ctlData.outCollect(:).meanS]));

for n=1:length([ctlData.outCollect(:).meanS])


[out_clt_mask(1,n), out_clt_mask(2, n)] = ztest([ctlData.outCollect(n).meanS], mean([ancData.outCollect(:).meanS]), std([ancData.outCollect(:).meanS]), "Alpha", 0.05/length([ctlData.outCollect(:).meanS]));


end

out_hyb_mask=zeros(2, length([ctlData.outCollect(:).meanS]));

for n=1:length([hybData.outCollect(:).meanS])


[out_hyb_mask(1,n), out_hyb_mask(2, n)] = ztest([hybData.outCollect(n).meanS], mean([ancData.outCollect(:).meanS]), std([ancData.outCollect(:).meanS]), "Alpha", 0.05/length([hybData.outCollect(:).meanS]));

end


out_ctl=ctlData.outCollect(logical(out_clt_mask(1, :)));
out_ctl_p= out_clt_mask(2, logical(out_clt_mask(1, :)));
rand_out_ctl = randsample(out_ctl, 2);

out_hyb=hybData.outCollect(logical(out_hyb_mask(1, :)));
out_hyb_p= out_hyb_mask(2, logical(out_hyb_mask(1, :)));
rand_out_hyb = randsample(out_hyb, 2);


%% choose random sample 

letter = ["A" "B" "C" "D" "E" "F" "G" "H"];
numbers = 1:12;

rand_letter=randsample(letter, 5, true);
rand_number=randsample(numbers , 5, true);
