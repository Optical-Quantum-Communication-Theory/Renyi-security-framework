%pick preset
qkdInput = RenyiDecoy46LossyPreset();

%run the QKDSolver with this input and store results
results = MainIteration(qkdInput);

%save the results and preset to a file.
save("RenyiDecoy46LossyResultsNewOpt.mat","results","qkdInput");

%% plot the result
keyRateFixed = arrayfun(@(x) x.debugInfo.keyRateModule.keyRateFixed, results);
fig = axes(figure());
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log","figAxis",fig)

hold(fig,"on")
trans = qkdInput.scanParameters.transmittance;

db = -10*log10(cell2mat(trans));

plot(fig,db,keyRateFixed,"+-","DisplayName","Fixed Rate")

%plot asymptotic rate
tdb = linspace(-10*log10(trans{1}),-10*log10(trans{end}),100);
t = 10.^(-tdb/10);

musig = qkdInput.fixedParameters.GROUP_decoys_1;
ptest = qkdInput.fixedParameters.probTest;
probsB = qkdInput.fixedParameters.probsB;
pzB = probsB(1);

plot(tdb,musig*exp(-musig)*(1-ptest)*pzB*t,"DisplayName","Asymptotic Rate")

lgd = legend();
lgd.FontSize = 10;
lgd.Location = 'best';

hold(fig,"off")