%pick preset
qkdInput = RenyiDecoyBB84ActiveLossyPreset();

%run the QKDSolver with this input and store results
results = MainIteration(qkdInput);

%save the results and preset to a file.
save("RenyiDecoyBB84LossyResults.mat","results","qkdInput");

%% plot the result
keyRateFixed = arrayfun(@(x) x.debugInfo.keyRateModule.keyRateFixed, results);
fig = axes(figure());
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log","figAxis",fig)

hold(fig,"on")
trans = qkdInput.scanParameters.transmittance;

db = -10*log10(cell2mat(trans));

plot(fig,db,keyRateFixed,"+-")

hold(fig,"off")
