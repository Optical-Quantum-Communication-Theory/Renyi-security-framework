%pick preset
qkdInput = RenyiBB84LossyPreset();

%run the QKDSolver with this input and store results
results = MainIteration(qkdInput);

%save the results and preset to a file.
% save("RenyiBB84LossyResults.mat","results","qkdInput");

%% plot the result
keyRateFixed = arrayfun(@(x) x.debugInfo.keyRateModule.keyRateFixed, results);
fig = axes(figure());
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log","figAxis",fig)

hold(fig,"on")
trans = qkdInput.scanParameters.transmittance;

db = -10*log10(cell2mat(trans));

plot(fig,db,keyRateFixed,"+-")

hold(fig,"off")


%% plot Frank Wolfe Progress
% fVals = cell2mat(results.debugInfo.keyRateModule.mathSolver.fValQueue.toCell);
% gaps = cell2mat(results.debugInfo.keyRateModule.mathSolver.gapsQueue.toCell);
% points = results.debugInfo.keyRateModule.mathSolver.pointsQueue.toCell;
% grads = results.debugInfo.keyRateModule.mathSolver.gradQueue.toCell;
% 
% renyiAlpha = results.currentParams.renyiAlpha;
% renyiAlphaHat = 1/(2-renyiAlpha);
% 
% penaltyTerms = cellfun(@(x)x(end),points(1:end-1))/(renyiAlphaHat-1);
% renyiTerms = fVals-penaltyTerms;
% 
% 
% iter = 1:numel(gaps);
% semilogy(iter,fVals,"x-","DisplayName","fValUpper");
% hold on
% semilogy(iter,fVals-gaps,"+-","DisplayName","fValLower");
% semilogy(iter,gaps,"o-","DisplayName","gap");
% semilogy(iter,fVals-max(fVals-gaps),"*-","DisplayName","fValUpper-max(fValLower)");
% semilogy(iter,penaltyTerms,"s-","DisplayName","penaltyTerm");
% semilogy(iter,fVals-penaltyTerms,"--","DisplayName","RenyiTerm")
% hold off; legend();

%%


% stepPoints = results.debugInfo.keyRateModule.mathSolver.stepQueue.toCell;


% index = 19;
% stepPointIn = cellfun(@(x,y)FrankWolfe.inProdR(stepPoints{index}-x,y)<=0,points(1:index),grads(1:index)).'

