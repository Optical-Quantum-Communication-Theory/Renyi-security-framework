function RenyiPhaseImpDecoy46batch(basePath,baseFileName,jobIndex,totalNumJobs)
% This is just a hacked together script to compute the points in parallel
% because they take a lot of time.
arguments
    basePath (1,1) string {mustBeFolder}
    baseFileName (1,1) string
    jobIndex (1,1) uint64 {mustBePositive}
    totalNumJobs (1,1) uint64 {mustBeGreaterThanOrEqual(totalNumJobs,jobIndex)}
end

qkdInput = RenyiPhaseImpDecoy46LossyPreset();

%total number of signals sent
Ntot = 1e10;
qkdInput.addFixedParameter("Ntot", Ntot); 

%Loss
%total array of loss values to iterate over
lossdB = linspace(0,40,21);
lossdB = lossdB(1:16);
transmittance = 10.^(-lossdB/10);

[startIndex,endIndex] = selectPoints(jobIndex,totalNumJobs,numel(transmittance));

scanIndexList = startIndex:endIndex;
numScanList= numel(scanIndexList);

%filestring for optimal values (specifically for us the logrenyiAlpha)
filestrOptVals = "optValsDecoy46_N=";
fileStrTemp = filestrOptVals + sprintf("%.2e",Ntot) +"_q=9.90e-01.csv";
optvals = readmatrix(fileStrTemp);


%initialize results
results = struct("debugInfo",cell(numScanList,1),"keyRate",cell(numScanList,1),"currentParams",cell(numScanList,1));

for indexList = 1:numScanList
    indexLoss = scanIndexList(indexList);

    qkdInput.addFixedParameter("transmittance", transmittance(indexLoss));

    %Renyi param
    %optimized alpha
    logAlpha = optvals(indexLoss,1);
    bndsLogAlpha = lowerUpperBnds_from_optvals(indexLoss,optvals(:,1),-5,-0.5);
    logrenyiAlpha.lowerBound = bndsLogAlpha(1);
    logrenyiAlpha.upperBound = bndsLogAlpha(2);
    logrenyiAlpha.initVal = logAlpha;
    qkdInput.addOptimizeParameter("logrenyiAlpha", logrenyiAlpha);
    
    %run the QKDSolver with this input and store results
    disp("starting Main Iteration. Good Luck!")
    results(indexList) = MainIteration(qkdInput);

end

fileName = sprintf("%s_%d.mat",baseFileName,jobIndex);
fileName = fullfile(basePath,fileName);


qkdInput.addScanParameter("transmittance", num2cell(transmittance));

save(fileName,"results","qkdInput","scanIndexList");

end


function [startIndex,endIndex] = selectPoints(jobIndex,totalNumJobs,numIterations)
arguments
    jobIndex (1,1) uint64
    totalNumJobs (1,1) uint64
    numIterations (1,1) uint64
end

basePointsPerJob = idivide(numIterations,totalNumJobs,"floor");
remainder = numIterations - basePointsPerJob*totalNumJobs;

if jobIndex <= remainder
    % the first jobs get the extra remainder
    startIndex = (basePointsPerJob + 1) * (jobIndex-1) + 1;
    endIndex   = startIndex + basePointsPerJob;
else
    startIndex = basePointsPerJob * (jobIndex -1) + remainder + 1;
    endIndex = startIndex + basePointsPerJob - 1;
end
end