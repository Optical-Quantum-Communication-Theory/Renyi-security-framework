baseName = "RenyiPhaseImpDecoy46_1e8Results";
fileNameFormat = baseName +"_%d.mat";
numJobs = 12;

tmpCell = cell(numJobs,1);
results = struct("debugInfo",tmpCell,"keyRate",tmpCell,"currentParams",tmpCell);

for index = 1:numJobs
    tmpFile = load(sprintf(fileNameFormat,index));

    if index == 1
        qkdInput = tmpFile.qkdInput;
    end
    
    results(index) = tmpFile.results;
end
save(baseName+".mat","qkdInput","results");

%%

QKDPlot.plotParameters({results},"transmittance",QKDPlot.keyRateTag,"xScaleStyle","dB","yScaleStyle","log")