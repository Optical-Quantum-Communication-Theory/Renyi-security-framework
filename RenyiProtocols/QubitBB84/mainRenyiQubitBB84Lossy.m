%pick preset
clear all
qkdInput = RenyiBB84LossyPreset();

%List of mutiple total signals sent
N_list = [1e4,1e5,1e6,1e7,1e8,1e9,1e10];

%filestring for optimal values
filestrOptVals = "optimalValues\optValsQubitBB84_N=";

%Loss stored in array with loss values to iterate over
lossdB = linspace(0,50,26);
transmittance = 10.^(-lossdB/10);

%list of maximal element of loss
lossList = [2,8,12,16,20,24,26];

% store qkdInput for later
qkdInputSave = qkdInput;

for indexSignals = 1:numel(N_list)
    %set transmittance list for current value of N
    transmittanceTemp = transmittance(1:lossList(indexSignals));

    %Load optimal values for current number of signal values
    fileStrTemp = filestrOptVals + sprintf("%.2e",N_list(indexSignals)) +".csv";
    %optimal values are sorted in coulmns as | logRenyiAlpha | probTest |
    optvals = readmatrix(fileStrTemp);

    for indexLoss = 1:numel(transmittanceTemp)
        fprintf("Iteration %.0f of %.0f for %.0e",indexLoss, numel(transmittanceTemp),N_list(indexSignals))

        %Add total signals sent from list above
        qkdInput.addFixedParameter("Ntot",N_list(indexSignals))
    
        %Add loss until element from list above
        qkdInput.addScanParameter("transmittance", num2cell(transmittance(indexLoss)));  
    
        %Add Renyi param from optimal values
        % fixed alpha
        logAlpha = optvals(indexLoss,1);
        qkdInput.addFixedParameter("logrenyiAlpha", logAlpha);

        %Add probTest from optimal values
        %fixed probTest
        probTest = optvals(indexLoss,2);
        qkdInput.addFixedParameter("probTest", probTest);

        % run the QKDSolver with this input
        results(indexLoss) = MainIteration(qkdInput);
    end
    %edit qkdinput to save correct results
    qkdInput.addScanParameter("transmittance", num2cell(transmittanceTemp));

    %filestring for saving
    filestr = sprintf("data/RenyiBB84LossyResults_%.2e",N_list(indexSignals)) + "_depol.mat";

    % save the results and preset to a file
    %save results
    results = results(1:numel(transmittanceTemp));
    save(filestr,"results","qkdInput");
end

%Plot result for N=1e10 (or last value in N_list)
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")

