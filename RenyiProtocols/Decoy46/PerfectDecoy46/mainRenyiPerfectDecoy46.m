%pick preset
clear all
qkdInput = RenyiDecoy46LossyPreset();

%List of mutiple total signals sent
N_list = [1e6,1e8,1e10];

%Loss
%total array of loss values to iterate over
lossdB = linspace(0,40,21);
transmittance = 10.^(-lossdB/10);

%list of maximal element of loss
lossList = [7,14,21];

%filestring for optimal values
filestrOptVals = "optValsDecoy46_N=";

% store qkdInput for later
qkdInputSave = qkdInput;

for indexSignals = 1:numel(N_list)
    %set transmittance list for current value of N
    transmittanceTemp = transmittance(1:lossList(indexSignals));

    %Load optimal values for current number of signal values
    fileStrTemp = filestrOptVals + sprintf("%.2e",N_list(indexSignals)) +"_q=1.00e+00.csv";
    %optimal values are sorted in coulmns as | logRenyiAlpha | ...
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
        
        % % optimize alpha
        % logrenyiAlpha.lowerBound = -5;
        % logrenyiAlpha.upperBound = -0.5;
        % logrenyiAlpha.initVal = logAlpha;
        % qkdInput.addOptimizeParameter("logrenyiAlpha", logrenyiAlpha);
        
        % run the QKDSolver with this input
        results(indexLoss) = MainIteration(qkdInput);
    end
    %edit qkdinput to save correct results
    qkdInput.addScanParameter("transmittance", num2cell(transmittanceTemp));

    %filestring for saving
    filestr = sprintf("../data/RenyiDecoy46Results_%.2e",N_list(indexSignals)) + "_q=1.00e+00_1decoy.mat";

    % save the results and preset to a file
    %save results
    results = results(1:numel(transmittanceTemp));
    save(filestr,"results","qkdInput");
end

%% plot the result
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")