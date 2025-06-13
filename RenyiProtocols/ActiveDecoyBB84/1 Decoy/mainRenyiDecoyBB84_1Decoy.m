%pick preset
clear all
qkdInput = RenyiDecoyBB84ActiveLossyPreset_1Decoy();

%List of mutiple total signals sent
N_list = [1e5,1e6,1e7,1e8,1e9,1e10,1e11];

%Loss
%total array of loss values to iterate over
lossdB = linspace(0,50,26);
transmittance = 10.^(-lossdB/10);

%list of maximal element of loss
lossList = [2,8,10,12,17,19,24];

%filestring for optimal values
filestrOptVals = "optimalValues\optValsActiveDecoyBB84_1_decoy_N=";

% store qkdInput for later
qkdInputSave = qkdInput;

for indexSignals = 1:numel(N_list)
    %set transmittance list for current value of N
    transmittanceTemp = transmittance(1:lossList(indexSignals));

    %Load optimal values for current number of signal values
    fileStrTemp = filestrOptVals + sprintf("%.2e",N_list(indexSignals)) +".csv";
    %optimal values are sorted in coulmns as | logRenyiAlpha | probTest |
    %intensity_1 | ...
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
        
        % optimize alpha
        % bndsLogAlpha = lowerUpperBnds_from_optvals(indexLoss,optvals(:,1),-5,-0.5);
        % logrenyiAlpha.lowerBound = bndsLogAlpha(1);
        % logrenyiAlpha.upperBound = bndsLogAlpha(2);
        % logrenyiAlpha.initVal = logAlpha;
        % qkdInput.addOptimizeParameter("logrenyiAlpha", logrenyiAlpha);

        %Add probTest from optimal values
        %fixed probTest
        pTest = optvals(indexLoss,2);
        qkdInput.addFixedParameter("probTest", pTest);

        % %optimize probTest
        % bndsProbTest = lowerUpperBnds_from_optvals(indexLoss,optvals(:,2),0.01,0.99);
        % probTest.lowerBound = bndsProbTest(1);
        % probTest.upperBound = bndsProbTest(2);
        % probTest.initVal = pTest;
        % qkdInput.addOptimizeParameter("probTest", probTest);
        
        %Add signal intensity from optimal values
        % %fixed signal intensity
        signalIntensity = optvals(indexLoss,3);
        qkdInput.addFixedParameter("GROUP_decoys_1", signalIntensity); %signal intensity
        
        % %Optimize signal intensity
        % bndsSignal = lowerUpperBnds_from_optvals(indexLoss,optvals(:,3),0.001,1);
        % qkdInput.addOptimizeParameter("GROUP_decoys_1", struct("lowerBound",bndsSignal(1),"initVal",signalIntensity,"upperBound",bndsSignal(2))); %signal intensity

        % run the QKDSolver with this input
        results(indexLoss) = MainIteration(qkdInput);
    end
    %edit qkdinput to save correct results
    qkdInput.addScanParameter("transmittance", num2cell(transmittanceTemp));

    %filestring for saving
    filestr = sprintf("data/RenyiDecoyBB84LossyResults_%.2e",N_list(indexSignals)) + "_1decoy.mat";

    % save the results and preset to a file
    results = results(1:numel(transmittanceTemp));
    save(filestr,"results","qkdInput");
end

%% plot the result
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")