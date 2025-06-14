%pick preset
clear all
qkdInput = RenyiPhaseImpDecoy46LossyPreset();

%List of mutiple total signals sent
N_list = [1e6,1e8,1e10];

%Loss
%total array of loss values to iterate over
lossdB = linspace(0,40,21);
transmittance = 10.^(-lossdB/10);

%list of maximal element of loss with q = 0.99
lossList = [7,12,16]; 

%filestring for optimal values
filestrOptVals = "optValsDecoy46_N=";

% store qkdInput for later
qkdInputSave = qkdInput;

for indexNumsig = 1:numel(N_list)
    %set transmittance list for current value of N
    transmittanceTemp = transmittance(1:lossList(indexNumsig));

    %Load optimal values for current number of signal values
    fileStrTemp = filestrOptVals + sprintf("%.2e",N_list(indexNumsig)) +"_q=9.90e-01.csv";
    %optimal values are sorted in coulmns as | logRenyiAlpha | ...
    optvals = readmatrix(fileStrTemp);

    for indexLoss = 1:numel(transmittanceTemp)
        fprintf("Iteration %.0f of %.0f for %.0e",indexLoss, numel(transmittanceTemp),N_list(indexNumsig))
        
        %Number of signals sent
        qkdInput.addFixedParameter("Ntot", N_list(indexNumsig));
        
        %loss
        qkdInput.addFixedParameter("transmittance", transmittanceTemp(indexLoss));
    
        %Renyi param
        %fixed alpha
        logAlpha = optvals(indexLoss,1);
        qkdInput.addFixedParameter("logrenyiAlpha", logAlpha);
    
        % optimize alpha
        % bndsLogAlpha = lowerUpperBnds_from_optvals(indexLoss,optvals(:,1),-5,-0.5);
        % logrenyiAlpha.lowerBound = bndsLogAlpha(1);
        % logrenyiAlpha.upperBound = bndsLogAlpha(2);
        % logrenyiAlpha.initVal = logAlpha;
        % qkdInput.addOptimizeParameter("logrenyiAlpha", logrenyiAlpha);
        
        %run the QKDSolver with this input and store results
        results(indexLoss) = MainIteration(qkdInput);
        qkdInputSave(indexLoss) = qkdInput;
    end
    %edit qkdinput to save correct results
    qkdInput.addScanParameter("transmittance", num2cell(transmittanceTemp));
    
    %filestring for saving
    filestr = sprintf("../data/RenyiDecoy46Results_%.2e",N_list(indexNumsig)) ...
    + sprintf("_q=%.2e",0.99) + "_1decoy.mat";

    % save the results and preset to a file
    %save results
    results = results(1:numel(transmittanceTemp));
    save(filestr,"results","qkdInput");
end

%% plot the result
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")