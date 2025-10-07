clear all

%List of mutiple total signals sent
N_list = [1e6,1e8,1e10];

%Loss
%total array of loss values to iterate over
lossdB = linspace(0,40,21);
transmittance = 10.^(-lossdB/10);

%list of maximal element of loss ordered epsilon_int = 0, 10%, 25%
%(named delta in description)
lossList = [[7,16,21]; [6,15,21] ;[4,13,17]];

% list of epsilon_int values
epsilonInt_List = [0, 0.1, 0.25];

%filestring for optimal values
filestrOptVals = "optimalValues\optValsPassiveDecoyBB84_2_decoy_N=";

for indexEps = 1:numel(epsilonInt_List)
    %pick preset (depends on deviation eps_int
    qkdInput = RenyiDecoyBB84PassivePreset_2Decoy(epsilonInt_List(indexEps));

    % store qkdInput for later
    qkdInputSave = qkdInput;

    for indexSignals = 1:numel(N_list)
        %set transmittance list for current value of N
        transmittanceTemp = transmittance(1:lossList(indexEps,indexSignals));
    
        %Load optimal values for current number of signal values
        fileStrTemp = filestrOptVals + sprintf("%.2e",N_list(indexSignals)) ...
            + sprintf("_epsilonInt=%.2e",epsilonInt_List(indexEps)) +".csv";
        %optimal values are sorted in coulmns as | logRenyiAlpha| ...
        optvals = readmatrix(fileStrTemp);
    
        for indexLoss = 1:numel(transmittanceTemp)
            fprintf("Iteration %.0f of %.0f for N=%.0e with eps_int=%.2e \n",indexLoss,...
                numel(transmittanceTemp),N_list(indexSignals),epsilonInt_List(indexEps))
    
            %Add total signals sent from list above
            qkdInput.addFixedParameter("Ntot",N_list(indexSignals))
        
            %Add loss until element from list above
            qkdInput.addScanParameter("transmittance", num2cell(transmittance(indexLoss))); 
        
            %Add Renyi param from optimal values
            % fixed alpha
            logAlpha = optvals(indexLoss,1);
            qkdInput.addFixedParameter("logrenyiAlpha", logAlpha);
            
            % optimize alpha
            % bndsLogAlpha = lowerUpperBnds_from_optvals(indexLoss,optvals(:,1),-4.5,-0.8);
            % logrenyiAlpha.lowerBound = bndsLogAlpha(1);
            % logrenyiAlpha.upperBound = bndsLogAlpha(2);
            % logrenyiAlpha.initVal = logAlpha;
            % qkdInput.addOptimizeParameter("logrenyiAlpha", logrenyiAlpha);
    
            % run the QKDSolver with this input
            results(indexLoss) = MainIteration(qkdInput);
        end
        % edit qkdinput to save correct results
        qkdInput.addScanParameter("transmittance", num2cell(transmittanceTemp));

        %filestring for saving
        filestr = sprintf("data/RenyiDecoyBB84PassiveResults_%.2e",N_list(indexSignals)) + ...
            sprintf("_epsInt=%.2e",epsilonInt_List(indexEps)) + ".mat";

        % save the results and preset to a file
        %save results
        results = results(1:numel(transmittanceTemp));
        save(filestr,"results","qkdInput");
    end
end

%% plot the last result
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")