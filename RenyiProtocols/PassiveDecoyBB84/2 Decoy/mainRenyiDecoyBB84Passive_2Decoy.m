clear all

%List of mutiple total signals sent
% N_list = [1e6,1e8,1e10];
N_list = [1e10];

%Loss
%total array of loss values to iterate over
lossdB = linspace(0,40,21);
transmittance = 10.^(-lossdB/10);

%list of maximal element of loss ordered epsilon_int = 0, 10%, 25%
%(named delta in description)
lossList = [[7,17,21]; [5,16,21] ;[4,13,21]];

% list of epsilon_int values
% epsilonInt_List = [0, 0.1, 0.25];
epsilonInt_List = [0.1];

%filestring for optimal values
filestrOptVals = "optimalValues\optValsPassiveDecoyBB84_N=";


for indexEps = 1:numel(epsilonInt_List)
    %pick preset (depends on deviation eps_int
    qkdInput = RenyiDecoyBB84PassivePreset_2Decoy(epsilonInt_List(indexEps));

    % store qkdInput for later
    qkdInputSave = qkdInput;

    for indexSignals = 1:numel(N_list)
        %set transmittance list for current value of N
        transmittanceTemp = transmittance(19:21);%transmittance(1:lossList(indexEps,indexSignals));
    
        %Load optimal values for current number of signal values
        fileStrTemp = filestrOptVals + sprintf("%.2e",N_list(indexSignals)) ...
            + sprintf("_epsilonInt=%.2e",epsilonInt_List(indexEps)) +".csv";
        %optimal values are sorted in coulmns as | logRenyiAlpha| ...
        optvals = readmatrix(fileStrTemp);
    
        for indexLoss = 19:21%1:numel(transmittanceTemp)
            fprintf("Iteration %.0f of %.0f for N=%.0e with eps_int=%.2e",indexLoss,...
                numel(transmittanceTemp),N_list(indexSignals),epsilonInt_List(indexEps))
    
            %Add total signals sent from list above
            qkdInput.addFixedParameter("Ntot",N_list(indexSignals))
        
            %Add loss until element from list above
            qkdInput.addScanParameter("transmittance", num2cell(transmittance(indexLoss))); 
        
            %Add Renyi param from optimal values
            % fixed alpha
            logAlpha = -2;%optvals(indexLoss,1);
            % qkdInput.addFixedParameter("logrenyiAlpha", logAlpha);
            
            % optimize alpha
            bndsLogAlpha = lowerUpperBnds_from_optvals(indexLoss,optvals(:,1),-2.5,-1.5);
            logrenyiAlpha.lowerBound = bndsLogAlpha(1);
            logrenyiAlpha.upperBound = bndsLogAlpha(2);
            logrenyiAlpha.initVal = logAlpha;
            qkdInput.addOptimizeParameter("logrenyiAlpha", logrenyiAlpha);
    
            % run the QKDSolver with this input
            results(indexLoss) = MainIteration(qkdInput);
        end
        % edit qkdinput to save correct results
        qkdInput.addScanParameter("transmittance", num2cell(transmittanceTemp));

        %filestring for saving
        filestr = sprintf("data/RenyiDecoyBB84PassiveResults_%.2e",N_list(indexSignals)) + ...
            sprintf("_epsInt=%.2e",epsilonInt_List(indexEps)) + "_new_tight.mat";

        % save the results and preset to a file
        %save results
        results = results(19:21);%results(1:numel(transmittanceTemp));
        save(filestr,"results","qkdInput");
    end
end

%% plot the last result
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")