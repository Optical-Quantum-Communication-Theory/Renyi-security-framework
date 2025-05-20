%pick preset
qkdInput = RenyiDecoyBB84PassivePreset();

%List of mutiple total signals sent
% N_list = [1e10];
% N_list = [1e8,1e10];
N_list = [1e6,1e8,1e10];

%Loss
%total array of loss values to iterate over
lossdB = linspace(0,40,21);
transmittance = 10.^(-lossdB/10);

%list of maximal element of loss
% lossList = [7,16,21]; %delta = 0
lossList = [6,16,21]; %delta = 1e-6
% lossList = [6,14,15]; %delta = 1e-4
% lossList = [4,4,4]; %delta =1e-2
% lossList = [8,12];


qkdInputSave = qkdInput;

for indexNumsig = 1:numel(N_list)
    %set transmittance list for current value of N
    transmittanceTemp = transmittance(1:lossList(indexNumsig));

    for indexLoss = 1:numel(transmittanceTemp)
        fprintf("Iteration %.0f of %.0f for %.0e",indexLoss, numel(transmittanceTemp),N_list(indexNumsig))
        
        %Number of signals sent
        qkdInput.addFixedParameter("Ntot", N_list(indexNumsig));
        
        %loss
        qkdInput.addFixedParameter("transmittance", transmittanceTemp(indexLoss));
    
        %Renyi param
        % %fixed alpha
        % logAlpha = logAlphaFunc(lossdB(index),Ntot);
        % qkdInput.addFixedParameter("logrenyiAlpha", logAlpha);
    
        % %optimzed alpha
        logrenyiAlpha.lowerBound = -5;
        logrenyiAlpha.upperBound = -0.5;
        logrenyiAlpha.initVal = logAlphaFunc1Decoy(lossdB(indexLoss),N_list(indexNumsig));
        qkdInput.addOptimizeParameter("logrenyiAlpha",logrenyiAlpha);
        
        %run the QKDSolver with this input and store results
        results(indexLoss) = MainIteration(qkdInput);
        qkdInputSave(indexLoss) = qkdInput;
    
        filestr = sprintf("RenyiDecoyBB84PassiveResults_%.2e",N_list(indexNumsig)) ...
            + sprintf("_delta%.2e",1e-6) + "_1decoy_Rerun.mat";
    
        %save results
        save(filestr,"results","qkdInputSave");
    end
    %edit qkdinput
    qkdInput.addScanParameter("transmittance", num2cell(transmittance));
    
    %save results
    save(filestr,"results","qkdInputSave","qkdInput");
end

%% plot the result
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")

function logAlpha = logAlphaFunc1Decoy(transmittance,numSig)
    switch numSig
        case 1e10
            p = [0.0460, -3.1688];
            logAlpha = polyval(p,transmittance);
        case 1e8
            p = [0.0493, -2.9120];
            logAlpha = polyval(p,transmittance);
        case 1e6
            p = [0.1067, -2.2323];
            logAlpha = polyval(p,transmittance);
        otherwise
            disp("numSig does not match")
    end
end