%pick preset
qkdInput = RenyiPhaseImpDecoy46LossyPreset();

%List of mutiple total signals sent
N_list = [1e10];
% N_list = [1e8,1e10];
% N_list = [1e6,1e8,1e10];

%Loss
%total array of loss values to iterate over
lossdB = linspace(0,40,21);
transmittance = 10.^(-lossdB/10);

%list of maximal element of loss
% lossList = [7,15,18]; %q = 1
% lossList = [7,15,18]; %q = 0.99
lossList = [18]; %q = 0.99

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
        % logAlpha = logAlphaFunc1Decoy46(lossdB(indexLoss),N_list(indexNumsig));
        % qkdInput.addFixedParameter("logrenyiAlpha", logAlpha);
    
        % %optimzed alpha
        logrenyiAlpha.lowerBound = -5;
        logrenyiAlpha.upperBound = -0.5;
        logrenyiAlpha.initVal = logAlphaFunc1Decoy46(lossdB(indexLoss),N_list(indexNumsig));
        qkdInput.addOptimizeParameter("logrenyiAlpha",logrenyiAlpha);
        
        %run the QKDSolver with this input and store results
        results(indexLoss) = MainIteration(qkdInput);
        qkdInputSave(indexLoss) = qkdInput;
    
        filestr = sprintf("RenyiDecoy46Results_%.2e",N_list(indexNumsig)) ...
            + sprintf("_q=%.2e",0.99) + "_1decoy.mat";
    
        %save results
        save(filestr,"results","qkdInputSave");
    end
    %edit qkdinput
    qkdInput.addScanParameter("transmittance", num2cell(transmittance));
    
    %save results
    results = results(1:numel(transmittanceTemp));
    save(filestr,"results","qkdInputSave","qkdInput");
end

%% plot the result
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")

function logAlpha = logAlphaFunc1Decoy46(transmittance,numSig)
    switch numSig
        case 1e10
            p = [0.041094329702541,-2.338102588945663];
            logAlpha = polyval(p,transmittance);
        case 1e8
            p = [0.034248368091830,-2.201454968462765];
            logAlpha = polyval(p,transmittance);
        case 1e6
            p = [0.042986307705319,-1.498134304600438];
            logAlpha = polyval(p,transmittance);
        otherwise
            disp("numSig does not match")
    end
end