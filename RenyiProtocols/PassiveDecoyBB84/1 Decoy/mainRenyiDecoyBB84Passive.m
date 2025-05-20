%pick preset
clear all

qkdInput = RenyiDecoyBB84PassivePreset();

%total number of signals sent
Ntot = 1e8;
qkdInput.addFixedParameter("Ntot", Ntot); 

%Loss
%total array of loss values to iterate over
lossdB = linspace(0,40,21);
transmittance = 10.^(-lossdB/10);

transmittance = transmittance(14);
lossdB = lossdB(14);

% transmittance = transmittance(16);
% lossdB = lossdB(16); %needs to be bigger than 7e-5


%initialize results
% results = struct();
qkdInputSave = qkdInput;

for indexLoss = 1:numel(transmittance)
    fprintf("Iteration %.0f of %.0f for %.1e",indexLoss, numel(transmittance),Ntot  )
    %loss
    qkdInput.addFixedParameter("transmittance", transmittance(indexLoss));

    %Renyi param
    %fixed alpha
    % logAlpha = logAlphaFunc1Decoy(lossdB(indexLoss),Ntot);
    qkdInput.addFixedParameter("logrenyiAlpha", -1.727614109167721); %-1.18 = 2e1-6 %-1.2 = 3 1e-7

    %optimzed alpha
    % logrenyiAlpha.lowerBound = -5;
    % logrenyiAlpha.upperBound = -0.5;
    % logrenyiAlpha.initVal = -1.65;%logAlphaFunc1Decoy(lossdB(indexLoss),Ntot);
    % qkdInput.addOptimizeParameter("logrenyiAlpha",logrenyiAlpha);
    
    %run the QKDSolver with this input and store results
    results(indexLoss) = MainIteration(qkdInput);
    qkdInputSave(indexLoss) = qkdInput;

    filestr = sprintf("RenyiDecoyBB84PassiveResults_%.2e",Ntot) ...
        + sprintf("_delta%.2e",0) + "14_1decoy.mat";

    %save results
    save(filestr,"results","qkdInputSave");
end

%edit qkdinput
qkdInput.addScanParameter("transmittance", num2cell(transmittance));

%save results
save(filestr,"results","qkdInputSave","qkdInput");

%% plot the result
keyRateFixed = arrayfun(@(x) x.debugInfo.keyRateModule.keyRateFixed, results);
fig = axes(figure());
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log","figAxis",fig)

hold(fig,"on")
trans = qkdInput.scanParameters.transmittance;

db = -10*log10(cell2mat(trans));

plot(fig,db,keyRateFixed,"+-")

hold(fig,"off")

function logAlpha = logAlphaFunc2Decoy(transmittancedB,numSig)
    switch numSig
        case 1e10
            % logAlpha = 0.056*transmittancedB -3.5;
            % logAlpha = 0.044*transmittancedB -3.5;
            % logAlpha = 0.063*transmittancedB -3.5;
            
            %linear
            % p = [0.045069089032499,-3.445907184885494];

            % p = [-3.346287989829687e-05,0.001872513491812,0.023077104380528,-3.391031373627249];
            % p = [4.689340987364580e-09,-6.184804552979916e-07,3.180259445765552e-05,-7.976709597381468e-04,0.009787240193045,-0.047392095402089,0.047370561161732,-3.241595818952964];
            % p = [2.336624946723978e-09,-3.120964582220572e-07,1.622244821316323e-05,-4.072192193209328e-04,0.004794432617370,-0.017252365845436,-0.018746568316996,-3.227567432635448];
            % p = [-2.682990367342843e-08,8.859000785921754e-06,-6.287452473535632e-04,0.016010003960142,-0.092409911281620,-3.212349813651729];
            
            %best degree 5
            % p = [-3.279345492347693e-07,3.800060955429918e-05,-0.001572374005146,0.027571697737400,-0.138701792578240,-3.187300573727054];
            
            % logAlpha = polyval(p,transmittancedB);
            
            lossdB = linspace(0,40,21);
            indexList = find(transmittancedB==lossdB);

            % best list
            logAlphaList = [-3.24899291992188,-3.21679687500000,-3.50195312500000,-3.24999237060547,-3.09375000000000,-3.37500000000000,-2.87500000000000,-2.74218750000000,-2.34765625000000,-2.31250000000000,-2.78131103515625,-2.62989098223079,-2.42853854294223,-2.33810610284449,-1.75616754824392,-2.17580000000000,-2.12949875503764,-2.00185179631295,-1.93386821150768,-1.50573414892714,-1.59882275851676];
            logAlpha = logAlphaList(indexList);
            disp(logAlpha)
        case 1e8
            logAlpha = 0.05*transmittancedB -3.03;
        case 1e6
            logAlpha = 0.047*transmittancedB -2.18;
        otherwise
            disp("numSig does not match")
    end
end

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