%% Load Results
%number of loss values used
numLoss = 21;
optVarNames = {'logRenyiAlpha'};

%% epsilonInt = 0
%N=1e6
matN6_eps0 = load("data/RenyiDecoyBB84PassiveResults_1.00e+06_epsInt=0.00e+00.mat");
[keyRatesN6_eps0,optvalsN6_eps0] = parseKeyRatesAndOptVals(matN6_eps0,numLoss,optVarNames);

%N=1e8
matN8_eps0 = load("data/RenyiDecoyBB84PassiveResults_1.00e+08_epsInt=0.00e+00.mat");
[keyRatesN8_eps0,optvalsN8_eps0] = parseKeyRatesAndOptVals(matN8_eps0,numLoss,optVarNames);

%N=1e10
matN10_eps0 = load("data/RenyiDecoyBB84PassiveResults_1.00e+10_epsInt=0.00e+00.mat");
[keyRatesN10_eps0,optvalsN10_eps0] = parseKeyRatesAndOptVals(matN10_eps0,numLoss,optVarNames);

%% epsilonInt = 10%
%N=1e6
matN6_eps10 = load("data/RenyiDecoyBB84PassiveResults_1.00e+06_epsInt=1.00e-01.mat");
[keyRatesN6_eps10, optvalsN6_eps10] = parseKeyRatesAndOptVals(matN6_eps10,numLoss,optVarNames);

%N=1e8
matN8_eps10 = load("data/RenyiDecoyBB84PassiveResults_1.00e+08_epsInt=1.00e-01.mat");
[keyRatesN8_eps10, optvalsN8_eps10] = parseKeyRatesAndOptVals(matN8_eps10,numLoss,optVarNames);

%N=1e10
matN10_eps10 = load("data/RenyiDecoyBB84PassiveResults_1.00e+10_epsInt=1.00e-01.mat");
[keyRatesN10_eps10, optvalsN10_eps10] = parseKeyRatesAndOptVals(matN10_eps10,numLoss,optVarNames);

%% epsilonInt = 25%
%N=1e6
matN6_eps25 = load("data/RenyiDecoyBB84PassiveResults_1.00e+06_epsInt=2.50e-01.mat");
[keyRatesN6_eps25, optvalsN6_eps25] = parseKeyRatesAndOptVals(matN6_eps25,numLoss,optVarNames);

%N=1e8
matN8_eps25 = load("data/RenyiDecoyBB84PassiveResults_1.00e+08_epsInt=2.50e-01.mat");
[keyRatesN8_eps25,optvalsN8_eps25] = parseKeyRatesAndOptVals(matN8_eps25,numLoss,optVarNames);

%N=1e10
matN10_eps25 = load("data/RenyiDecoyBB84PassiveResults_1.00e+10_epsInt=2.50e-01.mat");
[keyRatesN10_eps25,optvalsN10_eps25] = parseKeyRatesAndOptVals(matN10_eps25,numLoss,optVarNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
% Loss etc. for plotting
tempeta = arrayfun(@(x) x.currentParams.transmittance, matN10_eps0.results);
tempptest = arrayfun(@(x) x.currentParams.probTest, matN10_eps0.results);

eta = tempeta(1:end);
etadB = -10*log10(eta);

etaAsymp = 10.^(-linspace(0,max(etadB,[],"all"),100));
etaAsympdB = -10*log10(etaAsymp);

%List of total signals sent
Nlist = 10.^([6,8,10]);

%Color list
colorList = ["#0072BD", "#D95319", "#77AC30"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Asymptotic key rates

results_eps0 = matN10_eps0.results;
fEC      = results_eps0(1).currentParams.fEC;
theta    = results_eps0(1).currentParams.misalignmentAngle;
musig    = results_eps0(1).currentParams.GROUP_decoys_1;
delta    = results_eps0(1).currentParams.GROUP_deltaDecoys_1;
probtest = results_eps0(1).currentParams.probTest;

pzA = 1 - probtest;
pzB = 1 - probtest;
pxB = probtest;

% 1D array instead of 2D matrix since there is only 1 epsilonInt
keyRatesAsymp = zeros(length(etaAsymp), 1);

for index = 1:length(etaAsymp)
    keyRatesAsymp(index) = asymptotic_decoy_BB84_Passive(pzA, pzB, pxB, musig, etaAsymp(index), theta, fEC);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Key rate
%Options for plots 
[x0,y0,width,height] = deal(50,100,600,500);
figure
set(gcf,'position',[x0,y0,width,height])

darkGreen = "#77AC30";

% epsilonInt = 0, i.e. perfect
semilogy(etadB,keyRatesN10_eps0,":p","Color","black","Markersize",10,"HandleVisibility","off")
hold on
semilogy(etadB,keyRatesN8_eps0,"--s","Color","black","HandleVisibility","off")
semilogy(etadB,keyRatesN6_eps0,"-.o","Color","black","HandleVisibility","off")

%epsilonInt = 10%
semilogy(etadB,keyRatesN10_eps10,":p","Color",colorList(1),"Markersize",10,"HandleVisibility","off")
semilogy(etadB,keyRatesN8_eps10,"--s","Color",colorList(1),"Markersize",8,"HandleVisibility","off")
semilogy(etadB,keyRatesN6_eps10,"-.o","Color",colorList(1),"HandleVisibility","off")

%epsilonInt = 25%
semilogy(etadB,keyRatesN10_eps25,":p","Color",colorList(2),"Markersize",10,"HandleVisibility","off")
semilogy(etadB,keyRatesN8_eps25,"--s","Color",colorList(2),"Markersize",8,"HandleVisibility","off")
semilogy(etadB,keyRatesN6_eps25,"-.o","Color",colorList(2),"HandleVisibility","off")
semilogy(etaAsympdB,keyRatesAsymp,"--","Color",darkGreen,"LineWidth",1.5,"HandleVisibility","off")

h_n10 = plot(NaN,NaN,":p","Color","black","Markersize",8,"DisplayName",sprintf("n = 10^{%.0f}",log10(Nlist(3))));
h_n8  = plot(NaN,NaN,"--s","Color","black","Markersize",8,"DisplayName",sprintf("n = 10^{%.0f}",log10(Nlist(2))));
h_n6  = plot(NaN,NaN,"-.o","Color","black","Markersize",8,"DisplayName",sprintf("n = 10^{%.0f}",log10(Nlist(1))));
h_inf = plot(NaN,NaN,"--","Color",darkGreen,"LineWidth",1.5,"DisplayName","Infinite Decoy");

h_e0  = plot(NaN,NaN,"-","Color","black",      "LineWidth",2,"DisplayName","\epsilon_{int} = 0");
h_e10 = plot(NaN,NaN,"-","Color",colorList(1), "LineWidth",2,"DisplayName","\epsilon_{int} = 10%");
h_e25 = plot(NaN,NaN,"-","Color",colorList(2), "LineWidth",2,"DisplayName","\epsilon_{int} = 25%");

lgd = legend([h_n10, h_n8, h_n6, h_inf, h_e0, h_e10, h_e25], 'NumColumns', 2);
lgd.FontSize = 10;
lgd.Location = 'northeast';
lgd.Box = 'on';

xlabel('Loss in dB',FontSize=14)
ylabel('Secret Key Rate',FontSize=14)
ylim([5*1e-6 0.5])
hold off

%save figure
f1=gca;
filestr1 = "PassiveDecoyBB84_2Decoy";
% exportgraphics(f1,filestr1 + ".pdf",'ContentType','vector')
% exportgraphics(f1,filestr1 + ".eps",'ContentType','vector')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rates,optValsTable] = parseKeyRatesAndOptVals(data,numElmts,optValNames)
    %extract results from data
    results = data.results;

    %extract key rates
    listKeyRate = [results(:).keyRate];

    %extract current parameters
    rates = zeros(1,numElmts);
    
    %get size
    numList = numel(listKeyRate);

    %populate with entries in list
    rates(1:numList) = listKeyRate;

    %extract current parameters
    listOptVals = [results(:).currentParams];

    % extract logrenyiAlpha, probTest, signal intensity
    listOptVals = [[listOptVals(:).logrenyiAlpha].'];

    %get size
    numCols = size(listOptVals,2);
    numList = size(listOptVals,1);

    %preallocate optimal values with 0's
    optVals = zeros(numElmts,numCols);

    %populate with entries in list
    optVals(1:numList,:) = listOptVals;

    %convert to table
    optValsTable = array2table(optVals);

    %assign header
    optValsTable.Properties.VariableNames(1:numel(optValNames)) = optValNames;
end

function R = asymptotic_decoy_BB84_Passive(pzA, pzB, pxB, mu, eta, theta, f_EC)
    % 1. Convert Bloch sphere misalignment angle to intrinsic error probability
    e_mis = sin(theta/2)^2;

    % Light intensity arriving at Bob's passive arms
    mu_z = eta .* pzB .* mu;
    mu_x = eta .* pxB .* mu;

    % 2. Strict Single-Click Probabilities in Z
    % Click in correct detector AND no click in wrong detector
    P_correct_only = exp(-mu_z .* e_mis) - exp(-mu_z);

    % Click in wrong detector AND no click in correct detector
    P_wrong_only = exp(-mu_z .* (1 - e_mis)) - exp(-mu_z);

    % Total strict single-click yield in Z (must also have NO clicks in X)
    P_no_X = exp(-mu_x);
    Q_mu_Z = (P_correct_only + P_wrong_only) .* P_no_X;

    % 3. QBER for the strict single-click Z events
    ez = P_wrong_only ./ (P_correct_only + P_wrong_only);

    % Single-photon phase error remains purely optical misalignment
    ex = e_mis; 

    % Safe binary entropy function
    h2 = @(x) -max(x, eps).*log2(max(x, eps)) - max(1-x, eps).*log2(max(1-x, eps));

    % 4. Single-photon detection probability per pulse
    % Single photons inherently only cause single clicks
    Q_1_Z = (mu .* exp(-mu)) .* (eta .* pzB);

    % 5. Asymptotic decoy-state BB84 secret key rate per pulse
    R = pzA .* ( Q_1_Z .* (1 - h2(ex)) - f_EC .* Q_mu_Z .* h2(ez) );

    % Ensure rate doesn't drop below zero
    R = max(R, 0);
end
