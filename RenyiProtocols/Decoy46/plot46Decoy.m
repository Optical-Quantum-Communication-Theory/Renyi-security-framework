%% Results GREAT
%number of loss values used
numLoss = 21;
optVarNames = {'logRenyiAlpha'};

%% Post-selection q=1

%N=1e8
matN8PS  = readmatrix("PS Data/AdaptiveDecoy46Protocol_PS_1.00e+08.csv");
keyRatesN8PS  = matN8PS(:,2);

%N=1e10
matN10PS  = readmatrix("PS Data/AdaptiveDecoy46Protocol_PS_1.00e+10.csv");
keyRatesN10PS  = matN10PS(:,2);

%N=1e12
matN12PS  = readmatrix("PS Data/AdaptiveDecoy46Protocol_PS_1.00e+12.csv");
keyRatesN12PS  = matN12PS(:,2);

%N=1e14
matN14PS  = readmatrix("PS Data/AdaptiveDecoy46Protocol_PS_1.00e+14.csv");
keyRatesN14PS  = matN14PS(:,2);

%% GREAT q=1
%N=1e6
matN6_q1  = load("data/RenyiDecoy46Results_1.00e+06_q=1.00e+00_1decoy.mat");
[keyRatesN6_q1, optvalsN6_q1] = parseKeyRatesAndOptVals(matN6_q1,numLoss,optVarNames);

%N=1e8
matN8_q1 = load("data/RenyiDecoy46Results_1.00e+08_q=1.00e+00_1decoy.mat");
[keyRatesN8_q1, optvalsN8_q1] = parseKeyRatesAndOptVals(matN8_q1,numLoss,optVarNames);

%N=1e10
matN10_q1 = load("data/RenyiDecoy46Results_1.00e+10_q=1.00e+00_1decoy.mat");
[keyRatesN10_q1, optvalsN10_q1] = parseKeyRatesAndOptVals(matN10_q1,numLoss,optVarNames);

%% GREAT q=0.99
%N=1e6
matN6_q099  = load("data/RenyiDecoy46Results_1.00e+06_q=9.90e-01_1decoy.mat");
[keyRatesN6_q099, optvalsN6_q099] = parseKeyRatesAndOptVals(matN6_q099 ,numLoss,optVarNames);

%N=1e8
matN8_q099 = load("data/RenyiDecoy46Results_1.00e+08_q=9.90e-01_1decoy.mat");
[keyRatesN8_q099, optvalsN8_q099] = parseKeyRatesAndOptVals(matN8_q099,numLoss,optVarNames);

%N=1e10
matN10_q099 = load("data/RenyiDecoy46Results_1.00e+10_q=9.90e-01_1decoy.mat");
[keyRatesN10_q099, optvalsN10_q099] = parseKeyRatesAndOptVals(matN10_q099,numLoss,optVarNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
% Loss etc. for plotting
tempeta = arrayfun(@(x) x.currentParams.transmittance, matN10_q1.results);
tempptest = arrayfun(@(x) x.currentParams.probTest, matN10_q1.results);

eta = tempeta(1:end);
etadB = -10*log10(eta);

etaAsymp = 10.^(-linspace(0,max(etadB,[],"all")/10,100));
etaAsympdB = -10*log10(etaAsymp);

%List of total signals sent
Nlist = 10.^([6,8,10,12,14]);

%Color list
colorList = ["#0072BD", "#D95319", "#77AC30"];
darkGreen = "#77AC30";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Asymptotic key rates
resultsN10 = matN10_q1.results;
fEC      = resultsN10(1).currentParams.fEC;
theta    = resultsN10(1).currentParams.misalignmentAngle;
musig    = resultsN10(1).currentParams.GROUP_decoys_1;
probtest = resultsN10(1).currentParams.probTest;

pzA = 1 - probtest;
pxA = probtest;
probs = resultsN10(1).currentParams.probsB;
pzB = probs(1);
pxB = probs(2);
pyB = probs(3);

keyRatesAsymp = arrayfun(@(x) asymptotic_decoy_4_6_Passive(pzA, pzB, pxB, pyB, musig, x, theta, fEC), etaAsymp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Options for plots 
[x0,y0,width,height] = deal(50,100,600,500);
figure
set(gcf,'position',[x0,y0,width,height])
%PS with q=1
semilogy(etadB,keyRatesN14PS,":p","Color","black",'LineWidth',1,"HandleVisibility","off")
hold on
semilogy(etadB,keyRatesN12PS,":s","Color","black",'LineWidth',1,"HandleVisibility","off")
semilogy(etadB,keyRatesN10PS,":o","Color","black",'LineWidth',1,"HandleVisibility","off")
% semilogy(etadB,keyRatesN8PS,":o","Color","black",'LineWidth',1,"DisplayName", sprintf("n =  10^{%.0f}",log10(Nlist(2))))
%GREAT q=1
semilogy(etadB,keyRatesN10_q1,"-p","Color",colorList(1),"Markersize",10,"HandleVisibility","off")
semilogy(etadB,keyRatesN8_q1,"--s","Color",colorList(1),"Markersize",8,"HandleVisibility","off")
semilogy(etadB,keyRatesN6_q1,"-.o","Color",colorList(1),"HandleVisibility","off")
%GREAT q=0.99
semilogy(etadB,keyRatesN10_q099,"-p","Color",colorList(2),"Markersize",10,"HandleVisibility","off")
semilogy(etadB,keyRatesN8_q099,"--s","Color",colorList(2),"Markersize",8,"HandleVisibility","off")
semilogy(etadB,keyRatesN6_q099,"-.o","Color",colorList(2),"HandleVisibility","off")
% Asymptotic infionite decoy
semilogy(etaAsympdB,keyRatesAsymp,"--","Color",darkGreen,"LineWidth",1.5,"HandleVisibility","off")

% Explicit handles for Column 1 (PS q=1)
h_ps14   = plot(NaN,NaN,":p","Color","black",'LineWidth',1,"Markersize",8,"DisplayName",sprintf("n =  10^{%.0f}",log10(Nlist(5))));
h_ps12   = plot(NaN,NaN,":s","Color","black",'LineWidth',1,"Markersize",8,"DisplayName",sprintf("n =  10^{%.0f}",log10(Nlist(4))));
h_ps10   = plot(NaN,NaN,":o","Color","black",'LineWidth',1,"Markersize",8,"DisplayName",sprintf("n =  10^{%.0f}",log10(Nlist(3))));
h_blank1 = plot(NaN,NaN,"LineStyle","none","Marker","none","DisplayName",""); % Fixed empty placeholder

% Explicit handles for Column 2 (GREAT q=1 + Infinite Decoy)
h_q1_10  = plot(NaN,NaN,"-p","Color",colorList(1),"LineWidth",1,"Markersize",10,"DisplayName",sprintf("n =  10^{%.0f}",log10(Nlist(3))));
h_q1_8   = plot(NaN,NaN,"--s","Color",colorList(1),"LineWidth",1,"Markersize",8,"DisplayName",sprintf("n =  10^{%.0f}",log10(Nlist(2))));
h_q1_6   = plot(NaN,NaN,"-.o","Color",colorList(1),"LineWidth",1,"Markersize",8,"DisplayName",sprintf("n =  10^{%.0f}",log10(Nlist(1))));
h_inf    = plot(NaN,NaN,"--","Color",darkGreen,"LineWidth",1.5,"DisplayName","Infinite Decoy");

% Explicit handles for Column 3 (GREAT q=0.99)
h_q99_10 = plot(NaN,NaN,"-p","Color",colorList(2),"LineWidth",1,"Markersize",10,"DisplayName",sprintf("n =  10^{%.0f}",log10(Nlist(3))));
h_q99_8  = plot(NaN,NaN,"--s","Color",colorList(2),"LineWidth",1,"Markersize",8,"DisplayName",sprintf("n =  10^{%.0f}",log10(Nlist(2))));
h_q99_6  = plot(NaN,NaN,"-.o","Color",colorList(2),"LineWidth",1,"Markersize",8,"DisplayName",sprintf("n =  10^{%.0f}",log10(Nlist(1))));

% Order handles down columns: [Col1_row1, Col1_row2... Col2_row1... Col3_row1...]
lgd = legend([h_ps14, h_ps12, h_ps10, h_blank1, ...
              h_q1_10, h_q1_8, h_q1_6, h_inf, ...
              h_q99_10, h_q99_8, h_q99_6], 'NumColumns', 3);

lgd.FontSize = 10;
lgd.Location = 'northeast';
lgd.Box = 'on';

xlabel('Loss in dB',FontSize=14)
ylabel('Secret key rate',FontSize=14)
ylim([4*1e-6 0.3])
titles = {'PS q=1','q=1', 'q=0.99'};
titles = sprintf('%-25s', titles{:});
lgd.Title.String = titles;
hold off

%save figure
f1=gca;
filestr1 = "Decoy46";
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

function R = asymptotic_decoy_4_6_Passive(pzA, pzB, pxB, pyB, mu, eta, theta, f_EC)
    % ASYMPTOTIC_DECOY_4_6_PASSIVE Calculates secret key rate for a 4-6 protocol
    % (Alice 4 states / 2 bases, Bob passive 6-state receiver) with infinite decoys.
    %
    % Inputs:
    %   pzA           - Alice's active Z-basis probability
    %   pzB, pxB, pyB - Bob's passive beam-splitter ratios (pzB + pxB + pyB = 1)
    %   mu            - Nominal signal state mean photon number
    %   eta           - Overall channel + receiver transmittance
    %   theta         - Optical misalignment angle on Bloch sphere (in radians)
    %   f_EC          - Error correction efficiency factor

    % 1. Convert Bloch sphere misalignment angle to intrinsic error probability
    e_mis = sin(theta/2)^2;
    
    % Light intensity arriving at Bob's three passive arms
    mu_z = eta .* pzB .* mu;
    mu_x = eta .* pxB .* mu;
    mu_y = eta .* pyB .* mu;
    
    % 2. Strict Single-Click Probabilities in Z
    % Click in correct detector AND no click in wrong detector
    P_correct_only = exp(-mu_z .* e_mis) - exp(-mu_z);
    
    % Click in wrong detector AND no click in correct detector
    P_wrong_only = exp(-mu_z .* (1 - e_mis)) - exp(-mu_z);
    
    % Strict single-click yield in Z requires ZERO clicks across BOTH X and Y arms
    P_no_X = exp(-mu_x);
    P_no_Y = exp(-mu_y);
    Q_mu_Z = (P_correct_only + P_wrong_only) .* P_no_X .* P_no_Y;
    
    % 3. QBER for the strict single-click Z events
    ez = P_wrong_only ./ (P_correct_only + P_wrong_only);
    
    % Single-photon error rates across all three bases (pure optical misalignment)
    ex = e_mis; 
    ey = e_mis;
    
    % 4. Entropy Functions
    % Safe binary entropy function (for classical Z-basis error correction)
    h2 = @(x) -max(x, eps).*log2(max(x, eps)) - max(1-x, eps).*log2(max(1-x, eps));
    
    % Three-parameter Six-State Protocol (SSP) entropy function (for privacy amplification)
    function H = h_ssp(ex_val, ey_val, ez_val)
        p0 = 1 - (ex_val + ey_val + ez_val)/2;
        p1 = (-ex_val + ey_val + ez_val)/2;
        p2 = ( ex_val - ey_val + ez_val)/2;
        p3 = ( ex_val + ey_val - ez_val)/2;
        
        p = [p0, p1, p2, p3];
        p = p(p > eps); % Filter zeros safely
        H = -sum(p .* log2(p));
    end
    
    % 5. Single-photon detection probability per pulse (in Z arm)
    % Single photons inherently only cause single clicks
    Q_1_Z = (mu .* exp(-mu)) .* (eta .* pzB);
    
    % 6. Asymptotic Decoy-State 4-6 secret key rate per pulse
    % Privacy amplification uses Six-State entropy: (1 - h_SSP(ex, ey, ez))
    % Classical reconciliation uses Z-basis binary entropy: f_EC * Q_mu_Z * h2(ez)
    R = pzA .* ( Q_1_Z .* (1 - h_ssp(ex, ey, ez)) - f_EC .* Q_mu_Z .* h2(ez) );
    
    % Ensure rate doesn't drop below zero
    R = max(R, 0);
end