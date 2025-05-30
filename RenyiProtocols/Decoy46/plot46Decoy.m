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
matN8_q099 = load("data/RenyiDecoy46Results_1.00e+08_q=9.90e-01_1decoy_Rerun.mat");
[keyRatesN8_q099, optvalsN8_q099] = parseKeyRatesAndOptVals(matN8_q099,numLoss,optVarNames);

%N=1e10
% matN10_q099 = load("data/RenyiDecoy46Results_1.00e+10_q=9.90e-01_1decoy_Rerun.mat");
matN10_q099 = load("data/decoy46Results.mat");
[keyRatesN10_q099, optvalsN10_q099] = parseKeyRatesAndOptVals(matN10_q099,numLoss,optVarNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
% Loss etc. for plotting
tempeta = arrayfun(@(x) x.currentParams.transmittance, matN10_q1.results);
tempptest = arrayfun(@(x) x.currentParams.probTest, matN10_q1.results);

eta = tempeta(1:end);
etadB = -10*log10(eta);

etaAsymp = 10.^(-linspace(0,max(etadB,[],"all"),100));
etaAsympdB = -10*log10(etaAsymp);

%List of total signals sent
Nlist = 10.^([6,8,10,12,14]);

%Color list
colorList = ["#0072BD", "#D95319", "#77AC30"];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Options for plots 
[x0,y0,width,height] = deal(50,100,600,500);


figure
set(gcf,'position',[x0,y0,width,height])

%PS with q=1
semilogy(etadB,keyRatesN14PS,":p","Color","black",'LineWidth',1,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(5))))

hold on
semilogy(etadB,keyRatesN12PS,":s","Color","black",'LineWidth',1,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(4))))
semilogy(etadB,keyRatesN10PS,":o","Color","black",'LineWidth',1,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(3))))
% semilogy(etadB,keyRatesN8PS,":o","Color","black",'LineWidth',1,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(2))))

%GREAT q=1
semilogy(etadB,keyRatesN10_q1,"-p","Color",colorList(1),"Markersize",10,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(3))))
semilogy(etadB,keyRatesN8_q1,"--s","Color",colorList(1),"Markersize",8, "DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(2))))
semilogy(etadB,keyRatesN6_q1,"-.o","Color",colorList(1),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(1))))

%GREAT q=0.99
semilogy(etadB,keyRatesN10_q099,"-p","Color",colorList(2),"Markersize",10, "DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(3))))
semilogy(etadB,keyRatesN8_q099,"--s","Color",colorList(2),"Markersize",8, "DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(2))))
semilogy(etadB,keyRatesN6_q099,"-.o","Color",colorList(2),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(1))))

lgd = legend('NumColumns',3);
lgd.FontSize = 10;
lgd.Location = 'northeast';
xlabel('transmittance in dB',FontSize=14)
ylabel('Secret key rate',FontSize=14)
ylim([1*1e-6 2])

titles = {'PS q=1','q=1', 'q=0.99'};
titles = sprintf('%-25s', titles{:});
lgd.Title.String = titles;
hold off

%save figure
f1=gca;
filestr1 = "Decoy46.pdf";
% exportgraphics(f1,filestr1,'ContentType','vector')

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