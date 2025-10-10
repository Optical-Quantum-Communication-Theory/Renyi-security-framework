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
%% Key rate
%Options for plots 
[x0,y0,width,height] = deal(50,100,600,500);


figure
set(gcf,'position',[x0,y0,width,height])

%epsilonInt = 0, i.e. perfect
semilogy(etadB,keyRatesN10_eps0,":p","Color","black",'LineWidth',1,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(3))))

hold on
semilogy(etadB,keyRatesN8_eps0,":s","Color","black",'LineWidth',1,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(2))))
semilogy(etadB,keyRatesN6_eps0,":o","Color","black",'LineWidth',1,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(1))))

%epsilonInt = 10%
semilogy(etadB,keyRatesN10_eps10,"-p","Color",colorList(2),"Markersize",10, "DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(3))))
semilogy(etadB,keyRatesN8_eps10,"--s","Color",colorList(2),"Markersize",8, "DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(2))))
semilogy(etadB,keyRatesN6_eps10,"-.o","Color",colorList(2),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(1))))

%epsilonInt = 25%
semilogy(etadB,keyRatesN10_eps25,"-p","Color",colorList(3),"Markersize",10,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(3))))
semilogy(etadB,keyRatesN8_eps25,"--s","Color",colorList(3),"Markersize",8,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(2))))
semilogy(etadB,keyRatesN6_eps25,"-.o","Color",colorList(3),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(1))))

lgd = legend('NumColumns',4);
lgd.FontSize = 10;
lgd.Location = 'northeast';
xlabel('Loss in dB',FontSize=14)
ylabel('Secret key rate',FontSize=14)
ylim([1*1e-6 5])

titles = {'\epsilon_{int} = 0','\epsilon_{int} = 10%', '\epsilon_{int} = 25%'};
titles = sprintf('%-30s', titles{:});
lgd.Title.String = titles;
hold off

%save figure
f1=gca;
filestr1 = "PassiveDecoyBB84_2Decoy.pdf";
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