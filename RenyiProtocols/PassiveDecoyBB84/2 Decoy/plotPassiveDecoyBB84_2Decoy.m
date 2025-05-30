%% Load Results
%number of loss values used
numLoss = 21;
optVarNames = {'logRenyiAlpha'};

%% delta = 0
%N=1e6
matN6d0 = load("data/RenyiDecoyBB84PassiveResults_1.00e+06_epsInt=0.00e+00.mat");
[keyRatesN6d0,optvalsN6d0] = parseKeyRatesAndOptVals(matN6d0,numLoss,optVarNames);

%N=1e8
matN8d0 = load("data/RenyiDecoyBB84PassiveResults_1.00e+08_epsInt=0.00e+00.mat");
[keyRatesN8d0,optvalsN8d0] = parseKeyRatesAndOptVals(matN8d0,numLoss,optVarNames);

%N=1e10
matN10d0 = load("data/RenyiDecoyBB84PassiveResults_1.00e+10_epsInt=0.00e+00.mat");
[keyRatesN10d0,optvalsN10d0] = parseKeyRatesAndOptVals(matN10d0,numLoss,optVarNames);

%% delta = 1e-6
%N=1e6
matN6d6 = load("data/RenyiDecoyBB84PassiveResults_1.00e+06_epsInt=1.00e-06.mat");
[keyRatesN6d6,optvalsN6d6] = parseKeyRatesAndOptVals(matN6d6,numLoss,optVarNames);

%N=1e8
matN8d6 = load("data/RenyiDecoyBB84PassiveResults_1.00e+08_epsInt=1.00e-06.mat");
[keyRatesN8d6,optvalsN8d6] = parseKeyRatesAndOptVals(matN8d6,numLoss,optVarNames);

%N=1e10
matN10d6 = load("data/RenyiDecoyBB84PassiveResults_1.00e+10_epsInt=1.00e-06.mat");
[keyRatesN10d6, optvalsN10d6] = parseKeyRatesAndOptVals(matN10d6,numLoss,optVarNames);

%% delta = 1e-4
matN6d4 = load("data/RenyiDecoyBB84PassiveResults_1.00e+06_epsInt=1.00e-04.mat");
[keyRatesN6d4, optvalsN6d4] = parseKeyRatesAndOptVals(matN6d4,numLoss,optVarNames);

matN8d4 = load("data/RenyiDecoyBB84PassiveResults_1.00e+08_epsInt=1.00e-04.mat");
[keyRatesN8d4, optvalsN8d4] = parseKeyRatesAndOptVals(matN8d4,numLoss,optVarNames);

matN10d4 = load("data/RenyiDecoyBB84PassiveResults_1.00e+10_epsInt=1.00e-04.mat");
[keyRatesN10d4, optvalsN10d4] = parseKeyRatesAndOptVals(matN10d4,numLoss,optVarNames);

%% delat = 1e-2
matN6d2 = load("data/RenyiDecoyBB84PassiveResults_1.00e+06_epsInt=1.00e-02.mat");
[keyRatesN6d2, optvalsN6d2] = parseKeyRatesAndOptVals(matN6d2,numLoss,optVarNames);

matN8d2 = load("data/RenyiDecoyBB84PassiveResults_1.00e+08_epsInt=1.00e-02.mat");
[keyRatesN8d2,optvalsN8d2] = parseKeyRatesAndOptVals(matN8d2,numLoss,optVarNames);

matN10d2 = load("data/RenyiDecoyBB84PassiveResults_1.00e+10_epsInt=1.00e-02.mat");
[keyRatesN10d2,optvalsN10d2] = parseKeyRatesAndOptVals(matN10d2,numLoss,optVarNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
% Loss etc. for plotting
tempeta = arrayfun(@(x) x.currentParams.transmittance, matN10d0.results);
tempptest = arrayfun(@(x) x.currentParams.probTest, matN10d0.results);

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

%delta = 0, i.e. perfect
semilogy(etadB,keyRatesN10d0,":p","Color","black",'LineWidth',1,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(3))))

hold on
semilogy(etadB,keyRatesN8d0,":s","Color","black",'LineWidth',1,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(2))))
semilogy(etadB,keyRatesN6d0,":o","Color","black",'LineWidth',1,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(1))))

%delta = 1e-6
semilogy(etadB,keyRatesN10d6,"-p","Color",colorList(1),"Markersize",10,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(3))))
semilogy(etadB,keyRatesN8d6,"--s","Color",colorList(1),"Markersize",8, "DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(2))))
semilogy(etadB,keyRatesN6d6,"-.o","Color",colorList(1),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(1))))

%delta = 1e-4
semilogy(etadB,keyRatesN10d4,"-p","Color",colorList(2),"Markersize",10, "DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(3))))
semilogy(etadB,keyRatesN8d4,"--s","Color",colorList(2),"Markersize",8, "DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(2))))
semilogy(etadB,keyRatesN6d4,"-.o","Color",colorList(2),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(1))))

%delta = 1e-2
semilogy(etadB,keyRatesN10d2,"-p","Color",colorList(3),"Markersize",10,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(3))))
semilogy(etadB,keyRatesN8d2,"--s","Color",colorList(3),"Markersize",8,"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(2))))
semilogy(etadB,keyRatesN6d2,"-.o","Color",colorList(3),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(1))))

lgd = legend('NumColumns',4);
lgd.FontSize = 10;
lgd.Location = 'northeast';
xlabel('transmittance in dB',FontSize=14)
ylabel('Secret key rate',FontSize=14)
ylim([1*1e-6 5])

titles = {'\epsilon_{int} = 0','\epsilon_{int} = 10^{-6}', '\epsilon_{int} = 10^{-4}', '\epsilon_{int} = 10^{-2}'};
titles = sprintf('%-30s', titles{:});
lgd.Title.String = titles;
hold off

%save figure
f1=gca;
filestr1 = "PassiveDecoyBB84.pdf";
exportgraphics(f1,filestr1,'ContentType','vector')

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