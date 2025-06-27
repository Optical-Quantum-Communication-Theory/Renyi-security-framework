%% Results GREAT
matN5 = load("data/RenyiDecoyBB84LossyResults_1.00e+05_1decoy.mat");
matN6 = load("data/RenyiDecoyBB84LossyResults_1.00e+06_1decoy.mat"); 
matN7 = load("data/RenyiDecoyBB84LossyResults_1.00e+07_1decoy.mat");
matN8 = load("data/RenyiDecoyBB84LossyResults_1.00e+08_1decoy.mat");
matN9 = load("data/RenyiDecoyBB84LossyResults_1.00e+09_1decoy.mat");
matN10 = load("data/RenyiDecoyBB84LossyResults_1.00e+10_1decoy.mat");
matN11 = load("data/RenyiDecoyBB84LossyResults_1.00e+11_1decoy.mat");

%number of loss values used
numLoss = 26;

%key rates
keyRatesN5 = parseKeyRates(matN5,numLoss);
keyRatesN6 = parseKeyRates(matN6,numLoss);
keyRatesN7 = parseKeyRates(matN7,numLoss);
keyRatesN8 = parseKeyRates(matN8,numLoss);
keyRatesN9 = parseKeyRates(matN9,numLoss);
keyRatesN10 = parseKeyRates(matN10,numLoss);
keyRatesN11 = parseKeyRates(matN11,numLoss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results EUR
matEURN6 = readmatrix("EUR Data/ActiveDecoyBB84_1Decoy_1.00e+06.csv"); 
matEURN7 = readmatrix("EUR Data/ActiveDecoyBB84_1Decoy_1.00e+07.csv");
matEURN8 = readmatrix("EUR Data/ActiveDecoyBB84_1Decoy_1.00e+08.csv");
matEURN9 = readmatrix("EUR Data/ActiveDecoyBB84_1Decoy_1.00e+09.csv");
matEURN10 = readmatrix("EUR Data/ActiveDecoyBB84_1Decoy_1.00e+10.csv");
matEURN11 = readmatrix("EUR Data/ActiveDecoyBB84_1Decoy_1.00e+11.csv");

keyRatesEURN6 = matEURN6(:,2);
keyRatesEURN7 = matEURN7(:,2);
keyRatesEURN8 = matEURN8(:,2);
keyRatesEURN9 = matEURN9(:,2);
keyRatesEURN10 = matEURN10(:,2);
keyRatesEURN11 = matEURN11(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
% Loss etc. for plotting
lossdB = linspace(0,50,numLoss);
tempeta = 10.^(-lossdB/10);
% logAlpha = arrayfun(@(x) x.currentParams.logrenyiAlpha, matN11.results);

eta = tempeta(1:end);
etadB = -10*log10(eta).';

etaAsymp = 10.^(-linspace(0,max(etadB,[],"all"),100));
etaAsympdB = -10*log10(etaAsymp);

%List of total signals sent
Nlist = 10.^(5:1:11);

%Color list
colorList = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options for plots 
[x0,y0,width,height] = deal(50,100,600,500);


figure
set(gcf,'position',[x0,y0,width,height])

%GREAT
semilogy(etadB,keyRatesN5,"-o","Color",colorList(1),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(1))))

hold on
semilogy(etadB,keyRatesN6,"-o","Color",colorList(2),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(2))))
semilogy(etadB,keyRatesN7,"-o","Color",colorList(3),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(3))))
semilogy(etadB,keyRatesN8,"-o","Color",colorList(4),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(4))))
semilogy(etadB,keyRatesN9,"-o","Color",colorList(5),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(5))))
semilogy(etadB,keyRatesN10,"-o","Color",colorList(6),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(6))))
semilogy(etadB,keyRatesN11,"-o","Color",colorList(7),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(7))))

%EUR
semilogy(etadB,zeros(1,26),"--x","Color",colorList(1),"DisplayName", sprintf("N = 10^{%.0f} (EUR)",log10(Nlist(1))))
semilogy(etadB,keyRatesEURN6,"--x","Color",colorList(2),"DisplayName", sprintf("N = 10^{%.0f} (EUR)",log10(Nlist(2))))
semilogy(etadB,keyRatesEURN7,"--x","Color",colorList(3),"DisplayName", sprintf("N = 10^{%.0f} (EUR)",log10(Nlist(3))))
semilogy(etadB,keyRatesEURN8,"--x","Color",colorList(4),"DisplayName", sprintf("N = 10^{%.0f} (EUR)",log10(Nlist(4))))
semilogy(etadB,keyRatesEURN9,"--x","Color",colorList(5),"DisplayName", sprintf("N = 10^{%.0f} (EUR)",log10(Nlist(5))))
semilogy(etadB,keyRatesEURN10,"--x","Color",colorList(6),"DisplayName", sprintf("N = 10^{%.0f} (EUR)",log10(Nlist(6))))
semilogy(etadB,keyRatesEURN11,"--x","Color",colorList(7),"DisplayName", sprintf("N = 10^{%.0f} (EUR)",log10(Nlist(7))))

lgd = legend('NumColumns',2);
lgd.FontSize = 10;
lgd.Location = 'northeast';
xlabel('transmittance in dB',FontSize=14)
ylabel('Secret key rate',FontSize=14)
ylim([2*1e-6 1/exp(1)])
xlim([0 45])
hold off

%save figure
f1=gca;
filestr1 = "ActiveDecoyBB84_1Decoy.pdf";
% exportgraphics(f1,filestr1,'ContentType','vector')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rates = parseKeyRates(data,numElmts)
    results = data.results;
    list = [results(:).keyRate];
    rates = zeros(1,numElmts);
    numList = numel(list);
    rates(1:numList) = list;
end