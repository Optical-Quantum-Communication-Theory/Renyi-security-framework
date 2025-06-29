%% Results GREAT
matN4 = load("data/RenyiBB84LossyResults_1.00e+04_depol.mat");
matN5 = load("data/RenyiBB84LossyResults_1.00e+05_depol.mat");
matN6 = load("data/RenyiBB84LossyResults_1.00e+06_depol.mat"); 
matN7 = load("data/RenyiBB84LossyResults_1.00e+07_depol.mat");
matN8 = load("data/RenyiBB84LossyResults_1.00e+08_depol.mat");
matN9 = load("data/RenyiBB84LossyResults_1.00e+09_depol.mat");
matN10 = load("data/RenyiBB84LossyResults_1.00e+10_depol.mat");

%number of loss values used
numLoss = 26;

%key rates
keyRatesN4 = parseKeyRates(matN4,numLoss);
keyRatesN5 = parseKeyRates(matN5,numLoss);
keyRatesN6 = parseKeyRates(matN6,numLoss);
keyRatesN7 = parseKeyRates(matN7,numLoss);
keyRatesN8 = parseKeyRates(matN8,numLoss);
keyRatesN9 = parseKeyRates(matN9,numLoss);
keyRatesN10 = parseKeyRates(matN10,numLoss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results EUR
matEURN5 = readmatrix("EUR Data/EURQubitBB84Results_1.00e+05CP_depol.csv");
matEURN6 = readmatrix("EUR Data/EURQubitBB84Results_1.00e+06CP_depol.csv"); 
matEURN7 = readmatrix("EUR Data/EURQubitBB84Results_1.00e+07CP_depol.csv");
matEURN8 = readmatrix("EUR Data/EURQubitBB84Results_1.00e+08CP_depol.csv");
matEURN9 = readmatrix("EUR Data/EURQubitBB84Results_1.00e+09CP_depol.csv");
matEURN10 = readmatrix("EUR Data/EURQubitBB84Results_1.00e+10CP_depol.csv");

keyRatesEURN5 = matEURN5(:,2);
keyRatesEURN6 = matEURN6(:,2);
keyRatesEURN7 = matEURN7(:,2);
keyRatesEURN8 = matEURN8(:,2);
keyRatesEURN9 = matEURN9(:,2);
keyRatesEURN10 = matEURN10(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
% Loss etc. for plotting
tempeta = arrayfun(@(x) x.currentParams.transmittance, matN10.results);
tempptest = arrayfun(@(x) x.currentParams.probTest, matN10.results);

eta = tempeta(1:end);
etadB = -10*log10(eta);

etaAsymp = 10.^(-linspace(0,max(etadB,[],"all"),100));
etaAsympdB = -10*log10(etaAsymp);

%List of total signals sent
Nlist = 10.^(4:1:10);

%Color list
colorList = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options for plots 
[x0,y0,width,height] = deal(50,100,600,500);


figure
set(gcf,'position',[x0,y0,width,height])

% semilogy(etadB,keyRatesN4,"-o","Color",colorList(1),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(1))))
semilogy(etadB,keyRatesN5,"-o","Color",colorList(2),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(2))))

hold on

semilogy(etadB,keyRatesN6,"-o","Color",colorList(3),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(3))))
semilogy(etadB,keyRatesN7,"-o","Color",colorList(4),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(4))))
semilogy(etadB,keyRatesN8,"-o","Color",colorList(5),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(5))))
semilogy(etadB,keyRatesN9,"-o","Color",colorList(6),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(6))))
semilogy(etadB,keyRatesN10,"-o","Color",colorList(7),"DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(7))))

semilogy(etadB,keyRatesEURN5,"--x","Color",colorList(2),"DisplayName", sprintf("N = 10^{%.0f} (EUR)",log10(Nlist(2))))
semilogy(etadB,keyRatesEURN6,"--x","Color",colorList(3),"DisplayName", sprintf("N = 10^{%.0f} (EUR)",log10(Nlist(3))))
semilogy(etadB,keyRatesEURN7,"--x","Color",colorList(4),"DisplayName", sprintf("N = 10^{%.0f} (EUR)",log10(Nlist(4))))
semilogy(etadB,keyRatesEURN8,"--x","Color",colorList(5),"DisplayName", sprintf("N = 10^{%.0f} (EUR)",log10(Nlist(5))))
semilogy(etadB,keyRatesEURN9,"--x","Color",colorList(6),"DisplayName", sprintf("N = 10^{%.0f} (EUR)",log10(Nlist(6))))
semilogy(etadB,keyRatesEURN10,"--x","Color",colorList(7),"DisplayName", sprintf("N = 10^{%.0f} (EUR)",log10(Nlist(7))))

lgd = legend('NumColumns',2);
lgd.FontSize = 10;
lgd.Location = 'northeast';
xlabel('transmittance in dB',FontSize=14)
ylabel('Secret key rate',FontSize=14)
ylim([1/2*1e-6 1])
hold off

%save figure
f1=gca;
filestr1 = "QubitBB84Depol.pdf";
% exportgraphics(f1,filestr1,'ContentType','vector')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rates = parseKeyRates(data,numElmts)
    results = data.results;
    list = [results(:).keyRate];
    rates = zeros(1,numElmts);
    numList = numel(list);
    rates(1:numList) = list;
end