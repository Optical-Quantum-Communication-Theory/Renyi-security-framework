%% Results GREAT
matN5 = load("data/RenyiDecoyBB84LossyResults_1.00e+05.mat");
matN6 = load("data/RenyiDecoyBB84LossyResults_1.00e+06.mat"); 
matN7 = load("data/RenyiDecoyBB84LossyResults_1.00e+07.mat");
matN8 = load("data/RenyiDecoyBB84LossyResults_1.00e+08.mat");
matN9 = load("data/RenyiDecoyBB84LossyResults_1.00e+09.mat");
matN10 = load("data/RenyiDecoyBB84LossyResults_1.00e+10.mat");
matN11 = load("data/RenyiDecoyBB84LossyResults_1.00e+11.mat");

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

%Optimal Intensities
optMuN5 = parseOptInt(matN5,numLoss);
optMuN6 = parseOptInt(matN6,numLoss);
optMuN7 = parseOptInt(matN7,numLoss);
optMuN8 = parseOptInt(matN8,numLoss);
optMuN9 = parseOptInt(matN9,numLoss);
optMuN10 = parseOptInt(matN10,numLoss);
optMuN11 = parseOptInt(matN11,numLoss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results EUR
matEURN6 = readmatrix("EUR Data/EURDecoyBB84Results_1.00e+06CP.csv"); 
matEURN7 = readmatrix("EUR Data/EURDecoyBB84Results_1.00e+07CP.csv");
matEURN8 = readmatrix("EUR Data/EURDecoyBB84Results_1.00e+08CP.csv");
matEURN9 = readmatrix("EUR Data/EURDecoyBB84Results_1.00e+09CP.csv");
matEURN10 = readmatrix("EUR Data/EURDecoyBB84Results_1.00e+10CP.csv");
matEURN11 = readmatrix("EUR Data/EURDecoyBB84Results_1.00e+11CP.csv");

keyRatesEURN6 = matEURN6(:,2);
keyRatesEURN7 = matEURN7(:,2);
keyRatesEURN8 = matEURN8(:,2);
keyRatesEURN9 = matEURN9(:,2);
keyRatesEURN10 = matEURN10(:,2);
keyRatesEURN11 = matEURN11(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results Postselection
matPSN6 = readmatrix("PS Data/DecoyBB84Adaptive_PS_1.00e+06.csv"); 
matPSN7 = readmatrix("PS Data/DecoyBB84Adaptive_PS_1.00e+07.csv");
matPSN8 = readmatrix("PS Data/DecoyBB84Adaptive_PS_1.00e+08.csv");
matPSN9 = readmatrix("PS Data/DecoyBB84Adaptive_PS_1.00e+09.csv");
matPSN10 = readmatrix("PS Data/DecoyBB84Adaptive_PS_1.00e+10.csv");
matPSN11 = readmatrix("PS Data/DecoyBB84Adaptive_PS_1.00e+11.csv");

%key rates
keyRatesPSN6 = matPSN6(:,2);
keyRatesPSN7 = matPSN7(:,2);
keyRatesPSN8 = matPSN8(:,2);
keyRatesPSN9 = matPSN9(:,2);
keyRatesPSN10 = matPSN10(:,2);
keyRatesPSN11 =  matPSN11(:,2);

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
rRenyi(1) = semilogy(etadB,keyRatesN11,"-o","Color",colorList(7),"DisplayName", sprintf("n = 10^{%.0f}",log10(Nlist(7))));

hold on
rRenyi(2) = semilogy(etadB,keyRatesN10,"-o","Color",colorList(6),"DisplayName", sprintf("n = 10^{%.0f}",log10(Nlist(6))));
rRenyi(3) = semilogy(etadB,keyRatesN9,"-o","Color",colorList(5),"DisplayName", sprintf("n = 10^{%.0f}",log10(Nlist(5))));
rRenyi(4) = semilogy(etadB,keyRatesN8,"-o","Color",colorList(4),"DisplayName", sprintf("n = 10^{%.0f}",log10(Nlist(4))));
rRenyi(5) = semilogy(etadB,keyRatesN7,"-o","Color",colorList(3),"DisplayName", sprintf("n = 10^{%.0f}",log10(Nlist(3))));
rRenyi(6) = semilogy(etadB,keyRatesN6,"-o","Color",colorList(2),"DisplayName", sprintf("n = 10^{%.0f}",log10(Nlist(2))));
rRenyi(7) = semilogy(etadB,keyRatesN5,"-o","Color",colorList(1),"DisplayName", sprintf("n = 10^{%.0f}",log10(Nlist(1))));

%EUR (All 0 key rates commented out)
% semilogy(etadB,zeros(1,26),"--x","Color",colorList(1),"DisplayName", sprintf("n = 10^{%.0f} (EUR)",log10(Nlist(1))))
% semilogy(etadB,keyRatesEURN6,"--x","Color",colorList(2),"DisplayName", sprintf("n = 10^{%.0f} (EUR)",log10(Nlist(2))))
% semilogy(etadB,keyRatesEURN7,"--x","Color",colorList(3),"DisplayName", sprintf("n = 10^{%.0f} (EUR)",log10(Nlist(3))))
rEUR(1) = semilogy(etadB,keyRatesEURN11,"--x","Color",colorList(7),"DisplayName", sprintf("n = 10^{%.0f} (EUR)",log10(Nlist(7))));
rEUR(2) = semilogy(etadB,keyRatesEURN10,"--x","Color",colorList(6),"DisplayName", sprintf("n = 10^{%.0f} (EUR)",log10(Nlist(6))));
rEUR(3) = semilogy(etadB,keyRatesEURN9,"--x","Color",colorList(5),"DisplayName", sprintf("n = 10^{%.0f} (EUR)",log10(Nlist(5))));
rEUR(4) = semilogy(etadB,keyRatesEURN8,"--x","Color",colorList(4),"DisplayName", sprintf("n = 10^{%.0f} (EUR)",log10(Nlist(4))));
rEUR(5) = plot(nan,nan,'LineStyle', 'none',DisplayName='');
rEUR(6) = plot(nan,nan,'LineStyle', 'none',DisplayName='');
rEUR(7) = plot(nan,nan,'LineStyle', 'none',DisplayName='');

% %PS (All 0 key rates commented out)
% rPS(1) = semilogy(etadB,keyRatesPSN11,"-.^","Color",colorList(7),"DisplayName", sprintf("n = 10^{%.0f} (PS)",log10(Nlist(7))));
% rPS(2) = semilogy(etadB,keyRatesPSN10,"-.^","Color",colorList(6),"DisplayName", sprintf("n = 10^{%.0f} (PS)",log10(Nlist(6))));
% rPS(3) = semilogy(etadB,keyRatesPSN9,"-.^","Color",colorList(5),"DisplayName", sprintf("n = 10^{%.0f} (PS)",log10(Nlist(5))));
% rPS(4) = plot(nan,nan,'LineStyle', 'none',DisplayName='');
% rPS(5) = plot(nan,nan,'LineStyle', 'none',DisplayName='');
% rPS(6) = plot(nan,nan,'LineStyle', 'none',DisplayName='');
% rPS(7) = plot(nan,nan,'LineStyle', 'none',DisplayName='');
% % semilogy(etadB,zeros(1,26),"-.^","Color",colorList(1),"DisplayName", sprintf("n = 10^{%.0f} (PS)",log10(Nlist(1))))
% % semilogy(etadB,zeros(1,26),"-.^","Color",colorList(2),"DisplayName", sprintf("n = 10^{%.0f} (PS)",log10(Nlist(2))))
% % semilogy(etadB,zeros(1,26),"-.^","Color",colorList(3),"DisplayName", sprintf("n = 10^{%.0f} (PS)",log10(Nlist(3))))
% % semilogy(etadB,zeros(1,26),"-.^","Color",colorList(4),"DisplayName", sprintf("n = 10^{%.0f} (PS)",log10(Nlist(4))))


lgd = legend('NumColumns',2);
lgd.FontSize = 10;
lgd.Location = 'northeast';
xlabel('Loss in dB',FontSize=14)
ylabel('Secret key rate',FontSize=14)

ylim([2*1e-6 1/exp(1)])
xlim([0 45])
hold off

%save figure
f1=gca;
filestr1 = "ActiveDecoyBB84_2Decoy";
% exportgraphics(f1,filestr1 + ".pdf",'ContentType','vector')
% exportgraphics(f1,filestr1 + ".eps",'ContentType','vector')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define distinct colors for each method
methodColors = lines(3); 
cRenyi = methodColors(1,:); % Color 1 for Standard/Renyi
cEUR   = methodColors(2,:); % Color 2 for EUR
cPS    = methodColors(3,:); % Color 3 for PS

figure
set(gcf,'position',[x0,y0,width,height])

% --- Renyi rates (Method 1: Solid Line, Circle Marker) ---
semilogy(etadB,keyRatesN6,"-^","Color",cRenyi,"DisplayName", sprintf("n = 10^{%.0f}",log10(Nlist(2))))
hold on
semilogy(etadB,keyRatesN11,"-o","Color",cRenyi,"DisplayName", sprintf("n = 10^{%.0f}",log10(Nlist(7))))


% --- EUR rates (Method 2: Dashed Line, Cross Marker) ---
semilogy(etadB,keyRatesEURN8,"--^","Color",cEUR,"DisplayName", sprintf("n = 10^{%.0f} (EUR)",log10(Nlist(4))))
semilogy(etadB,keyRatesEURN11,"--o","Color",cEUR,"DisplayName", sprintf("n = 10^{%.0f} (EUR)",log10(Nlist(7))))


% --- PS rates (Method 3: Dot-Dash Line, Triangle Marker) ---
semilogy(etadB,keyRatesPSN9,"-.^","Color",cPS,"DisplayName", sprintf("n = 10^{%.0f} (PS)",log10(Nlist(5))))
semilogy(etadB,keyRatesPSN11,"-.o","Color",cPS,"DisplayName", sprintf("n = 10^{%.0f} (PS)",log10(Nlist(7))))

% Formatting
lgd = legend('NumColumns',3);
lgd.FontSize = 10;
lgd.Location = 'northeast';
xlabel('Loss in dB',FontSize=14)
ylabel('Secret Key Rate',FontSize=14)
ylim([1/2*1e-6 1])
hold off

% Save figure
f1=gca;
filestr1 = "ActiveDecoyBB84_2Decoy_Comparison.pdf";
% exportgraphics(f1,filestr1,'ContentType','vector')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot optimal intensity

figure
set(gcf,'position',[x0,y0,width,height])
plot(etadB(1),optMuN5(1),"--o","DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(1))), "Color",colorList(1))
hold on

plot(etadB(1:3),optMuN6(1:3),"--o","DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(2))), "Color",colorList(2))
plot(etadB(1:6),optMuN7(1:6),"--o","DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(3))), "Color",colorList(3))
plot(etadB(1:12),optMuN8(1:12),"--o","DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(4))), "Color",colorList(4))
plot(etadB(1:17),optMuN9(1:17),"--o","DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(5))), "Color",colorList(5))
plot(etadB(1:20),optMuN10(1:20),"--o","DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(6))), "Color",colorList(6))
plot(etadB(1:22),optMuN11(1:22),"--o","DisplayName", sprintf("N = 10^{%.0f}",log10(Nlist(7))), "Color",colorList(7))


lgd = legend(NumColumns=2);
lgd.FontSize = 10;
lgd.Location = 'best';
xlabel('Loss in dB',FontSize=14)
ylabel('Optimal Signal Intensity',FontSize=14)
% ylim([0.5 1])
hold off

f2=gca;
filestr2 = "Renyi_DecoyBB84_optimal_mu.pdf";
exportgraphics(f2,filestr2,'ContentType','vector')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rates = parseKeyRates(data,numElmts)
    results = data.results;
    list = [results(:).keyRate];
    rates = zeros(1,numElmts);
    numList = numel(list);
    rates(1:numList) = list;
end

function optInt = parseOptInt(data,numElmts)
    results = data.results;
    list = arrayfun(@(x) x.currentParams.GROUP_decoys_1, results);
    optInt = zeros(1,numElmts);
    numList = numel(list);
    optInt(1:numList) = list;
end