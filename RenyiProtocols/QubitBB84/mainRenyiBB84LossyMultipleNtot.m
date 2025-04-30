%pick preset
qkdInput = RenyiBB84LossyPreset();

%List of mutiple total signals sent
% N_list = [1e4,1e5,1e6,1e7,1e8,1e9,1e10];
N_list = [1e5,1e6,1e7,1e8,1e9,1e10];
% N_list = [1e10];
% N_list = [1e9,1e10,1e11];

%Loss
%total array of loss values to iterate over
lossdB = linspace(0,50,26);
transmittance = 10.^(-lossdB/10);

%list of maximal element of loss
% lossList = [2,8,12,16,20,24,26];
lossList = [8,12,16,20,24,26];

% lossList = [1];

for index = 1:numel(N_list)
    %Add total signals sent from list above
    qkdInput.addFixedParameter("Ntot",N_list(index))

    %Add loss until element from list above
    qkdInput.addScanParameter("transmittance", num2cell(transmittance(1:lossList(index))));

    % run the QKDSolver with this input
    results = MainIteration(qkdInput);
    filestr = sprintf("RenyiBB84LossyResults_%.2e",N_list(index)) + "_depol.mat";

    % save the results and preset to a file.;
    save(filestr,"results","qkdInput");
end

%% plot the result
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")
