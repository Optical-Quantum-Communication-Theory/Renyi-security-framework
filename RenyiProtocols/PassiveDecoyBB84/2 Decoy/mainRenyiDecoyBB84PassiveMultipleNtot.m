%pick preset
qkdInput = RenyiDecoyBB84PassivePreset();

%List of mutiple total signals sent
N_list = [1e6,1e8,1e10];

%Loss
%total array of loss values to iterate over
lossdB = linspace(0,30,16);
transmittance = 10.^(-lossdB/10);

%list of maximal element of loss
lossList = [8,12,16];

for index = 1:numel(N_list)
    %Add total signals sent from list above
    qkdInput.addFixedParameter("Ntot",N_list(index))

    %Add loss until element from list above
    qkdInput.addScanParameter("transmittance", num2cell(transmittance(1:lossList(index))));

    %run the QKDSolver with this input and store results
    results = MainIteration(qkdInput);

    filestr = sprintf("RenyiDecoyBB84PassiveResults_%.2e",N_list(index)) ...
        + sprintf("_%.2e",0.1) + ".mat";

    % save the results and preset to a file.;
    % save(filestr,"results","qkdInput");
end

%% plot the result
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")