clear all
% get optimal values
optVarNames = {'logrenyiAlpha','probTest','GROUP_decoys_1'};

% number of points
numPoints = 20;

% results from data
mat = load("data/RenyiDecoyBB84LossyResults_1.00e+11_1decoy_17_end.mat");


% Parse key rates and optimal values
[keyRates, optvals] = parseKeyRatesAndOptVals(mat,numPoints,optVarNames);

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
    resultsCurrentParams = [results(:).currentParams];

    % extract logrenyiAlpha, probTest, signal intensity
    % listOptVals = [[resultsCurrentParams(:).logrenyiAlpha].'];
    listOptVals = [];
    for index = 1:numel(optValNames)
        elemnt = optValNames{index};
        listOptVals = [listOptVals,[resultsCurrentParams(:).(elemnt)].'];
    end

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