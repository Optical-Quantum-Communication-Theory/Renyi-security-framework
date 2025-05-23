function allCellsMust(cellArray,func)
% allCellsMust Validation function that applies a validation function to every
% member of a cell array.
%
% Inputs:
% * cellArray: The cell array to apply the validation function on
% * func: the function to apply to each member of the cell array
try
    cellfun(func,cellArray);
catch err
    err.throwAsCaller();
end
end