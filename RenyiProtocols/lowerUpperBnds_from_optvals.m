%% helper function
function bounds = lowerUpperBnds_from_optvals(index, values, default_lower, default_upper)
% f returns [lower_bound, upper_bound] using neighbours and a signed 10% margin
%
% If the selected bound is negative, the 10% margin is flipped:
%   - Lower bound: add 10% if negative, subtract 10% if positive
%   - Upper bound: subtract 10% if negative, add 10% if positive
%
% No margin is applied if the value is a default.

    n = length(values);

    if index < 1 || index > n
        error('Index out of bounds.');
    end

    % Get values and whether they are default
    if index == 1
        val_prev = default_lower;
        is_prev_default = true;
    else
        val_prev = values(index - 1);
        is_prev_default = false;
    end

    val_curr = values(index);
    is_curr_default = false;

    if index == n
        val_next = default_upper;
        is_next_default = true;
    else
        val_next = values(index + 1);
        is_next_default = false;
    end

    vals = [val_prev, val_curr, val_next];
    is_default = [is_prev_default, is_curr_default, is_next_default];

    % Lower bound
    [lower_bound, idx_min] = min(vals);
    if ~is_default(idx_min)
        if lower_bound < 0
            lower_bound = lower_bound * 1.1; % Add 10% if negative
        else
            lower_bound = lower_bound * 0.9; % Subtract 10% if positive
        end
    end

    % Upper bound
    [upper_bound, idx_max] = max(vals);
    if ~is_default(idx_max)
        if upper_bound < 0
            upper_bound = upper_bound * 0.9; % Subtract 10% if negative
        else
            upper_bound = upper_bound * 1.1; % Add 10% if positive
        end
    end

    bounds = [lower_bound, upper_bound];
end
