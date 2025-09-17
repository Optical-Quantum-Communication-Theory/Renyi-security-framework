function [minProbVec,maxProbVec] = minmaxPoisVec(mu,deltaMu,photonNumVec)
    arguments
        mu (1,1) double {mustBeNonnegative}
        deltaMu (1,1) double {mustBeNonnegative}
        photonNumVec (:,1) double {mustBeInteger}
    end
    %Number of elements 
    numPh = numel(photonNumVec);

    %Initialize vectors containing mins and maxs
    minProbVec = zeros(numPh,1);
    maxProbVec = zeros(numPh,1);

    %iterate over photon number
    for index = 1:numPh
        [minProbVec(index),maxProbVec(index)] = minmaxPois(mu,deltaMu,photonNumVec(index));
    end
end

function [minProb,maxProb] = minmaxPois(mu,deltaMu,photonNum)
    arguments
        mu (1,1) double {mustBeNonnegative}
        deltaMu (1,1) double {mustBeNonnegative}
        photonNum (1,1) double {mustBeInteger}
    end

    if photonNum == 0
        %For 0 photons min and max of Poisson dist on interval [mu,mu+delta]
        % are always the end points
        minProb = Poisson(mu+deltaMu,0);
        maxProb = Poisson(mu,0);
    else
        %For photonNum =/= 0 min of Poisson dist on interval [mu,mu+delta]
        % is achieved at either end points
        minProb = min(Poisson(mu,photonNum), ...
                      Poisson(mu + deltaMu,photonNum));
        
        %For photonNum =/= 0 max of Poisson dist on interval [mu,mu+delta]
        % is achieved at either end points or at mu = photonNum if
        % photonNum in [mu,mu+delta]

        %Check if photonNum in [mu,mu+delta]
        insideRange = discretize(photonNum,[mu,mu + deltaMu]) == 1;

        if insideRange
            %if photonNum in [mu,mu+delta] 3 possible values for mu
            maxProb = max([Poisson(mu,photonNum), ...
                      Poisson(mu + deltaMu,photonNum), ...
                      Poisson(photonNum,photonNum)],[],'all');
        else
            %if photonNum not in [mu,mu+delta] only end points need
            %to be considered
            maxProb = max([Poisson(mu,photonNum), ...
                      Poisson(mu + deltaMu,photonNum)],[],'all');
        end
    end
end

function prob = Poisson(mu,m)
    prob = mu.^m.*exp(-mu)./factorial(m);
end