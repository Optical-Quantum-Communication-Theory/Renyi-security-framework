classdef EqualityConstraintDecoy < BaseConstraint

    properties
        operator (:,:) double {mustBeNonempty,mustBeFinite,mustBeHermitian} = 0; % Hermitian operators inner producted with rho in the constraint.
        vector (1,:) double {mustBeNonNan,mustBeFinite,mustBeReal} = 0; % Real vector of values the inner product should be equal to.
    end

    properties (Dependent = true)
        rhoDim % Size of the input system's Hilbert space for this constraint.
        numDecoy
    end
    
    methods
        function obj = EqualityConstraintDecoy(operator,vector)
            arguments
                operator (:,:) double = 0;
                vector  (1,:) double = 0;
            end
            obj.operator = operator;
            obj.vector = vector;
        end

        %% getters and setters
        function operatorDim = get.rhoDim(obj)
            operatorDim = size(obj.operator,1);
        end

        function numDecoy = get.numDecoy(obj)
            numDecoy = numel(obj.vector);
        end
    end
end
