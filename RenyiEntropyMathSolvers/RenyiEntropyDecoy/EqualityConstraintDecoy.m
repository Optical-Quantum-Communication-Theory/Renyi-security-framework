classdef EqualityConstraintDecoy < BaseConstraint
    %EQUALITYCONSTRAINT A simple class to model equality constraints for
    %   the Hilbert space of hermitian operators.
    %   Constraints are of the form tr[operator'*rho] = scalar, up to some
    %   error tolerance. This is an affine constraint.
    %
    % Properties:
    % * operator: n x n complex hermitian operator. It must have at least
    %   one element and be finite in value.
    % * vector: Finite, real valued vector (also not nan).
    % * rhoDim: Size of the input system's Hilbert space for this
    %   constraint. (n x n = rhoDim x rhoDim)
    %
    % See also: BaseConstraint
    properties
        operator (:,:) double {mustBeNonempty,mustBeFinite,mustBeHermitian} = 0; % Hermitian operator inner producted with rho in the constraint.
        vector (1,:) double {mustBeNonNan,mustBeFinite,mustBeReal} = 0; % Real vector of values the inner product should be equal to.
    end

    properties (Dependent = true)
        rhoDim % Size of the input system's Hilbert space for this constraint.
        numDecoy
    end
    
    methods
        function obj = EqualityConstraintDecoy(operator,vector)
            %EQUALITYCONSTRAINT Construct an instance of this class.
            % See class description above.
            %
            % Input:
            % * operator (0): nxn complex hermitian operator. It must have
            %   at least one element and finite in value.
            % * vector (0): Finite, real valued vector (also not nan).
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
