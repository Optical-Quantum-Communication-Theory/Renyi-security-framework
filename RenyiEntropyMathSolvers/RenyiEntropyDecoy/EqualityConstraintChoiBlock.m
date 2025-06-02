classdef EqualityConstraintChoiBlock < BaseConstraint
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
        operatorCell (1,:) cell {mustBeNonempty,allCellsMustBeFinite(operatorCell),mustBeCellOf(operatorCell,"double"),allCellsMustBeHermitian} = {0}; % Hermitian operator inner producted with rho in the constraint.
        scalar (1,1) double {mustBeNonNan,mustBeFinite,mustBeReal} = 0; % Real scalar value the inner product should be equal to.
    end

    properties (Dependent = true)
        rhoDim % Size of the input systems for this constraint.
        numBlocks
    end
    
    methods
        function obj = EqualityConstraintChoiBlock(operatorCell,scalar)
            %EQUALITYCONSTRAINT Construct an instance of this class.
            % See class description above.
            %
            % Input:
            % * operatorCell (0): cell of nxn complex hermitian operator. It must have
            %   at least one element and finite in value.
            % * scalar (0): Finite, real valued scalar (also not nan).
            arguments
                operatorCell (:,:) cell = {0};
                scalar  (1,:) double = 0;
            end
            obj.operatorCell = operatorCell;
            obj.scalar = scalar;
        end

        %% getters and setters
        function operatorDims = get.rhoDim(obj)
            operatorDims = cellfun(@(x) size(x,1),obj.operatorCell);
        end

        function numBlocks = get.numBlocks(obj)
            numBlocks = numel(obj.operatorCell);
        end
    end
end

%Validation function
function allCellsMustBeFinite(operatorCell)
    allCellsMust(operatorCell,@(x) mustBeFinite(x))
end