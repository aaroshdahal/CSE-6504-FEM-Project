function [displacements, reactionForces] = solution(K, R, bcIndex, bcValues)
% Imposes essential BCs, solves K u = R, and computes reactions

dofTotal  = length(R);
dofActive  = setdiff(1:dofTotal, bcIndex);

% partition
Kff = K(dofActive, dofActive);
Rf  = R(dofActive) - K(dofActive, bcIndex)*bcValues;

% solve for free disp
u_free = Kff \ Rf;

% assemble global displacement vector
displacements        = zeros(dofTotal,1);
displacements(dofActive)   = u_free;
displacements(bcIndex)    = bcValues;

% reaction forces at fixed dofs:  
reactionForces = zeros(dofTotal,1);
reactionForces(bcIndex) = K(bcIndex,:) * displacements - R(bcIndex);
end


