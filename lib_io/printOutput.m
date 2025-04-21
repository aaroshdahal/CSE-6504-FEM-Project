function printOutput(Results)
    displacements = Results.displacements;    % nNode × nDof
    reactionForces = Results.reactions;        % (nNode*nDof) × 1
    stresses = Results.strains;        % 3 × nPoints × nElem
    strains = Results.stresses;       % 3 × nPoints × nElem

    [nNode, nDof] = size(displacements);
    [~, nPts, nElem] = size(strains);

    %----------------------------------------------------------------------
    fprintf('\n===================================================================\n');
    fprintf('                     O U T P U T   S U M M A R Y\n');
    fprintf('===================================================================\n\n');

    % 1) Nodal Displacements
    fprintf('------------------------\n');
    fprintf('  Nodal Displacements\n');
    fprintf('------------------------\n');
    for i = 1:nNode
        if nDof==2
            fprintf('[%2d]  %10.5f   %10.5f\n', i, displacements(i,1), displacements(i,2));
        else
            fprintf('[%2d]  ', i);
            fprintf('%10.5f  ', displacements(i,:));
            fprintf('\n');
        end
    end
    fprintf('\n');

    % 2) Reaction Forces (only non-zero entries)
    fprintf('------------------------\n');
    fprintf('   Reaction Forces\n');
    fprintf('------------------------\n');
    for dof = 1:length(reactionForces)
        if abs(reactionForces(dof))>0  % prints only DOFs with reactions
            fprintf('[%2d]  %10.5f\n', dof, reactionForces(dof));
        end
    end
    fprintf('\n');

    % 3) Strains & Stresses at Integration Points
    fprintf('-------------------------------------------------------------\n');
    fprintf('   Strains and Stresses at Gauss Points (εx, εy, γxy → σx, σy, τxy)\n');
    fprintf('-------------------------------------------------------------\n\n');

    for e = 1:nElem
        fprintf(' Element Number = %d\n', e);
        for ip = 1:nPts
            eps = stresses(:,ip,e);
            sig = strains(:,ip,e);
            fprintf('   Pt %d:  ε = [ %8.5f  %8.5f  %8.5f ]   σ = [ %8.5f  %8.5f  %8.5f ]\n', ...
                    ip, eps(1), eps(2), eps(3), sig(1), sig(2), sig(3));
        end
        fprintf('\n');
    end
end

