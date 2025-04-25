function printOutput(Results)
    displacements = Results.displacements;    % nNode × nDof
    reactionForces = Results.reactions;        % (nNode*nDof) × 1
    stresses = Results.strains;        % 3 × nPoints × nElem
    strains = Results.stresses;       % 3 × nPoints × nElem

    coords    = Results.coordinates;          % nNode×nDim
    elems     = Results.elements;             % nElem×nElemNode
    nDim      = Results.nDim;
    nElemNode = Results.nElemNode;
    [xiList, ~] = integrationPoints(nDim, nElemNode);

    [nNode, nDof] = size(displacements);
    [~, nPts, nElem] = size(strains);

    %----------------------------------------------------------------------
    fprintf('\n\n================================================================================\n');
    fprintf('                             O U T P U T   S U M M A R Y\n');
    fprintf('================================================================================\n');

    % 1) Nodal Displacements
    fprintf('%s\n', repmat('_',1,80));
    fprintf('                                 Nodal Displacements\n');
    fprintf('--------------------------------------------------------------------------------\n');
    for i = 1:nNode
        if nDof==2
            fprintf('[%2d]     %10.5f      %10.5f\n', i, displacements(i,1), displacements(i,2));
        else
            fprintf('[%2d]  ', i);
            fprintf('%10.5f  ', displacements(i,:));
            fprintf('\n');
        end
    end
    fprintf('\n');

    % 2) Reaction Forces (only non-zero entries)
    fprintf('%s\n', repmat('_',1,80));
    fprintf('                                    Reaction Forces\n');
    fprintf('--------------------------------------------------------------------------------\n');
    counter = 1;
    for dof = 1:length(reactionForces)
        if abs(reactionForces(dof))>0  % prints only DOFs with reactions
            fprintf('[%2d]     %10.5f\n', counter, reactionForces(dof));
            counter = counter + 1;
        end
    end
    fprintf('\n');

    % 3) Strains & Stresses at Integration Points
    fprintf('%s\n', repmat('_',1,124));
    fprintf('                                         Strains and Stresses at Gauss Points\n');
    fprintf('----------------------------------------------------------------------------------------------------------------------------\n');

    for e = 1:nElem

    fprintf('                                                  Element Number = %d\n', e);
    fprintf('----------------------------------------------------------------------------------------------------------------------------\n');
    % column titles
    fprintf(' Integration        Location             Location                         Strains                        Stresses\n');
    fprintf('    Point          (Natural)             (Global)               e11        e22       e12         s11       s22       s12\n');
    fprintf('----------------------------------------------------------------------------------------------------------------------------\n');

        for ip = 1:nPts
            xi = xiList(:,ip);
            [N, ~] = shapeFunctions(nDim,nElemNode, xi); 
            conn = elems(e,:);
            xy   = N' * coords(conn,:);

            eps = stresses(:,ip,e);
            sig = strains(:,ip,e);
            fprintf('  %d           (%8.4f, %8.4f)   (%8.4f, %8.4f)    %8.4f  %8.4f   %8.4f     %8.4f  %8.4f  %8.4f \n', ...
                    ip, xi(1), xi(1), xy(1), xy(2), eps(1), eps(2), eps(3), sig(1), sig(2), sig(3));
        end
        fprintf('\n');
    end
    fprintf('%s\n', repmat('=',1,124));
end

