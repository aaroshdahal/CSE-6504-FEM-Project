%%                           MATLAB PROJECT                   
%__________________________________________________________________________
% 
%                    CEE-6504 Finite Element Methods
%                            Dr. Rafi Muhanna
%                   Georgia Institute of Technology
%                              Spring 2024                      
%__________________________________________________________________________
% 
%           
%__________________________________________________________________________

function Results = fem2D(Model)
% The main FEM function to analyze a 2D plane stress/strain structure
% defined by the structure array "Model"

%% (A) Extract parameters from Model (do NOT modify)
problem     = Model.analysisType;
nDim        = Model.nDim;
nDof        = Model.nDof;
nNode       = Model.nNode;
nElem       = Model.nElem;
nElemNode   = Model.nElemNode;
coordinates = Model.geometry.coordinates;
elements    = Model.geometry.elements;
materials   = Model.sections(Model.elemSectId,:);

%% (B) Assemble global stiffness matrix K
K = zeros(nDof * nNode);
for e = 1:nElem
    conn      = elements(e,:);
    elemCoord = coordinates(conn,:);
    k_e       = elemStiff(elemCoord, materials(e,:), problem);
    % assemble into K
    idx = zeros(1, nElemNode*nDof);
    for i = 1:nElemNode
        base        = nDim*(conn(i)-1);
        idx(2*i-1:2*i) = base + (1:nDim);
    end
    K(idx,idx) = K(idx,idx) + k_e;
end

%% (C) Assemble global force vector F

F = zeros(nDof * nNode,1);

% 1) Body force b over the entire domain
%b = [0; -11000];  % body force vector [bx; by]
b = [0; 0];
[xiList, wList] = integrationPoints(nDim, nElemNode);
nGauss = length(wList);
for e = 1:nElem
    conn      = elements(e,:);
    elemCoord = coordinates(conn,:);
    f_loc     = zeros(nDim*nElemNode,1);
    for ip = 1:nGauss
        xi = xiList(:,ip);
        w  = wList(ip);
        [N, dNdxi] = shapeFunctions(nDim, nElemNode, xi);
        J    = dNdxi' * elemCoord;
        detJ = det(J);
        for iNode = 1:nElemNode
            idx_loc = nDim*(iNode-1)+(1:nDim);
            f_loc(idx_loc) = f_loc(idx_loc) + N(iNode)*b*detJ*w;
        end
    end
    % assemble local f into global F
    idx_glob = zeros(1, nDim*nElemNode);
    for iNode = 1:nElemNode
        base         = nDim*(conn(iNode)-1);
        idx_glob(2*iNode-1:2*iNode) = base + (1:nDim);
    end
    F(idx_glob) = F(idx_glob) + f_loc;
end

% 2) Hard‐coded traction on any element edge lying on y = 2
tol  = 1e-6;       % tolerance for y==2 check
%tval = 0
tval = 1e8;      % traction magnitude (global y direction)
% determine number of edges per element
if any(nElemNode == [3,6])
    nEdges = 3;
elseif any(nElemNode == [4,8])
    nEdges = 4;
else
    nEdges = 2;  % bar/truss
end

for e = 1:nElem
    conn = elements(e,:);
    for edgeID = 1:nEdges
        loc   = edgeNodes(nDim, nElemNode, edgeID);
        gIDs  = conn(loc);
        pts   = coordinates(gIDs,:);
        if all(abs(pts(:,2) - 20) < tol)
            fe = elemEquivalentLoad(pts, [0, tval]);
            % assemble fe into F
            idx_glob = zeros(1, nDim*length(loc));
            for j = 1:length(loc)
                idx_glob(2*j-1:2*j) = nDim*(gIDs(j)-1) + (1:nDim);
            end
            F(idx_glob) = F(idx_glob) + fe;
        end
    end
end

% 3) Optional concentrated nodal loads from input file (if any)
if Model.nNodalForce > 0
    for iN = 1:Model.nNodalForce
        nd = Model.loading.nodalForces(iN,1);
        fx = Model.loading.nodalForces(iN,2);
        fy = Model.loading.nodalForces(iN,3);
        F(nDim*(nd-1)+1) = F(nDim*(nd-1)+1) + fx;
        F(nDim*(nd-1)+2) = F(nDim*(nd-1)+2) + fy;
    end
end

%% (D) Impose hard‐coded Dirichlet BCs at nodes with y = 2 and solve
tol      = 1e-6;
bcIndex  = [];
bcValues = [];

for n = 1:nNode
    x = coordinates(n,1);
    y = coordinates(n,2);

    % --- Example: fix Ux = 0 on the left edge x=0
    % if (abs(x - 0) < tol) && (abs(y - 0) < tol)
    %     dof_x = nDim*(n-1) + 1;    % the x-DOF index of node n
    %     bcIndex  = [bcIndex,  dof_x];
    %     bcValues = [bcValues, 0.0];
    % end
    % 
    % % --- Example: prescribe Uy = 0.1 on the top edge y=2
    % if abs(y - 0) < tol
    %     dof_y = nDim*(n-1) + 2;    % the y-DOF index of node n
    %     bcIndex  = [bcIndex,  dof_y];
    %     bcValues = [bcValues, 0.0];
    % end

    if abs(x - 0) < tol
        dof_x = nDim*(n-1) + 1;    % the x-DOF index of node n
        bcIndex  = [bcIndex,  dof_x];
        bcValues = [bcValues, 0.0];
    end

    % --- Example: prescribe Uy = 0.1 on the top edge y=2
    if abs(y - 0) < tol
        dof_y = nDim*(n-1) + 2;    % the y-DOF index of node n
        bcIndex  = [bcIndex,  dof_y];
        bcValues = [bcValues, 0.0];
    end
end

% make sure it’s a column
[disp_vec, reactionForces] = solution(K, F, bcIndex, bcValues(:));
displacements = reshape(disp_vec, nDim, nNode)';

%% (E) Compute strains and stresses at integration points
[xiList, ~] = integrationPoints(nDim, nElemNode);
nPts    = size(xiList,2);
stresses = zeros(3, nPts, nElem);
strains  = zeros(3, nPts, nElem);

for e = 1:nElem
    conn      = elements(e,:);
    elemCoord = coordinates(conn,:);
    elemDisp  = displacements(conn,:);
    E  = materials(e,1);
    nu = materials(e,2);
    D  = strainStressMatrix(E, nu, problem);
    for ip = 1:nPts
        xi = xiList(:,ip);
        [~, dNdxi] = shapeFunctions(nDim, nElemNode, xi);
        J    = dNdxi' * elemCoord;
        dNdx = J \ dNdxi';
        % build B matrix
        B = zeros(3, nDim*nElemNode);
        for iNode = 1:nElemNode
            Bi = [ dNdx(1,iNode),       0;
                   0,       dNdx(2,iNode);
                   dNdx(2,iNode), dNdx(1,iNode) ];
            B(:,2*iNode-1:2*iNode) = Bi;
        end
        ue      = reshape(elemDisp', [], 1);
        epsilon = B * ue;      % [εx; εy; γxy]
        sigma   = D * epsilon;
        strains(:,ip,e)       = epsilon;
        stresses(:,ip,e)      = sigma;
    end
end

%% (F) Store results into output structure
Results.KGlobal       = K;
Results.FGlobal       = F;
Results.displacements = displacements;
Results.reactions     = reactionForces;
Results.strains       = strains;
Results.stresses      = stresses;

% for any downstream post-processing
Results.coordinates = coordinates;
Results.elements    = elements;
Results.nDim        = nDim;
Results.nElemNode   = nElemNode;
Results.materials   = materials;
Results.problem     = problem;

end
