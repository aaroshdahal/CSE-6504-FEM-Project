%%                           MATLAB PROJECT                  
%__________________________________________________________________________
% 
%                    CEE-6504 Finite Element Methods
%                            Dr. Rafi Muhanna
%                   Georgia Institute of Technology
%                              Spring 2024                      
%__________________________________________________________________________
% 
%                    Developed by  Shahrokh Shahi 
%__________________________________________________________________________

function Results = fem2D(Model)
% The main FEM function to analyze a 2d plane stress/strain structure
% defined by structure array "Model"


%% (A) Parameters (You do NOT need to modify this section):
%{
In the following, all the required parameters will be extracted from the
Model structure (all these data have already been read from the input
file):

(1) problem    : specifies the problem type ('plane stress' or 'plane strain' 
(2) nDim       : coordinates dimension (in this problem it is 2, since it is 2D problem)
(3) nDof       : number of degrees of freedom per each node (in this problem it is 2: u,v)
(4) nNode      : total number of nodes
(5) nElem      : total number of elements
(6) nElemNode  : number of nodes per each element e.g. 4 for a first-order quadrilateral element
(7) coordinates: a matrix including all nodes coordinates            (nNode x nDim     )
(8) elements   : a matrix including all elements connectivity data   (nElem x nElemNode)
(9) materials  : a matrix including all element material properties  (nElem x 3        )   (3: E, nu, t)
%}

problem     = Model.analysisType;

nDim        = Model.nDim ;
nDof        = Model.nDof ;
nNode       = Model.nNode;
nElem       = Model.nElem;
nElemNode   = Model.nElemNode;

coordinates = Model.geometry.coordinates;
elements    = Model.geometry.elements;
materials   = Model.sections(Model.elemSectId,:);

%% (B) Forming Global Stiffness Matrix

% initializing the variable to store global stiffness matrix
K = zeros(nDof * nNode);

% a loop over all elements to calculate elements' stiffness matrices and
% assemble the global stiffness matrix
for iElem = 1 : nElem
    connections = elements(iElem,:);           % a vector of element node numbers 
    elemCoord   = coordinates(connections,:);  % a matrix of element nodal coordinates    
    
    % call "elemStiff" function to calculate the element stiffness matrix
    k  = elemStiff(elemCoord, materials(iElem,:), problem);

    % assembly: finding the indices of corresponding entries in global 
    % matrix and assemble the global stiffness matrix

    % global‐index vector
    index = zeros(1, nElemNode*nDof);
    for i = 1:nElemNode
        base = nDof*(connections(i)-1);
        index(2*i-1:2*i) = base + (1:nDof);
    end

    K(index,index) = K(index,index) + k;
end


%% (C) Forming Global Equivalent Nodal Forces Vector

% initializing the variable (vector) to store global equivalent nodal force
F = zeros(nDof * nNode,1);

% (a) traction forces (finding their equivalent nodal forces)

if Model.nTractionForce > 0  % if there is at least one traction force, then:
    % extracting the tractions data from the Model structure
    tractions = Model.loading.tractionForces;
end
% loop over all the traction force data
for iTraction = 1 : Model.nTractionForce
    elem = tractions(iTraction,1);          % element number of i-th traction
    edge = tractions(iTraction,2);          % edge    number of i-th traction
    load = tractions(iTraction,3:end);      % the components of i-th traction [Tx, Ty]
    
    connections     = elements(elem,:);                 % a vector of element node numbers 
    edgeNodeIndex   = edgeNodes(nDim, nElemNode, edge); % a vector of THE edge nodes Ids
    edgeNodeNumbers = connections(edgeNodeIndex);       % a vector of THE edge nodes numbers 
    edgeCoord       = coordinates(edgeNodeNumbers,:);   % a matrix of coordinates of edge nodes
    
    % calling the function "elemEquivalentLoad" to calculate the element
    % equivalent nodal force
    f = elemEquivalentLoad(edgeCoord, load);
    
    % assembly: finding the indices of corresponding entries in global 
    % nodal forces and assemble the global nodal forces vector
    index = zeros(1, nDof*length(edgeNodeIndex));
    for i = 1:length(edgeNodeIndex)
        base = nDof*(edgeNodeNumbers(i)-1);
        index(2*i-1:2*i) = base + (1:nDof);
    end
    F(index) = F(index) + f;
    %
    % F(index,:) = ...
end

% (b) (concentrated) nodal forces 
if Model.nNodalForce > 0 % if there is at least one concentrated nodal force, then:
    % extracting the nodal forces data from the Model structure 
    nodalForces = Model.loading.nodalForces;
end

% loop over all the nodal forces data
% for iNodal = 1 : Model.nNodalForce
%
for iNodal = 1:Model.nNodalForce
    nd = Model.loading.nodalForces(iNodal,1);
    fx = Model.loading.nodalForces(iNodal,2);
    fy = Model.loading.nodalForces(iNodal,3);
    F(nDim*(nd-1)+1) = F(nDim*(nd-1)+1) + fx;
    F(nDim*(nd-1)+2) = F(nDim*(nd-1)+2) + fy;
end
%
% end

%% (D) Imposing Boundary Conditions and Solving the System

%{ 
boundary conditions components:
(1) bcIndex: a list (vector) of the degrees of freedom with boundary conditions
(2) bcValues: the prescribed displacements values
%}

bcIndex = (nDof*(Model.boundary(:,1)-1)+Model.boundary(:,2))';
bcValues = Model.boundary(:,3);


% calling the solution function to impose the boundary conditions and
% calculating the displacements and reaction forces
[displacements, reactionForces] = solution(K, F, bcIndex, bcValues);
displacements = reshape(displacements, nDim, nNode)';



%% (E) Calculating Strains and Stresses at Integration Points
[xiList, ~] = integrationPoints(nDim, nElemNode);
nPts = size(xiList,2);

stresses = zeros(3, nPts, nElem);
strains  = zeros(3, nPts, nElem);
% loop over all elements to calculate the elements' strains and stresses at
% integration points
for iElem = 1 : nElem
    connections = elements(iElem,:);           % a vector of element node numbers  
    elemCoord   = coordinates(connections,:);  % a matrix of element nodal coordinates 
    elemDisp    = displacements(connections,:);% a matrix of element nodal displacements  
    
    % generating a list of gauss quadrature integration points and weights
    % (in natural coordinates)
    
    % extracting the corresponding material properties
    E  = materials(iElem,1);
    nu = materials(iElem,2);
    D  = strainStressMatrix(E, nu, problem);
    
    for ip = 1:nPts
        xi = xiList(:,ip);
        [~, dNdxi] = shapeFunctions(nDim, nElemNode, xi);
        J    = dNdxi' * elemCoord;
        dNdx = J \ dNdxi';   % nDim×nElemNode

        % build B
        B = zeros(3, nDim*nElemNode);
        for i = 1:nElemNode
            Bi = [ dNdx(1,i),       0;
                   0,       dNdx(2,i);
                   dNdx(2,i), dNdx(1,i) ];
            B(:,2*i-1:2*i) = Bi;
        end

        % nodal disp → strain
        ue = reshape(elemDisp', [], 1);
        epsilon = B * ue;                   % [εx; εy; γxy]
        sigma  = D * epsilon;                   % stress

        strains(:,ip,iElem)  = epsilon;
        stresses(:,ip,iElem) = sigma;
    end
    
    
end


%% (F) Storing All the Results in an Output Structure Array Called "Results"


Results.KGlobal       = K;
Results.FGlobal       = F;
Results.displacements = displacements;
Results.reactions     = reactionForces;
Results.strains       = strains;
Results.stresses      = stresses;

end


