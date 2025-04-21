function r = elemEquivalentLoad(edgeCoord, load)
%{
This function takes the coordinates of an edge of an element, and the
componenets of the applied traction force [Tx, Ty], and calculates the
equivalent nodal forces by conducting numerical integration
%}

[nEdgeNode, nDim] = size(edgeCoord);
nDof = length(load);

% initializing the variable to store the equivalent nodal forces
r = zeros(nDof*nEdgeNode, 1);

% calling the "integrationPoints" function to generate a list of
% Guass quadrature integration points and weights
[xiList, wList] = integrationPoints(nDim-1,nEdgeNode); %nDim-1 cause we just want the shape function matrix not the derivative matrix

% loop over integration points to conduct the numerical integration and
% calculate the equivalent nodal forces
for ip = 1:length(wList)
    xi = xiList(ip);
    w  = wList(ip);

    [N, dNdxi] = shapeFunctions(1, nEdgeNode, xi);
    % Jacobian for the curve
    Jline = dNdxi' * edgeCoord;           % 1×nDim
    detJ  = norm(Jline);                  % length element

    % assemble
    for i = 1:nEdgeNode
        idx = nDim*(i-1) + (1:nDim); % calculating index for global orce vector
        r(idx) = r(idx) + N(i) * load(:) * detJ * w; % summing over diff gauss points
    end
end

end

