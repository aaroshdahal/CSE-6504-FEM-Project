function k = elemStiff(elemCoord, materials, problem)
%{
This function takes the coordinates of the element (elemCoord matrix), a
vector of the element material properties [E, nu, t], and the problem type
('plane stress / plane strain) as the input arguments and the conduct the
numerical integration to calculate of the element stiffness matrix
%}

[nElemNode,nDim] = size(elemCoord);

% extracting the associated material properties
E  = materials(1);
nu = materials(2);
t  = materials(3);

% calling the "strainStressMatrix" function to calculate the [D] matrix
D = strainStressMatrix (E, nu, problem);


% generating a list of gauss quadrature integration points and weights
% (in natural coordinates)
[xiList, wList] = integrationPoints(nDim,nElemNode);
nPoints = size(xiList,2);

% initializing the variable to store the [B] and [k] matrices values
B = zeros(3, nDim*nElemNode);
k = zeros(nDim*nElemNode);

% loop over all integration points to conduct the numerical integration and
% obtaining the element stiffness matrix
for ip = 1:nPoints
    xi = xiList(:,ip);
    w  = wList(ip);

    [~, dNdxi] = shapeFunctions(nDim, nElemNode, xi); % gets dNdxi (number of shape functions (nElemNode) x ndim) at xi

    % Jacobian
    J    = dNdxi' * elemCoord; % J =dX/d eta = dNiXi/d etaj

    detJ = det(J);
    dNdx = J \ dNdxi';             % physical‐space derivatives : dNdX = inv(J) * dNdxi'(same thing)

    % build B (building it for each shapefunction separately, then compiling it)
    for i = 1:nElemNode
        Bi = [ dNdx(1,i),        0;
               0,        dNdx(2,i);
               dNdx(2,i), dNdx(1,i) ];
        B(:, 2*i-1:2*i) = Bi;
    end

    k = k + (B' * D * B) * t * detJ * w;

end

end
