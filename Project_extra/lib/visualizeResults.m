function nodeVM = visualizeResults(Results, scaleFactor)
%VISUALIZERESULTS  plot undeformed & deformed mesh with Von Mises contour
%
%   nodeVM = visualizeResults(Results, scaleFactor)
%
%   - Results.coordinates: nNode×2
%   - Results.elements:    nElem×nElemNode
%   - Results.displacements: nNode×2
%   - Results.stresses:    3×nPts×nElem  (holds ε; your code already swapped names)
%
coords    = Results.coordinates;
elems     = Results.elements;
U         = Results.displacements;

% infer element geometry/order
[nNode, nDim]    = size(coords);
nElemNode        = size(elems,2);
[~, nPts, nElem] = size(Results.stresses);

% accumulate nodal Von Mises by weighted averaging
nodeVM    = zeros(nNode,1);
connCount = zeros(nNode,1);
[xiList,~] = integrationPoints(nDim,nElemNode);    % :contentReference[oaicite:0]{index=0}&#8203;:contentReference[oaicite:1]{index=1}

for e = 1:nElem
  conn    = elems(e,:);
  elemC   = coords(conn,:);
  elemU   = U(conn,:);
  % build D for this element (you can also cache D in Results if you like)
  E  = Results.materials(e,1);  nu = Results.materials(e,2);
  D  = strainStressMatrix(E, nu, Results.problem);

  for ip = 1:nPts
    xi = xiList(:,ip);
    [N, dNdxi] = shapeFunctions(nDim,nElemNode, xi);  % :contentReference[oaicite:2]{index=2}&#8203;:contentReference[oaicite:3]{index=3}
    J    = dNdxi' * elemC;
    dNdx = J \ dNdxi';
    % build B
    B = zeros(3, nDim*nElemNode);
    for i=1:nElemNode
      Bi = [ dNdx(1,i),        0;
             0,        dNdx(2,i);
             dNdx(2,i), dNdx(1,i) ];
      B(:,2*i-1:2*i) = Bi;
    end

    ue = reshape(elemU',[],1);
    eps = B*ue;           % engineering shear
    eps(3)=eps(3)/2;      % true shear
    sig = D*eps;          % σ

    % Von Mises
    s11=sig(1); s22=sig(2); s12=sig(3);
    vm = sqrt(s11^2 - s11*s22 + s22^2 + 3*s12^2);

    % scatter to nodes
    for i=1:nElemNode
      nodeVM(conn(i))    = nodeVM(conn(i))    + N(i)*vm;
      connCount(conn(i)) = connCount(conn(i)) + N(i);
    end
  end
end

nodeVM = nodeVM ./ connCount;

% now plot
figure; hold on; axis equal off
patch('Faces',elems,'Vertices',coords, ...
      'FaceColor','none','EdgeColor',[.8 .8 .8]);
coords_def = coords + scaleFactor*U;
patch('Faces',elems,'Vertices',coords_def, ...
      'FaceVertexCData',nodeVM,'FaceColor','interp','EdgeColor','k');
colorbar; title('Von Mises & deformed mesh');
end
