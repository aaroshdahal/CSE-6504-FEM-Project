function writeVTKall(nodes, elems, U, strains, stresses, filename)
%WRITEVTKALL  Export mesh + displacement, strain & stress to ASCII VTK
%
%   nodes      : nNode×2 array of (x,y) coords
%   elems      : nElem×nElemNode connectivity (1-based)
%   U          : nNode×2 displacement vectors
%   strains    : 3×nPts×nElem (ε11, ε22, ε12 at each Gauss pt)
%   stresses   : 3×nPts×nElem (σ11, σ22, σ12 at each Gauss pt)
%   filename   : name of the .vtk file to write

% ensure 3D coords
[nNode, nDim] = size(nodes);
if nDim==2, nodes(:,3)=0; end

nElem     = size(elems,1);
nElemNode = size(elems,2);

% average strains & stresses per element
strainAvg = zeros(nElem,3);
stressAvg = zeros(nElem,3);
for e = 1:nElem
  strainAvg(e,:) = mean(strains(:,:,e), 2)';
  stressAvg(e,:) = mean(stresses(:,:,e),2)';
end

fid = fopen(filename,'w');
fprintf(fid,'# vtk DataFile Version 2.0\nMATLAB FEM Results\nASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');

% points
fprintf(fid,'POINTS %d float\n', nNode);
fprintf(fid,'%f %f %f\n', nodes');

% cells
cellSize = nElem*(nElemNode+1);
fprintf(fid,'\nCELLS %d %d\n', nElem, cellSize);
for e=1:nElem
  fprintf(fid,'%d %s\n', nElemNode, sprintf('%d ', elems(e,:)-1));
end

% cell types (9 = VTK_QUAD)
fprintf(fid,'\nCELL_TYPES %d\n', nElem);
fprintf(fid,'%d\n', 9*ones(nElem,1));

% point data: displacement vectors
fprintf(fid,'\nPOINT_DATA %d\n', nNode);
fprintf(fid,'VECTORS displacement float\n');
for i=1:nNode
  fprintf(fid,'%f %f %f\n', U(i,1), U(i,2), 0);
end

% cell data: strains then stresses
fprintf(fid,'\nCELL_DATA %d\n', nElem);
fprintf(fid,'SCALARS strain_e11 float 1\nLOOKUP_TABLE default\n');
fprintf(fid,'%f\n', strainAvg(:,1));
fprintf(fid,'SCALARS strain_e22 float 1\nLOOKUP_TABLE default\n');
fprintf(fid,'%f\n', strainAvg(:,2));
fprintf(fid,'SCALARS strain_e12 float 1\nLOOKUP_TABLE default\n');
fprintf(fid,'%f\n', strainAvg(:,3));
fprintf(fid,'SCALARS stress_s11 float 1\nLOOKUP_TABLE default\n');
fprintf(fid,'%f\n', stressAvg(:,1));
fprintf(fid,'SCALARS stress_s22 float 1\nLOOKUP_TABLE default\n');
fprintf(fid,'%f\n', stressAvg(:,2));
fprintf(fid,'SCALARS stress_s12 float 1\nLOOKUP_TABLE default\n');
fprintf(fid,'%f\n', stressAvg(:,3));

fclose(fid);
end
