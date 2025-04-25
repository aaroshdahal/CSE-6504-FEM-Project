// This code was created by pygmsh v0.7.5.
SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMin = 1e-8;
Mesh.CharacteristicLengthMax = 10.0;
Rectangle(1) = {0,0,0, 1, 1};
//+

//+
Field[1] = Box;
Field[1].VIn = .1;
Field[1].VOut = .1;
Field[1].XMin = 0.0;
Field[1].XMax = 100;
Field[1].YMin = 0.0;
Field[1].YMax = 100;
Field[1].ZMin = 0.0;
Field[1].ZMax = 0.0;
Field[1].Thickness = 5;
//+
Field[6] = Min;
Field[6].FieldsList = {1};
Background Field = 6;
//+
Mesh.ElementOrder = 1;
Mesh 2;
Mesh.MshFileVersion = 2.0;
//+
Save "tofu.msh";

