// This code was created by pygmsh v0.7.5.
SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMin = 1e-8;
Mesh.CharacteristicLengthMax = 10.0;
Rectangle(1) = {0,0,0, 20, 20};
Disk(2) = {0, 0, 0, 1}; 
BooleanDifference{Surface{1}; Delete;}{Surface{2}; Delete;}
//+
Field[4] = Box;
Field[4].VIn = .02;
Field[4].VOut = 1;
Field[4].XMin = 0.0;
Field[4].XMax = 1.5;
Field[4].YMin = 0.0;
Field[4].YMax = 1.5;
Field[4].ZMin = 0.0;
Field[4].ZMax = 0.0;
Field[4].Thickness = .1;
//+
Field[1] = Box;
Field[1].VIn = .05;
Field[1].VOut = 1;
Field[1].XMin = 0.0;
Field[1].XMax = 2;
Field[1].YMin = 0.0;
Field[1].YMax = 2;
Field[1].ZMin = 0.0;
Field[1].ZMax = 0.0;
Field[1].Thickness = .1;
//+
Field[2] = Box;
Field[2].VIn = .15;
Field[2].VOut = 1;
Field[2].XMin = 0.0;
Field[2].XMax = 4;
Field[2].YMin = 0.0;
Field[2].YMax = 4;
Field[2].ZMin = 0.0;
Field[2].ZMax = 0.0;
Field[2].Thickness = .1;
//+
Field[3] = Box;
Field[3].VIn = .4;
Field[3].VOut = 1;
Field[3].XMin = 0.0;
Field[3].XMax = 8;
Field[3].YMin = 0.0;
Field[3].YMax = 8;
Field[3].ZMin = 0.0;
Field[3].ZMax = 0.0;
Field[3].Thickness = .1;
//+
Field[6] = Min;
Field[6].FieldsList = {1,2,3,4};
Background Field = 6;
//+
Mesh.ElementOrder = 1;
Mesh 2;
Mesh.MshFileVersion = 2.0;
//+
Save "stress_conc.msh";

