lx = 0.20;
ly = 0.02;
lz = 0.02;

nyz = 3; // number of celss in y and z directions
nx = 20;// number of cells in x direction
//+
Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {0.0, ly, 0.0, 1.0};
Point(3) = {0.0, ly, lz, 1.0};
Point(4) = {0.0, 0.0, lz, 1.0};


//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


//+
Curve Loop(1) = {2, 3, 4, 1}; Plane Surface(1) = {1};
//+
Transfinite Curve {1, 2, 3, 4} = nyz Using Progression 1;
//+
Transfinite Surface {1} = {4, 1, 2, 3}; Recombine Surface {1};
//+
Extrude {lx, 0, 0} {
  Surface{1}; Layers {nx}; Recombine;
}
//+
Physical Surface("XMIN", 27) = {1}; // Left boundary at x=0
//+
Physical Surface("XMAX", 28) = {26}; // right boundary at x = L

// Meshing 
// Open this file in gmsh and go to Mesh then press 3d , then press file 
// on top left corner and Expoert the mesh make sure to selelct all 
// elements option and choose the latest version 
 
