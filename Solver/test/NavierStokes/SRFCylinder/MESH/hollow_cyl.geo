SetFactory("OpenCASCADE");

// Parameters
R_inner = 1.0;  
R_outer = 4.0;    
L = 0.2;        // height  
zmin = 0.0;

n_radial = 12;    // Points in radial direction
n_circum = 12;   // Points in circumferential direction
n_axial = 1;     // Points in axial direction

// Create 8 points for inner circle (at 45° intervals)
Point(1) = {R_inner*Cos(0*Pi/4), R_inner*Sin(0*Pi/4), zmin, 1.0};
Point(2) = {R_inner*Cos(1*Pi/4), R_inner*Sin(1*Pi/4), zmin, 1.0};
Point(3) = {R_inner*Cos(2*Pi/4), R_inner*Sin(2*Pi/4), zmin, 1.0};
Point(4) = {R_inner*Cos(3*Pi/4), R_inner*Sin(3*Pi/4), zmin, 1.0};
Point(5) = {R_inner*Cos(4*Pi/4), R_inner*Sin(4*Pi/4), zmin, 1.0};
Point(6) = {R_inner*Cos(5*Pi/4), R_inner*Sin(5*Pi/4), zmin, 1.0};
Point(7) = {R_inner*Cos(6*Pi/4), R_inner*Sin(6*Pi/4), zmin, 1.0};
Point(8) = {R_inner*Cos(7*Pi/4), R_inner*Sin(7*Pi/4), zmin, 1.0};

// Create 8 points for outer circle
Point(11) = {R_outer*Cos(0*Pi/4), R_outer*Sin(0*Pi/4), zmin, 1.0};
Point(12) = {R_outer*Cos(1*Pi/4), R_outer*Sin(1*Pi/4), zmin, 1.0};
Point(13) = {R_outer*Cos(2*Pi/4), R_outer*Sin(2*Pi/4), zmin, 1.0};
Point(14) = {R_outer*Cos(3*Pi/4), R_outer*Sin(3*Pi/4), zmin, 1.0};
Point(15) = {R_outer*Cos(4*Pi/4), R_outer*Sin(4*Pi/4), zmin, 1.0};
Point(16) = {R_outer*Cos(5*Pi/4), R_outer*Sin(5*Pi/4), zmin, 1.0};
Point(17) = {R_outer*Cos(6*Pi/4), R_outer*Sin(6*Pi/4), zmin, 1.0};
Point(18) = {R_outer*Cos(7*Pi/4), R_outer*Sin(7*Pi/4), zmin, 1.0};

// Center point
Point(100) = {0, 0, zmin, 1.0};

// Inner circle arcs
Circle(1) = {1, 100, 2};
Circle(2) = {2, 100, 3};
Circle(3) = {3, 100, 4};
Circle(4) = {4, 100, 5};
Circle(5) = {5, 100, 6};
Circle(6) = {6, 100, 7};
Circle(7) = {7, 100, 8};
Circle(8) = {8, 100, 1};

// Outer circle arcs
Circle(11) = {11, 100, 12};
Circle(12) = {12, 100, 13};
Circle(13) = {13, 100, 14};
Circle(14) = {14, 100, 15};
Circle(15) = {15, 100, 16};
Circle(16) = {16, 100, 17};
Circle(17) = {17, 100, 18};
Circle(18) = {18, 100, 11};

// Radial lines
Line(21) = {1, 11};
Line(22) = {2, 12};
Line(23) = {3, 13};
Line(24) = {4, 14};
Line(25) = {5, 15};
Line(26) = {6, 16};
Line(27) = {7, 17};
Line(28) = {8, 18};

// Create 8 surfaces
Curve Loop(1) = {1, 22, -11, -21};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 23, -12, -22};
Plane Surface(2) = {2};
Curve Loop(3) = {3, 24, -13, -23};
Plane Surface(3) = {3};
Curve Loop(4) = {4, 25, -14, -24};
Plane Surface(4) = {4};
Curve Loop(5) = {5, 26, -15, -25};
Plane Surface(5) = {5};
Curve Loop(6) = {6, 27, -16, -26};
Plane Surface(6) = {6};
Curve Loop(7) = {7, 28, -17, -27};
Plane Surface(7) = {7};
Curve Loop(8) = {8, 21, -18, -28};
Plane Surface(8) = {8};

// Transfinite meshing
Transfinite Curve{1,2,3,4,5,6,7,8} = n_circum;
Transfinite Curve{11,12,13,14,15,16,17,18} = n_circum;
Transfinite Curve{21,22,23,24,25,26,27,28} = n_radial;
Transfinite Surface{1,2,3,4,5,6,7,8};
Recombine Surface{1,2,3,4,5,6,7,8};

// Extrude
Extrude {0, 0, L} {
  Surface{1,2,3,4,5,6,7,8};
  Layers{n_axial};
  Recombine;
}

// Groups
Physical Surface("inner_wall") = {9,14,18,22,26,30,34,38};
Physical Surface("outer_wall") = {11,16,20,24,28,32,36,39};
Physical Surface("top") = {13,17,21,25,29,33,37,40};
Physical Surface("bottom") = {1,2,3,4,5,6,7,8};
Physical Volume("fluid") = {1,2,3,4,5,6,7,8};


//Mesh.ElementOrder = 2;
//Mesh.SecondOrderLinear = 0;
//Mesh.HighOrderOptimize = 2;


