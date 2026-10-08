// Parameters
r1 = 1.0;  // Inner arc radius
r2 = 2.0;  // Middle arc radius
r3 = 4.0;  // Outer arc radius
n_radial = 18;    // Number of nodes in radial direction
n_circumferential = 10;  // Number of nodes in circumferential direction
n_special = 26;   // Number of nodes in special region (-15° to +15°)
angle = 15 * Pi / 180;  // 15 degrees in radians

// Center point
Point(1) = {0, 0, 0};

// Points at ±15° on inner arc (radius 1)
Point(2) = {r1 * Cos(angle), r1 * Sin(angle), 0};
Point(3) = {r1 * Cos(angle), -r1 * Sin(angle), 0};

// Points at ±15° on middle arc (radius 2)
Point(4) = {r2 * Cos(angle), r2 * Sin(angle), 0};
Point(5) = {r2 * Cos(angle), -r2 * Sin(angle), 0};

// Points at ±15° on outer arc (radius 3)
Point(6) = {r3 * Cos(angle), r3 * Sin(angle), 0};
Point(7) = {r3 * Cos(angle), -r3 * Sin(angle), 0};

// Points at 180° (negative x-axis)
Point(8) = {-r1, 0, 0};  // Inner arc at 180°
Point(9) = {-r2, 0, 0};  // Middle arc at 180°
Point(10) = {-r3, 0, 0}; // Outer arc at 180°

// Inner circle arcs
//Circle(1) = {2, 1, 3};   // +15° to -15° (short arc, 30°)
Circle(2) = {3, 1, 8};   // -15° to 180° (first half of long arc)
Circle(3) = {8, 1, 2};   // 180° to +15° (second half of long arc)

// Middle circle arcs
Circle(4) = {4, 1, 5};   // +15° to -15° (short arc)
Circle(5) = {5, 1, 9};   // -15° to 180°
Circle(6) = {9, 1, 4};   // 180° to +15°

// Outer circle arcs
Circle(7) = {6, 1, 7};   // +15° to -15° (short arc)
Circle(8) = {7, 1, 10};  // -15° to 180°
Circle(9) = {10, 1, 6};  // 180° to +15°

// Radial lines - inner to middle
Line(10) = {2, 4};   // +15°
Line(11) = {3, 5};   // -15°
Line(12) = {8, 9};   // 180°

// Radial lines - middle to outer
Line(13) = {4, 6};   // +15°
Line(14) = {5, 7};   // -15°
Line(15) = {9, 10};  // 180°

// Define surfaces - outer ring only (middle to outer)
Curve Loop(1) = {5, 15, -8, -14};  // Special region part 1: -15° to 180°
Plane Surface(1) = {1};

Curve Loop(2) = {6, 13, -9, -15};  // Special region part 2: 180° to +15°
Plane Surface(2) = {2};

Curve Loop(3) = {4, 14, -7, -13};  // Outer region: +15° to -15°
Plane Surface(3) = {3};

// Inner ring surface: 180° to +15° (between inner and middle arcs)
Curve Loop(4) = {3, 10, -6, -12};
Plane Surface(4) = {4};

// Inner ring surface: -15° to 180° (between inner and middle arcs)
Curve Loop(5) = {2, 12, -5, -11};
Plane Surface(5) = {5};


// Transfinite meshing - special region
Transfinite Curve {5, 6, 8, 9} = n_special;
Transfinite Curve {2, 3} = n_special;  // Add inner arc curves
Transfinite Curve {13, 14, 15} = n_radial;
Transfinite Curve {10, 11, 12} = n_radial;  // Add inner-to-middle radial lines

// Transfinite meshing - outer region
Transfinite Curve {4, 7} = n_circumferential;
Transfinite Curve {1} = n_circumferential;  // Add inner short arc

Transfinite Surface {1, 2, 3, 4, 5};
Recombine Surface {1, 2, 3, 4, 5};

// Extrusion parameters
h = 0.2;  // Height of extrusion

// Extrude all surfaces in z-direction with 1 element
Extrude {0, 0, h} {
  Surface{1, 2, 3, 4, 5}; Layers{1}; Recombine;
}

// Physical groups for boundary conditions
Physical Surface("bottom") = {1, 2, 3, 4, 5}; // Bottom surfaces (z=0)
Physical Surface("top") = {37,59,81,103,125};  // Top surfaces (z=h)
Physical Surface("inner_wall") = {112,124,68,94,90};  
Physical Surface("outer_wall") = {76,32,54};  
Physical Volume("fluid") = {1, 2, 3, 4, 5};  // All volumes


Mesh.RecombineAll = 1;


