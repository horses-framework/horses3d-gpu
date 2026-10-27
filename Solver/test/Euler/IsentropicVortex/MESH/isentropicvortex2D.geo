Mesh.MshFileVersion = 4.1;

// ----------- BEGIN INPUT ------------
Nx = 24;  
Ny = 24;
Nz = 1;
lx = 30.0;
ly = 30.0;
lz = 2.0;
// ----------- END INPUT --------------


// Centered domain
Point(1) = {-lx/2, -ly/2, -lz/2, 0.1};

line[] = Extrude {lx, 0, 0} {
  Point{1}; Layers{Nx};
};

surface[] = Extrude {0.0, ly, 0.0} {
  Line{line[1]}; Layers{Ny}; Recombine;
};

volume[] = Extrude {0.0, 0.0, lz} {
  Surface{surface[1]}; Layers{Nz}; Recombine;
};


// Physical entities (NOTE: IDs may change after shifting, see comment below)
Physical Surface ("left")   = {26};
Physical Surface ("right")  = {18};
Physical Surface ("top")    = {22};
Physical Surface ("bottom") = {14};
Physical Surface ("front")  = {27};
Physical Surface ("back")   = {5};
Physical Volume  ("fluid")  = {volume[1]};

Mesh.RecombineAll = 1;
