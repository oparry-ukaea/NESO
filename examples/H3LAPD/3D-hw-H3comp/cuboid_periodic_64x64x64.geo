//=============================== Parameters ==================================
// Lengths and resolutions in each dimension
one_over_rho_s0 = 773.7422539156483;
xsize = 18.2e-3*one_over_rho_s0;
ysize = 18.2e-3*one_over_rho_s0;
zsize = 10.0*one_over_rho_s0;
nx = 64;
ny = 64;
nz = 64;
//=============================================================================

// Create a line in the x-direction of length <xsize>, with <nx> divisions
Point(1) = {0, 0, 0, 0.01};
Point(2) = {xsize, 0, 0, 0.01};
Line(1) = {1, 2};
Transfinite Line(1) = nx+1;

// Extrude split line into meshed square/rectangle
sq = Extrude {0,ysize,0} {Curve{1}; Layers{ny}; Recombine;};

// Extrude square/rectangle into a cuboid
cbd = Extrude {0,0,zsize} {Surface{sq[1]}; Layers{nz}; Recombine;};

// Define physical volume, surfaces for BCs
// Domain
Physical Volume(0) = {cbd[1]};
// Low-x side
Physical Surface(1) = {cbd[5]};
// High-x side
Physical Surface(2) = {cbd[3]};
// Low-y side
Physical Surface(3) = {cbd[2]};
// High-y side
Physical Surface(4) = {cbd[4]};
// Low-z side
Physical Surface(5) = {sq[1]};
// High-z side
Physical Surface(6) = {cbd[0]};