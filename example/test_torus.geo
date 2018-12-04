//
// Mesh parameters
lc = 0.0000001;
rafpol = 3;  // Default poloidal refinement
raftor = 5; // Default toroidal refinement
// Geometry parameters
h=0.2;
R0def = 3;
adef = 2;
//
DefineConstant[ n = {rafpol, Name "RPol", Min 2} ];
DefineConstant[ l = {raftor, Name "RTor", Min 1} ];
DefineConstant[ R0 = {R0def, Name "R", Min 2} ];
DefineConstant[ a = {adef, Name "a", Min 1} ];
L = h*a;




// Geometry
Point(1) = {R0,0,0,lc};
Point(2) = {R0+a,0,0,lc};
Point(3) = {R0,a,0,lc};
Point(4) = {R0,-a,0,lc};
Point(5) = {R0-a,0,0,lc};
Point(6) = {R0+L,0,0};
Point(7) = {R0,L,0};
Point(8) = {R0,-L,0};
Point(9) = {R0-L,0,0};
Circle(1) = {2,1,3};
Circle(2) = {3,1,5};
Circle(3) = {5,1,4};
Circle(4) = {4,1,2};
Line(5) = {6,7};
Line(6) = {7,9};
Line(7) = {9,8};
Line(8) = {8,6};
Line(9) = {6,2};
Line(10) = {7,3};
Line(11) = {8,4};
Line(12) = {9,5};

Transfinite Line "*" = n;

Line Loop(1) = {1,-10,-5,9};
Plane Surface(1) = {1};
Line Loop(2) = {2,-12,-6,10};
Plane Surface(2) = {2};
Line Loop(3) = {3,-11,-7,12};
Plane Surface(3) = {3};
Line Loop(4) = {4,-9,-8,11};
Plane Surface(4) = {4};
Line Loop(5) = {5,6,7,8};
Plane Surface(5) = {5};

Transfinite Surface "*";
Recombine Surface "*";

// Extrusion to torus
//out1[] = Extrude {0,0,1} {Surface{1,2,3,4,5}                                 ;Layers{l};Recombine;};*/

out1[] = Extrude {{0,1,0}, {0,0,0}, Pi/2.0} {Surface{1,2,3,4,5}                                 ;Layers{l};Recombine;};
out2[] = Extrude {{0,1,0}, {0,0,0}, Pi/2} {Surface{out1[0],out1[6],out1[12],out1[18],out1[24]};Layers{l};Recombine;};
out3[] = Extrude {{0,1,0}, {0,0,0}, Pi/2} {Surface{out2[0],out2[6],out2[12],out2[18],out2[24]};Layers{l};Recombine;};
out4[] = Extrude {{0,1,0}, {0,0,0}, Pi/2} {Surface{out3[0],out3[6],out3[12],out3[18],out3[24]};Layers{l};Recombine;};

/*out1[] = Extrude {{0,1,0}, {0,0,0}, 2*Pi/3} {Surface{1,2,3,4,5}                                 ;Layers{l};Recombine;};*/
/*out2[] = Extrude {{0,1,0}, {0,0,0}, 2*Pi/3} {Surface{out1[0],out1[6],out1[12],out1[18],out1[24]};Layers{l};Recombine;};*/
/*out3[] = Extrude {{0,1,0}, {0,0,0}, 2*Pi/3} {Surface{out2[0],out2[6],out2[12],out2[18],out2[24]};Layers{l};Recombine;};*/

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

//Physical Surface("in") = {1,2,3,4,5};

// Meshing...
Mesh 3;
SetOrder 2;
Mesh.ElementOrder=2 ;
//
Mesh.SecondOrderIncomplete=1;
