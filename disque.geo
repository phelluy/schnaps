rafxy =3;
rafz =1;
h=0.7;
Point(1) = {-h, -h, 0};
Point(2) = {h, -h, 0};
Point(3) = {h, h, 0};
Point(4) = {-h, h, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};
Point(5) = {0, 0, 0};
Point(6) = {1.414, 1.414, 0};
Point(7) = {-1.414, 1.414, 0};
Point(8) = {-1.414, -1.414, 0};
Point(9) = {1.414, -1.414, 0};
Circle(10) = {9, 5, 6};
Circle(11) = {6, 5, 7};
Circle(12) = {7, 5, 8};
Circle(13) = {8, 5, 9};
Line(14) = {2, 9};
Line(15) = {3, 6};
Line(16) = {4, 7};
Line(17) = {1, 8};
// est
Line Loop(18) = {14, 10, -15, -2};
Plane Surface(19) = {18};
// nord
Line Loop(20) = {15, 11, -16, -3};
Plane Surface(21) = {20};
// ouest
Line Loop(22) = {16, 12, -17, -4};
Plane Surface(23) = {22};
//sud
Line Loop(24) = {17, 13, -14, -1};
Plane Surface(25) = {24};
Extrude {0, 0, 1} {
  Surface{21, 23, 6, 19, 25};
}
ll[] = Line "*";
ss[] = Surface "*";
vv[] = Volume "*"; 
Transfinite Line {ll[]} = rafxy+1 Using Progression 1;
Transfinite Line{37,41,59,63,99,77,33,32} = rafz+1 Using Progression 1;
Transfinite Surface{ss[]};
Recombine Surface{ss[]};
Transfinite Volume{vv[]};
Recombine Volume{vv[]};
//Physical Volume(136) = {1, 2, 3, 4, 5};
Mesh.ElementOrder=2 ;
Mesh.SecondOrderIncomplete = 1 ;
