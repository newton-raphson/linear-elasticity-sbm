Point(1) = {0, 0, 0, cl};
Point(2) = {1, 0, 0, cl};
Point(3) = {1, 1, 0, cl};
Point(4) = {0, 1, 0, cl};
Line(1) = {1, 2};
Line(2) = {3, 2};
Line(3) = {4, 3};
Line(4) = {4, 1};
Line Loop(5) = {3, 2, -1, -4};
Plane Surface(6) = {5};
Point(5) = {1, 1, 1, cl};
Point(6) = {0, 1, 1, cl};
Point(7) = {1, 0, 1, cl};
Point(8) = {0, 0, 1, cl};
Line(7) = {6, 5};
Line(8) = {5, 7};
Line(9) = {7, 8};
Line(10) = {8, 6};
Line(11) = {8, 1};
Line(12) = {4, 6};
Line(13) = {5, 3};
Line(14) = {7, 2};
Line Loop(15) = {11, 1, -14, 9};
Plane Surface(16) = {15};
Line Loop(17) = {4, -11, 10, -12};
Plane Surface(18) = {17};
Line Loop(19) = {7, 8, 9, 10};
Plane Surface(20) = {19};
Line Loop(21) = {13, 2, -14, -8};
Plane Surface(22) = {21};
Line Loop(23) = {12, 7, 13, -3};
Plane Surface(24) = {23};
Surface Loop(25) = {24, 18, 6, 22, 16, 20};
Volume(26) = {25};

Physical Surface(1) = {18};
Physical Surface(2) = {22};
Physical Surface(3) = {16};
Physical Surface(4) = {24};
Physical Surface(5) = {6};
Physical Surface(6) = {20};
Physical Volume(7) = {26};

//Transfinite Line "*" = 1;
//Transfinite Surface "*";
//Recombine Surface "*";
//Transfinite Volume "*";
