cl__1 = 0.04;
cl__2 = 0.01;
cl__3 = 0.004;
Point(1) = {0, 0, 0, 0.04};
Point(2) = {1, 0, 0, 0.04};
Point(3) = {0, 1, 0, 0.04};
Point(4) = {0, 0, 1, 0.04};
Point(5) = {1, 1, 0, 0.04};
Point(6) = {0, 1, 1, 0.04};
Point(7) = {1, 0, 1, 0.04};
Point(8) = {1, 1, 1, 0.04};
Line(1) = {3, 5};
Line(2) = {5, 8};
Line(3) = {8, 6};
Line(4) = {6, 3};
Line(5) = {3, 1};
Line(6) = {1, 2};
Line(7) = {2, 7};
Line(8) = {7, 4};
Line(9) = {4, 1};
Line(10) = {6, 4};
Line(11) = {8, 7};
Line(12) = {5, 2};
Line Loop(14) = {1, 2, 3, 4};
Plane Surface(14) = {14};
Line Loop(16) = {2, 11, -7, -12};
Plane Surface(16) = {16};
Line Loop(18) = {1, 12, -6, -5};
Plane Surface(18) = {18};
Line Loop(20) = {4, 5, -9, -10};
Plane Surface(20) = {20};
Line Loop(22) = {3, 10, -8, -11};
Plane Surface(22) = {22};
Line Loop(24) = {8, 9, 6, 7};
Plane Surface(24) = {24};
Surface Loop(26) = {14, 16, 18, 20, 22, 24};
Volume(26) = {26};
Physical Surface(1) = {20};
Physical Surface(2) = {16};
Physical Surface(3) = {24};
Physical Surface(4) = {14};
Physical Surface(5) = {18};
Physical Surface(6) = {22};
Physical Volume(27) = {26};
