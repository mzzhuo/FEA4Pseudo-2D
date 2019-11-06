
L_a = 100; 
L_s = 25;
L_c = 100;
h   = 10;

lcar = L_a/10;


Point(1) = {0, 0, 0, lcar};
Point(2) = {L_a, 0, 0, lcar};
Point(3) = {L_a, h, 0, lcar};
Point(4) = {0, h, 0, lcar};

Line(1) = {1, 2};
Line(2) = {4, 3};
Line(3) = {1, 4};
Line(4) = {2, 3};

Point(5) = {L_a+L_s, 0, 0, lcar};
Point(6) = {L_a+L_s, h, 0, lcar};

Line(5) = {2, 5};
Line(6) = {3, 6};
Line(7) = {5, 6};

Point(7) = {L_a+L_s+L_c, 0, 0, lcar};
Point(8) = {L_a+L_s+L_c, h, 0, lcar};

Line(8) = {5, 7};
Line(9) = {6, 8};
Line(10) = {7, 8};

Transfinite Line { 1,2 } = 21 Using Progression 0.95;
Transfinite Line { 5,6 } = 6  Using Bump 0.5;
Transfinite Line { 8,9 } = 21 Using Progression 1.05;
Transfinite Line { 3,10 } = 3;

Line Loop(1) = {1, 4, -2, -3};
Plane Surface(1) = {1};

Line Loop(2) = {5, 7, -6, -4};
Plane Surface(2) = {2};

Line Loop(3) = {8, 10, -9, -7};
Plane Surface(3) = {3};

// comment these two lines for Tri3; otherwises Quad4
Transfinite Surface {1,2,3};
Recombine Surface {1,2,3};

Physical Surface("anode")  = {1};
Physical Surface("separator")  = {2};
Physical Surface("cathode")  = {3};


Physical Line("left")  = {3};
Physical Line("right") = {10};


Mesh 2;