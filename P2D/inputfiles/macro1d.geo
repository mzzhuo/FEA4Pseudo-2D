
L_a = 100; 
L_s = 25;
L_c = 100;
lcar = L_a/10;

Point(1) = {0, 0, 0, lcar};
Point(2) = {L_a, 0, 0, lcar};
Point(3) = {L_a+L_s, 0, 0, lcar};
Point(4) = {L_a+L_s+L_c, 0, 0, lcar};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};

Transfinite Line { 1 } = 21 Using Progression 0.95;
Transfinite Line { 2 } = 6  Using Bump 0.5;
Transfinite Line { 3 } = 21 Using Progression 1.05;

Physical Line("anode")  = {1};
Physical Line("separator")  = {2};
Physical Line("cathode")  = {3};

Physical Point("left")  = {1};
Physical Point("right")  = {4};

Mesh 1;
//RefineMesh;
