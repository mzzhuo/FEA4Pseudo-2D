
R = 1; 
lcar = R/10;

Point(1) = {0, 0, 0, lcar};
Point(2) = {R, 0, 0, lcar};

Line(1) = {1, 2};

Transfinite Line { 1 } = 101 Using Progression 1.0;

Physical Line("particle1d")  = {1};

Mesh 2;
//RefineMesh;
