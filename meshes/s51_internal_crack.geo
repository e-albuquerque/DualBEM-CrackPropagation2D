// =============================================================================
// s51_internal_crack.geo - Placa Retangular 3x1 com Trinca Central (Isida)
// Seção 5.1 da Dissertação (E = 30000, nu = 0, sigma = 10, 2a0 = 0.1)
// =============================================================================

lc = 0.1;

// --- Vértices do Contorno Externo ---
Point(1) = {0.0, 0.00, 0.0, lc};
Point(2) = {3.0, 0.00, 0.0, lc};
Point(3) = {3.0, 1.00, 0.0, lc};
Point(4) = {0.0, 1.00, 0.0, lc};

// --- Pontas da Trinca Central Vertical (x = 1.5, y de 0.45 a 0.55) ---
Point(5) = {1.5, 0.45, 0.0, lc}; // Ponta 2 (inferior)
Point(6) = {1.5, 0.55, 0.0, lc}; // Ponta 1 (superior)

// --- Arestas do Contorno (sentido anti-horário) ---
Line(1) = {1, 2}; // Base livre
Line(2) = {2, 3}; // Direita tracionada (tx = 10.0)
Line(3) = {3, 4}; // Topo livre
Line(4) = {4, 1}; // Esquerda engastada (ux = 0, uy = 0)

// --- Trinca Central ---
Line(5) = {5, 6};

// --- Discretização ---
Transfinite Curve {1, 3} = 24 + 1 Using Progression 1;
Transfinite Curve {2, 4} = 16 + 1 Using Progression 1;
Transfinite Curve {5}    = 4 + 1 Using Progression 1;

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// --- Grupos Físicos ---
Physical Curve("free_horiz", 1) = {1, 3};
Physical Curve("right_load", 2) = {2};
Physical Curve("left_fixed", 3) = {4};
Physical Curve("crack_1", 4)    = {5};

Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
