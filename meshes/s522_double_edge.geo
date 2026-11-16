// =============================================================================
// s522_double_edge.geo - Placa com Duas Trincas Simétricas de Borda (AISI 4340)
// Seção 5.2.2 da Dissertação (R = 1/3, sigma_max = 526 MPa, a0 = 1.0 mm)
// Dimensões: largura 2b = 0.076 m, altura 2h = 0.304 m
// =============================================================================

lc = 0.01;

// --- Vértices do Contorno Externo ---
Point(1)  = {0.000, 0.000, 0.0, lc};
Point(2)  = {0.076, 0.000, 0.0, lc};
Point(3)  = {0.076, 0.152, 0.0, lc}; // Boca da trinca 2 (direita)
Point(4)  = {0.076, 0.304, 0.0, lc};
Point(5)  = {0.000, 0.304, 0.0, lc};
Point(6)  = {0.000, 0.152, 0.0, lc}; // Boca da trinca 1 (esquerda)

// --- Pontas das Duas Trincas de Borda (a0 = 1.0 mm) ---
Point(7)  = {0.001, 0.152, 0.0, lc}; // Ponta da trinca 1 (esquerda)
Point(8)  = {0.075, 0.152, 0.0, lc}; // Ponta da trinca 2 (direita)

// --- Arestas do Contorno (sentido anti-horário) ---
Line(1) = {1, 2}; // Base engastada
Line(2) = {2, 3}; // Direita inferior livre
Line(3) = {3, 4}; // Direita superior livre
Line(4) = {4, 5}; // Topo tracionado (ty = +526 MPa)
Line(5) = {5, 6}; // Esquerda superior livre
Line(6) = {6, 1}; // Esquerda inferior livre

// --- Trincas de Borda ---
Line(7) = {6, 7}; // Trinca 1 (borda esquerda -> interior)
Line(8) = {3, 8}; // Trinca 2 (borda direita  -> interior)

// --- Discretização ---
Transfinite Curve {1, 2, 3, 4, 5, 6} = 38 + 1 Using Progression 1;
Transfinite Curve {7, 8}             = 10 + 1 Using Progression 1;

Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

// --- Grupos Físicos ---
Physical Curve("bottom", 1)     = {1};
Physical Curve("free_sides", 2) = {2, 3, 5, 6};
Physical Curve("top", 3)        = {4};
Physical Curve("crack_1", 4)    = {7};
Physical Curve("crack_2", 5)    = {8};

Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
