// =============================================================================
// s521_single_edge.geo - Placa com Uma Trinca de Borda (Aço AISI 4340)
// Seção 5.2.1 da Dissertação (R = 1/3, sigma_max = 526 MPa, a0 = 2.0 mm)
// Dimensões: largura W = 0.076 m, altura H = 0.152 m
// =============================================================================

lc = 0.01;

// --- Vértices do Contorno Externo ---
// Nota: A aresta esquerda é dividida na boca da trinca (0, 0.076) e a aresta
// direita possui um segmento central engastado ao redor de y = 0.076 m para
// impedir movimento de corpo rígido sob tração simétrica superior/inferior.
Point(1) = {0.000, 0.000000, 0.0, lc};
Point(2) = {0.076, 0.000000, 0.0, lc};
Point(3) = {0.076, 0.073625, 0.0, lc};
Point(4) = {0.076, 0.078375, 0.0, lc};
Point(5) = {0.076, 0.152000, 0.0, lc};
Point(6) = {0.000, 0.152000, 0.0, lc};
Point(7) = {0.000, 0.076000, 0.0, lc}; // Boca da trinca na borda esquerda
Point(8) = {0.002, 0.076000, 0.0, lc}; // Ponta inicial da trinca (a0 = 2 mm)

// --- Arestas do Contorno (sentido anti-horário) ---
Line(1) = {1, 2}; // Base tracionada (ty = -526 MPa)
Line(2) = {2, 3}; // Direita inferior livre (15 elementos)
Line(3) = {3, 4}; // Direita central engastada (2 elementos)
Line(4) = {4, 5}; // Direita superior livre (15 elementos)
Line(5) = {5, 6}; // Topo tracionado (ty = +526 MPa)
Line(6) = {6, 7}; // Esquerda superior livre (16 elementos)
Line(7) = {7, 1}; // Esquerda inferior livre (16 elementos)

// --- Trinca de Borda ---
Line(8) = {7, 8};

// --- Discretização ---
Transfinite Curve {1, 5, 6, 7} = 16 + 1 Using Progression 1;
Transfinite Curve {2, 4}       = 15 + 1 Using Progression 1;
Transfinite Curve {3}          = 2 + 1 Using Progression 1;
Transfinite Curve {8}          = 20 + 1 Using Progression 1;

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(1) = {1};

// --- Grupos Físicos ---
Physical Curve("bottom", 1)        = {1};
Physical Curve("free_sides", 2)    = {2, 4, 6, 7};
Physical Curve("right_support", 3) = {3};
Physical Curve("top", 4)           = {5};
Physical Curve("crack_1", 5)       = {8};

Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
