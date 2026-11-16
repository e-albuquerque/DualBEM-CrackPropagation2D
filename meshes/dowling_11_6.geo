// =============================================================================
// dowling_11_6.geo - Placa com Trinca Central sob Amplitude Variável (Rainflow)
// Seção 5.2.3 da Dissertação (Exemplo 11.6 de Dowling, Aço AISI 4340)
// Dimensões: largura W = 0.076 m, altura H = 0.152 m, trinca 2a0 = 0.002 m
// =============================================================================

lc = 0.01;

// --- Vértices do Contorno Externo ---
Point(1) = {0.000, 0.000, 0.0, lc};
Point(2) = {0.076, 0.000, 0.0, lc};
Point(3) = {0.076, 0.152, 0.0, lc};
Point(4) = {0.000, 0.152, 0.0, lc};

// --- Pontas da Trinca Interna Central (2a0 = 2.0 mm em y = 0.076 m) ---
Point(5) = {0.037, 0.076, 0.0, lc}; // Ponta 2 (esquerda)
Point(6) = {0.039, 0.076, 0.0, lc}; // Ponta 1 (direita)

// --- Arestas do Contorno (sentido anti-horário) ---
Line(1) = {1, 2}; // Base engastada (ux = 0, uy = 0)
Line(2) = {2, 3}; // Lateral direita livre
Line(3) = {3, 4}; // Topo tracionado (ty = 219.298245614035 MPa)
Line(4) = {4, 1}; // Lateral esquerda livre

// --- Trinca Central ---
Line(5) = {5, 6};

// --- Discretização da Malha ---
Transfinite Curve {1, 3} = 16 + 1 Using Progression 1;
Transfinite Curve {2, 4} = 32 + 1 Using Progression 1;
Transfinite Curve {5}    = 20 + 1 Using Progression 1;

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// --- Grupos Físicos ---
Physical Curve("bottom", 1)  = {1};
Physical Curve("right", 2)   = {2};
Physical Curve("top", 3)     = {3};
Physical Curve("left", 4)    = {4};
Physical Curve("crack_1", 5) = {5};

Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
