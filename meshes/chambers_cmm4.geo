// =============================================================================
// chambers_cmm4.geo - Corpo de Prova Compact Mixed Mode (CMM4, alpha=110 deg)
// Seção 5.3.2 da Dissertação de Beatriz Borges Faria Fonseca (UFMG, 2023)
// Unidades: mm, N, MPa
// =============================================================================

lc = 1.0;

// --- Pontos do Contorno Externo e Entalhe ---
Point(1)  = { 0.000,  0.000, 0.0, lc};
Point(2)  = {18.150,  0.000, 0.0, lc};
Point(3)  = {18.150, 33.300, 0.0, lc};
Point(4)  = {19.140, 34.700, 0.0, lc}; // Raiz do entalhe / boca da pré-trinca
Point(5)  = {20.130, 33.300, 0.0, lc};
Point(6)  = {20.130,  0.000, 0.0, lc};
Point(7)  = {39.105,  0.000, 0.0, lc};
Point(8)  = {39.105, 20.740, 0.0, lc};
Point(9)  = {39.105, 25.740, 0.0, lc}; // Centro do furo direito (R = 5.0 mm)
Point(10) = {34.105, 25.740, 0.0, lc};
Point(11) = {39.105, 30.740, 0.0, lc};
Point(12) = {39.105, 58.740, 0.0, lc};
Point(13) = { 0.000, 58.740, 0.0, lc};

// --- Pontos do Furo Esquerdo de Carregamento (R = 5.0 mm) ---
Point(14) = {10.065, 25.740, 0.0, lc}; // Centro do furo esquerdo
Point(15) = {15.065, 25.740, 0.0, lc}; // Direita (0 deg)
Point(16) = {10.065, 20.740, 0.0, lc}; // Inferior (270 deg)
Point(17) = { 5.065, 25.740, 0.0, lc}; // Esquerda (180 deg)
Point(18) = {10.065, 30.740, 0.0, lc}; // Superior (90 deg)

// --- Ponta Inicial da Pré-Trinca (theta = -4 deg, L0 = 6.262 mm) ---
Point(19) = {19.304, 40.960, 0.0, lc};

// --- Curvas do Contorno Externo (sentido anti-horário) ---
Line(1)   = {1, 2};
Line(2)   = {2, 3};
Line(3)   = {3, 4};
Line(4)   = {4, 5};
Line(5)   = {5, 6};
Line(6)   = {6, 7};
Line(7)   = {7, 8};
Circle(8) = {8, 9, 10};   // Arco inferior do furo direito (sentido horário)
Circle(9) = {10, 9, 11};  // Arco superior do furo direito (sentido horário)
Line(10)  = {11, 12};
Line(11)  = {12, 13};
Line(12)  = {13, 1};

// --- Curvas do Furo Esquerdo (sentido horário: normal aponta para o furo) ---
Circle(13) = {15, 14, 16}; // 1º quadrante inferior-direito
Circle(14) = {16, 14, 17}; // 2º quadrante inferior-esquerdo
Circle(15) = {17, 14, 18}; // 3º quadrante superior-esquerdo (carregado)
Circle(16) = {18, 14, 15}; // 4º quadrante superior-direito  (carregado)

// --- Curva da Trinca Inicial ---
Line(17) = {4, 19};

// --- Discretização exata por curva (N elementos -> N + 1 pontos extremos) ---
Transfinite Curve {1, 6} = 5 + 1 Using Progression 1;
Transfinite Curve {3, 4} = 1 + 1 Using Progression 1;
Transfinite Curve {2, 5, 7, 8, 9, 10, 11, 12} = 10 + 1 Using Progression 1;
Transfinite Curve {13, 14, 15, 16} = 10 + 1 Using Progression 1;
Transfinite Curve {17} = 20 + 1 Using Progression 1;

// --- Superfície 2D para visualização no Gmsh ---
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Curve Loop(2) = {13, 14, 15, 16};
Plane Surface(1) = {1, 2};

// --- Grupos Físicos de Contorno e Trinca (lidos pelo preproc.py) ---
Physical Curve("free_boundary", 1)   = {1, 2, 3, 4, 5, 6, 9, 11, 12, 13, 14};
Physical Curve("right_symmetry", 2)  = {7, 10};
Physical Curve("right_hole_supp", 3) = {8};
Physical Curve("left_hole_q3", 4)    = {15};
Physical Curve("left_hole_q4", 5)    = {16};
Physical Curve("crack_1", 6)         = {17};

// --- Configuração de Elementos Quadráticos (3 nós: line3) ---
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
