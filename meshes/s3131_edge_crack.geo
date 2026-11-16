// =============================================================================
// s3131_edge_crack.geo - Placa Quadrada 10x10 com Trinca de Borda (Modo I)
// Seção 3.1.3.1 da Dissertação (Validação de KI contra Civelek & Erdogan)
// =============================================================================

DefineConstant[ aw = {0.4, Min 0.1, Max 0.8, Step 0.1, Name "Razao a/w"} ];

w  = 10.0;
a  = aw * w;
lc = 1.0;

// --- Vértices do Contorno Externo ---
Point(1) = { 0.0,  0.0000, 0.0, lc};
Point(2) = {10.0,  0.0000, 0.0, lc};
Point(3) = {10.0,  4.6875, 0.0, lc};
Point(4) = {10.0,  5.3125, 0.0, lc};
Point(5) = {10.0, 10.0000, 0.0, lc};
Point(6) = { 0.0, 10.0000, 0.0, lc};
Point(7) = { 0.0,  5.0000, 0.0, lc}; // Boca da trinca na borda esquerda

// --- Pontos da Trinca Graduada (8 elementos: 20%, 20%, 15%, 15%, 10%, 10%, 5%, 5%) ---
Point(11) = {0.20 * a, 5.0, 0.0, lc};
Point(12) = {0.40 * a, 5.0, 0.0, lc};
Point(13) = {0.55 * a, 5.0, 0.0, lc};
Point(14) = {0.70 * a, 5.0, 0.0, lc};
Point(15) = {0.80 * a, 5.0, 0.0, lc};
Point(16) = {0.90 * a, 5.0, 0.0, lc};
Point(17) = {0.95 * a, 5.0, 0.0, lc};
Point(18) = {1.00 * a, 5.0, 0.0, lc}; // Ponta da trinca

// --- Curvas do Contorno (sentido anti-horário) ---
Line(1) = {1, 2}; // Base (ty = -1.0)
Line(2) = {2, 3}; // Direita inferior livre (7 elementos)
Line(3) = {3, 4}; // Direita central engastada (2 elementos)
Line(4) = {4, 5}; // Direita superior livre (7 elementos)
Line(5) = {5, 6}; // Topo (ty = +1.0)
Line(6) = {6, 7}; // Esquerda superior livre (8 elementos)
Line(7) = {7, 1}; // Esquerda inferior livre (8 elementos)

// --- Segmentos da Trinca (1 elemento quadrático por subsegmento) ---
Line(11) = {7,  11};
Line(12) = {11, 12};
Line(13) = {12, 13};
Line(14) = {13, 14};
Line(15) = {14, 15};
Line(16) = {15, 16};
Line(17) = {16, 17};
Line(18) = {17, 18};

Transfinite Curve {1, 5, 6, 7} = 8 + 1 Using Progression 1;
Transfinite Curve {2, 4}       = 7 + 1 Using Progression 1;
Transfinite Curve {3}          = 2 + 1 Using Progression 1;
Transfinite Curve {11, 12, 13, 14, 15, 16, 17, 18} = 1 + 1 Using Progression 1;

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(1) = {1};

// --- Grupos Físicos ---
Physical Curve("bottom", 1)        = {1};
Physical Curve("free_sides", 2)    = {2, 4, 6, 7};
Physical Curve("right_support", 3) = {3};
Physical Curve("top", 4)           = {5};
Physical Curve("crack_1", 5)       = {11, 12, 13, 14, 15, 16, 17, 18};

Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
