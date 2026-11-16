// =============================================================================
// s3132_inclined_crack.geo - Placa 10x20 com Trinca Central Inclinada a 45 deg
// Seção 3.1.3.2 da Dissertação (Validação de KI e KII em Modo Misto)
// =============================================================================

DefineConstant[ aw = {0.4, Min 0.1, Max 0.8, Step 0.1, Name "Razao 2a/w"} ];

a  = aw * 5.0; // Semi-comprimento da trinca (w = 10, meia-largura = 5)
dx = a / Sqrt(2.0);
dy = a / Sqrt(2.0);
lc = 1.0;

x_ini = 5.0 - dx;
y_ini = 10.0 - dy;
x_fim = 5.0 + dx;
y_fim = 10.0 + dy;

// --- Vértices do Contorno Externo ---
Point(1) = { 0.0,  0.0, 0.0, lc};
Point(2) = {10.0,  0.0, 0.0, lc};
Point(3) = {10.0, 20.0, 0.0, lc};
Point(4) = { 0.0, 20.0, 0.0, lc};

// --- Pontos da Trinca Graduada (6 elementos: 10%, 15%, 25%, 25%, 15%, 10%) ---
Point(10) = {x_ini,                           y_ini,                           0.0, lc};
Point(11) = {x_ini + 0.10 * (x_fim - x_ini),  y_ini + 0.10 * (y_fim - y_ini),  0.0, lc};
Point(12) = {x_ini + 0.25 * (x_fim - x_ini),  y_ini + 0.25 * (y_fim - y_ini),  0.0, lc};
Point(13) = {x_ini + 0.50 * (x_fim - x_ini),  y_ini + 0.50 * (y_fim - y_ini),  0.0, lc};
Point(14) = {x_ini + 0.75 * (x_fim - x_ini),  y_ini + 0.75 * (y_fim - y_ini),  0.0, lc};
Point(15) = {x_ini + 0.90 * (x_fim - x_ini),  y_ini + 0.90 * (y_fim - y_ini),  0.0, lc};
Point(16) = {x_fim,                           y_fim,                           0.0, lc};

// --- Curvas do Contorno (sentido anti-horário) ---
Line(1) = {1, 2}; // Base engastada (ux = 0, uy = 0)
Line(2) = {2, 3}; // Direita livre (8 elementos)
Line(3) = {3, 4}; // Topo tracionado (ty = +1.0)
Line(4) = {4, 1}; // Esquerda livre (8 elementos)

// --- Segmentos da Trinca ---
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 13};
Line(14) = {13, 14};
Line(15) = {14, 15};
Line(16) = {15, 16};

Transfinite Curve {1, 3} = 4 + 1 Using Progression 1;
Transfinite Curve {2, 4} = 8 + 1 Using Progression 1;
Transfinite Curve {11, 12, 13, 14, 15, 16} = 1 + 1 Using Progression 1;

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// --- Grupos Físicos ---
Physical Curve("bottom", 1)     = {1};
Physical Curve("free_sides", 2) = {2, 4};
Physical Curve("top", 3)        = {3};
Physical Curve("crack_1", 4)    = {11, 12, 13, 14, 15, 16};

Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
