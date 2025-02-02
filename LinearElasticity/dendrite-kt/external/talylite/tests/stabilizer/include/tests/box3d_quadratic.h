/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //
#pragma once

#include <stabilizertest.h>

class Box3DQuadraticStabilizerTest : public SUPGTest {
 public:
  std::string name() override { return "3D Box (Quadratic)"; }
  ElemType elm_type() override { return kElem3dHexahedral; }
  GridType grid_type() override { return kGrid3dBox; }
  kBasisFunction basis_function() override { return BASIS_QUADRATIC; }
  int nsd() override { return 3; }

  ZEROPTV node_position(int i) override {
    double A0[27][3] = {};
    A0[0][0] = -1.0;
    A0[0][2] = 1.0;
    A0[1][0] = 1.0;
    A0[1][2] = 1.0;
    A0[2][0] = 1.0;
    A0[2][1] = 2.0;
    A0[2][2] = 1.0;
    A0[3][0] = -1.0;
    A0[3][1] = 2.0;
    A0[3][2] = 1.0;
    A0[4][0] = -1.0;
    A0[4][2] = 3.0;
    A0[5][0] = 1.0;
    A0[5][2] = 3.0;
    A0[6][0] = 1.0;
    A0[6][1] = 2.0;
    A0[6][2] = 3.0;
    A0[7][0] = -1.0;
    A0[7][1] = 2.0;
    A0[7][2] = 3.0;
    A0[8][0] = -1.0;
    A0[8][2] = 2.0;
    A0[9][0] = 1.0;
    A0[9][2] = 2.0;
    A0[10][0] = 1.0;
    A0[10][1] = 2.0;
    A0[10][2] = 2.0;
    A0[11][0] = -1.0;
    A0[11][1] = 2.0;
    A0[11][2] = 2.0;
    A0[12][2] = 1.0;
    A0[13][0] = 1.0;
    A0[13][1] = 1.0;
    A0[13][2] = 1.0;
    A0[14][1] = 2.0;
    A0[14][2] = 1.0;
    A0[15][0] = -1.0;
    A0[15][1] = 1.0;
    A0[15][2] = 1.0;
    A0[16][1] = 1.0;
    A0[16][2] = 1.0;
    A0[17][2] = 3.0;
    A0[18][0] = 1.0;
    A0[18][1] = 1.0;
    A0[18][2] = 3.0;
    A0[19][1] = 2.0;
    A0[19][2] = 3.0;
    A0[20][0] = -1.0;
    A0[20][1] = 1.0;
    A0[20][2] = 3.0;
    A0[21][1] = 1.0;
    A0[21][2] = 3.0;
    A0[22][2] = 2.0;
    A0[23][0] = 1.0;
    A0[23][1] = 1.0;
    A0[23][2] = 2.0;
    A0[24][1] = 2.0;
    A0[24][2] = 2.0;
    A0[25][0] = -1.0;
    A0[25][1] = 1.0;
    A0[25][2] = 2.0;
    A0[26][1] = 1.0;
    A0[26][2] = 2.0;

    assert(i >= 0 && i < 27);
    return ZEROPTV(A0[i][0], A0[i][1], A0[i][2]);
  }

  double SUPG(int itg_pt, int bf) override {
    double A0[27][27] = {};
    double t2 = sqrt(1.5E1);
    double t26 = t2*2.582644628099174E-3;
    double t3 = -t26+2.7E1/9.68E2;
    double t27 = t2*2.221074380165289E-2;
    double t4 = -t27+8.7E1/9.68E2;
    double t32 = t2*(1.5E1/9.68E2);
    double t5 = -t32+4.1E1/9.68E2;
    double t33 = t2*(1.3E1/9.68E2);
    double t6 = -t33-2.9E1/9.68E2;
    double t7 = t2*(3.0/1.21E2);
    double t8 = t7-5.7E1/4.84E2;
    double t9 = t2*(1.01E2/9.68E2);
    double t10 = t9-3.89E2/9.68E2;
    double t11 = t2*(2.3E1/9.68E2);
    double t12 = t11-7.9E1/9.68E2;
    double t13 = t2*(7.0/2.42E2);
    double t14 = t13-3.0/2.42E2;
    double t21 = t2*(3.1E1/2.42E2);
    double t15 = -t21+1.17E2/2.42E2;
    double t16 = t2*(9.0/4.84E2);
    double t17 = t16+1.2E1/1.21E2;
    double t18 = -t7+5.7E1/4.84E2;
    double t20 = t2*(1.05E2/4.84E2);
    double t19 = -t20+1.02E2/1.21E2;
    double t22 = -t13+3.0/2.42E2;
    double t23 = t21-1.17E2/2.42E2;
    double t25 = t2*(8.0/1.21E2);
    double t24 = -t25+3.8E1/1.21E2;
    double t28 = t27-8.7E1/9.68E2;
    double t29 = t26-2.7E1/9.68E2;
    double t30 = t2*(3.7E1/9.68E2);
    double t31 = t2*(1.73E2/9.68E2);
    double t34 = -t16-1.2E1/1.21E2;
    double t35 = t33+2.9E1/9.68E2;
    double t36 = t32-4.1E1/9.68E2;
    double t37 = -t11+7.9E1/9.68E2;
    double t38 = t20-1.02E2/1.21E2;
    double t39 = -t9+3.89E2/9.68E2;
    double t40 = t2*(5.0/7.26E2);
    double t41 = t2*(4.3E1/7.26E2);
    double t42 = t25-3.8E1/1.21E2;
    double t43 = t2*(1.0/4.4E1);
    double t44 = t43+9.0/4.4E1;
    double t45 = t2*(5.0/4.4E1);
    double t46 = t45-2.1E1/4.4E1;
    double t47 = t2*(1.0/1.1E1);
    double t48 = t2*(1.0/1.76E2);
    double t49 = -t43-9.0/4.4E1;
    double t50 = -t45+2.1E1/4.4E1;
    double t51 = -t47+2.0/1.1E1;
    double t52 = t2*(1.0/9.68E2);
    double t53 = t2*(1.7E1/9.68E2);
    double t54 = t30-1.503099173553719E-1;
    double t55 = t31-6.926652892561983E-1;
    double t56 = t41-2.9E1/1.21E2;
    double t57 = t40-9.0/1.21E2;
    double t58 = t53+7.283057851239669E-2;
    double t59 = t52+2.634297520661157E-2;
    double t60 = -t31+6.926652892561983E-1;
    double t61 = -t30+1.503099173553719E-1;
    double t62 = -t40+9.0/1.21E2;
    double t63 = -t41+2.9E1/1.21E2;
    double t64 = t2*(5.0/1.76E2);
    double t65 = t48-1.3E1/1.76E2;
    double t66 = t2*(1.7E1/1.76E2);
    double t67 = t2*(3.0/1.76E2);
    double t68 = t2*(2.0/3.3E1);
    double t69 = t2*(1.0/3.3E1);
    double t70 = t47-2.0/1.1E1;
    double t71 = t2+2.0;
    double t72 = 1.0/t71;
    double t73 = -t48+1.3E1/1.76E2;
    double t74 = t64+2.3E1/1.76E2;
    double t75 = t67-1.7E1/1.76E2;
    double t76 = t66-6.7E1/1.76E2;
    double t77 = t69+1.0/4.4E1;
    double t78 = t68-9.0/4.4E1;
    double t79 = t2*(1.0/2.4E1);
    double t80 = -t67+1.7E1/1.76E2;
    double t81 = -t66+6.7E1/1.76E2;
    double t82 = -t64-2.3E1/1.76E2;
    double t83 = -t68+9.0/4.4E1;
    double t84 = -t69-1.0/4.4E1;
    double t85 = -t53-7.283057851239669E-2;
    double t86 = -t52-2.634297520661157E-2;
    A0[0][0] = t2*(-1.7E1/9.68E2)-7.283057851239669E-2;
    A0[0][1] = t2*(-1.0/9.68E2)-2.634297520661157E-2;
    A0[0][2] = t4;
    A0[0][3] = t3;
    A0[0][4] = t3;
    A0[0][5] = t4;
    A0[0][6] = t55;
    A0[0][7] = t54;
    A0[0][8] = t6;
    A0[0][9] = t5;
    A0[0][10] = t10;
    A0[0][11] = t12;
    A0[0][12] = t17;
    A0[0][13] = t5;
    A0[0][14] = t8;
    A0[0][15] = t6;
    A0[0][16] = t14;
    A0[0][17] = t8;
    A0[0][18] = t10;
    A0[0][19] = t19;
    A0[0][20] = t12;
    A0[0][21] = t15;
    A0[0][22] = t14;
    A0[0][23] = t56;
    A0[0][24] = t15;
    A0[0][25] = t57;
    A0[0][26] = t24;
    A0[1][0] = t34;
    A0[1][1] = t17;
    A0[1][2] = t8;
    A0[1][3] = t18;
    A0[1][4] = t18;
    A0[1][5] = t8;
    A0[1][6] = t19;
    A0[1][7] = t38;
    A0[1][8] = t22;
    A0[1][9] = t14;
    A0[1][10] = t15;
    A0[1][11] = t23;
    A0[1][13] = t14;
    A0[1][15] = t22;
    A0[1][18] = t15;
    A0[1][20] = t23;
    A0[1][23] = t24;
    A0[1][25] = t42;
    A0[2][0] = t59;
    A0[2][1] = t58;
    A0[2][2] = t29;
    A0[2][3] = t28;
    A0[2][4] = t28;
    A0[2][5] = t29;
    A0[2][6] = t61;
    A0[2][7] = t60;
    A0[2][8] = t36;
    A0[2][9] = t35;
    A0[2][10] = t37;
    A0[2][11] = t39;
    A0[2][12] = t34;
    A0[2][13] = t35;
    A0[2][14] = t18;
    A0[2][15] = t36;
    A0[2][16] = t22;
    A0[2][17] = t18;
    A0[2][18] = t37;
    A0[2][19] = t38;
    A0[2][20] = t39;
    A0[2][21] = t23;
    A0[2][22] = t22;
    A0[2][23] = t62;
    A0[2][24] = t23;
    A0[2][25] = t63;
    A0[2][26] = t42;
    A0[3][13] = t65;
    A0[3][15] = t2*(-5.0/1.76E2)-2.3E1/1.76E2;
    A0[3][16] = t44;
    A0[3][18] = t2*(-1.7E1/1.76E2)+6.7E1/1.76E2;
    A0[3][20] = t2*(-3.0/1.76E2)+1.7E1/1.76E2;
    A0[3][21] = t46;
    A0[3][23] = t2*(-2.0/3.3E1)+9.0/4.4E1;
    A0[3][25] = t2*(-1.0/3.3E1)-1.0/4.4E1;
    A0[3][26] = t70;
    A0[4][13] = t44;
    A0[4][15] = t49;
    A0[4][18] = t46;
    A0[4][20] = t50;
    A0[4][23] = t72;
    A0[4][25] = t51;
    A0[5][13] = t74;
    A0[5][15] = t73;
    A0[5][16] = t49;
    A0[5][18] = t75;
    A0[5][20] = t76;
    A0[5][21] = t50;
    A0[5][23] = t77;
    A0[5][25] = t78;
    A0[5][26] = t51;
    A0[6][0] = t3;
    A0[6][1] = t4;
    A0[6][2] = t86;
    A0[6][3] = t85;
    A0[6][4] = t54;
    A0[6][5] = t55;
    A0[6][6] = t4;
    A0[6][7] = t3;
    A0[6][8] = t12;
    A0[6][9] = t10;
    A0[6][10] = t5;
    A0[6][11] = t6;
    A0[6][12] = t8;
    A0[6][13] = t5;
    A0[6][14] = t17;
    A0[6][15] = t6;
    A0[6][16] = t14;
    A0[6][17] = t19;
    A0[6][18] = t10;
    A0[6][19] = t8;
    A0[6][20] = t12;
    A0[6][21] = t15;
    A0[6][22] = t15;
    A0[6][23] = t56;
    A0[6][24] = t14;
    A0[6][25] = t57;
    A0[6][26] = t24;
    A0[7][0] = t18;
    A0[7][1] = t8;
    A0[7][2] = t17;
    A0[7][3] = t34;
    A0[7][4] = t38;
    A0[7][5] = t19;
    A0[7][6] = t8;
    A0[7][7] = t18;
    A0[7][8] = t23;
    A0[7][9] = t15;
    A0[7][10] = t14;
    A0[7][11] = t22;
    A0[7][13] = t14;
    A0[7][15] = t22;
    A0[7][18] = t15;
    A0[7][20] = t23;
    A0[7][23] = t24;
    A0[7][25] = t42;
    A0[8][0] = t28;
    A0[8][1] = t29;
    A0[8][2] = t58;
    A0[8][3] = t59;
    A0[8][4] = t60;
    A0[8][5] = t61;
    A0[8][6] = t29;
    A0[8][7] = t28;
    A0[8][8] = t39;
    A0[8][9] = t37;
    A0[8][10] = t35;
    A0[8][11] = t36;
    A0[8][12] = t18;
    A0[8][13] = t35;
    A0[8][14] = t34;
    A0[8][15] = t36;
    A0[8][16] = t22;
    A0[8][17] = t38;
    A0[8][18] = t37;
    A0[8][19] = t18;
    A0[8][20] = t39;
    A0[8][21] = t23;
    A0[8][22] = t23;
    A0[8][23] = t62;
    A0[8][24] = t22;
    A0[8][25] = t63;
    A0[8][26] = t42;
    A0[9][8] = t82;
    A0[9][9] = t65;
    A0[9][10] = t81;
    A0[9][11] = t80;
    A0[9][22] = t44;
    A0[9][23] = t83;
    A0[9][24] = t46;
    A0[9][25] = t84;
    A0[9][26] = t70;
    A0[10][8] = t49;
    A0[10][9] = t44;
    A0[10][10] = t46;
    A0[10][11] = t50;
    A0[10][23] = t72;
    A0[10][25] = t51;
    A0[11][8] = t73;
    A0[11][9] = t74;
    A0[11][10] = t75;
    A0[11][11] = t76;
    A0[11][22] = t49;
    A0[11][23] = t77;
    A0[11][24] = t50;
    A0[11][25] = t78;
    A0[11][26] = t51;
    A0[12][23] = t79-1.0/4.0;
    A0[12][25] = -t79-1.0/4.0;
    A0[12][26] = 1.0/2.0;
    A0[13][23] = 1.0/2.0;
    A0[13][25] = -1.0/2.0;
    A0[14][23] = t79+1.0/4.0;
    A0[14][25] = -t79+1.0/4.0;
    A0[14][26] = -1.0/2.0;
    A0[15][8] = t80;
    A0[15][9] = t81;
    A0[15][10] = t65;
    A0[15][11] = t82;
    A0[15][22] = t46;
    A0[15][23] = t83;
    A0[15][24] = t44;
    A0[15][25] = t84;
    A0[15][26] = t70;
    A0[16][8] = t50;
    A0[16][9] = t46;
    A0[16][10] = t44;
    A0[16][11] = t49;
    A0[16][23] = t72;
    A0[16][25] = t51;
    A0[17][8] = t76;
    A0[17][9] = t75;
    A0[17][10] = t74;
    A0[17][11] = t73;
    A0[17][22] = t50;
    A0[17][23] = t77;
    A0[17][24] = t49;
    A0[17][25] = t78;
    A0[17][26] = t51;
    A0[18][0] = t3;
    A0[18][1] = t4;
    A0[18][2] = t55;
    A0[18][3] = t54;
    A0[18][4] = t85;
    A0[18][5] = t86;
    A0[18][6] = t4;
    A0[18][7] = t3;
    A0[18][8] = t6;
    A0[18][9] = t5;
    A0[18][10] = t10;
    A0[18][11] = t12;
    A0[18][12] = t8;
    A0[18][13] = t10;
    A0[18][14] = t19;
    A0[18][15] = t12;
    A0[18][16] = t15;
    A0[18][17] = t17;
    A0[18][18] = t5;
    A0[18][19] = t8;
    A0[18][20] = t6;
    A0[18][21] = t14;
    A0[18][22] = t14;
    A0[18][23] = t56;
    A0[18][24] = t15;
    A0[18][25] = t57;
    A0[18][26] = t24;
    A0[19][0] = t18;
    A0[19][1] = t8;
    A0[19][2] = t19;
    A0[19][3] = t38;
    A0[19][4] = t34;
    A0[19][5] = t17;
    A0[19][6] = t8;
    A0[19][7] = t18;
    A0[19][8] = t22;
    A0[19][9] = t14;
    A0[19][10] = t15;
    A0[19][11] = t23;
    A0[19][13] = t15;
    A0[19][15] = t23;
    A0[19][18] = t14;
    A0[19][20] = t22;
    A0[19][23] = t24;
    A0[19][25] = t42;
    A0[20][0] = t28;
    A0[20][1] = t29;
    A0[20][2] = t61;
    A0[20][3] = t60;
    A0[20][4] = t59;
    A0[20][5] = t58;
    A0[20][6] = t29;
    A0[20][7] = t28;
    A0[20][8] = t36;
    A0[20][9] = t35;
    A0[20][10] = t37;
    A0[20][11] = t39;
    A0[20][12] = t18;
    A0[20][13] = t37;
    A0[20][14] = t38;
    A0[20][15] = t39;
    A0[20][16] = t23;
    A0[20][17] = t34;
    A0[20][18] = t35;
    A0[20][19] = t18;
    A0[20][20] = t36;
    A0[20][21] = t22;
    A0[20][22] = t22;
    A0[20][23] = t62;
    A0[20][24] = t23;
    A0[20][25] = t63;
    A0[20][26] = t42;
    A0[21][13] = t81;
    A0[21][15] = t80;
    A0[21][16] = t46;
    A0[21][18] = t65;
    A0[21][20] = t82;
    A0[21][21] = t44;
    A0[21][23] = t83;
    A0[21][25] = t84;
    A0[21][26] = t70;
    A0[22][13] = t46;
    A0[22][15] = t50;
    A0[22][18] = t44;
    A0[22][20] = t49;
    A0[22][23] = t72;
    A0[22][25] = t51;
    A0[23][13] = t75;
    A0[23][15] = t76;
    A0[23][16] = t50;
    A0[23][18] = t74;
    A0[23][20] = t73;
    A0[23][21] = t49;
    A0[23][23] = t77;
    A0[23][25] = t78;
    A0[23][26] = t51;
    A0[24][0] = t54;
    A0[24][1] = t55;
    A0[24][2] = t4;
    A0[24][3] = t3;
    A0[24][4] = t3;
    A0[24][5] = t4;
    A0[24][6] = t86;
    A0[24][7] = t85;
    A0[24][8] = t12;
    A0[24][9] = t10;
    A0[24][10] = t5;
    A0[24][11] = t6;
    A0[24][12] = t19;
    A0[24][13] = t10;
    A0[24][14] = t8;
    A0[24][15] = t12;
    A0[24][16] = t15;
    A0[24][17] = t8;
    A0[24][18] = t5;
    A0[24][19] = t17;
    A0[24][20] = t6;
    A0[24][21] = t14;
    A0[24][22] = t15;
    A0[24][23] = t56;
    A0[24][24] = t14;
    A0[24][25] = t57;
    A0[24][26] = t24;
    A0[25][0] = t38;
    A0[25][1] = t19;
    A0[25][2] = t8;
    A0[25][3] = t18;
    A0[25][4] = t18;
    A0[25][5] = t8;
    A0[25][6] = t17;
    A0[25][7] = t34;
    A0[25][8] = t23;
    A0[25][9] = t15;
    A0[25][10] = t14;
    A0[25][11] = t22;
    A0[25][13] = t15;
    A0[25][15] = t23;
    A0[25][18] = t14;
    A0[25][20] = t22;
    A0[25][23] = t24;
    A0[25][25] = t42;
    A0[26][0] = t60;
    A0[26][1] = t61;
    A0[26][2] = t29;
    A0[26][3] = t28;
    A0[26][4] = t28;
    A0[26][5] = t29;
    A0[26][6] = t58;
    A0[26][7] = t59;
    A0[26][8] = t39;
    A0[26][9] = t37;
    A0[26][10] = t35;
    A0[26][11] = t36;
    A0[26][12] = t38;
    A0[26][13] = t37;
    A0[26][14] = t18;
    A0[26][15] = t39;
    A0[26][16] = t23;
    A0[26][17] = t18;
    A0[26][18] = t35;
    A0[26][19] = t34;
    A0[26][20] = t36;
    A0[26][21] = t22;
    A0[26][22] = t23;
    A0[26][23] = t62;
    A0[26][24] = t22;
    A0[26][25] = t63;
    A0[26][26] = t42;

    assert(itg_pt >= 0 && itg_pt < 27);
    assert(bf >= 0 && bf < 27);
    return A0[itg_pt][bf];
  }
};

