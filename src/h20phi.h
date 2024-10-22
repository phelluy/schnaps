/* phi[0] = (-1 + z) * (-1 + y) * (-1 + x) * (2 * x + 2 * y + 2 * z - 1); */
/* phi[1] = x * (-1 + z) * (-1 + y) * (-2 * y - 2 * z + 2 * x - 1); */
/* phi[2] = -x * y * (-1 + z) * (2 * y - 2 * z - 3 + 2 * x); */
/* phi[3] = -y * (-1 + z) * (-1 + x) * (2 * x + 2 * z + 1 - 2 * y); */
/* phi[4] = -z * (-1 + y) * (-1 + x) * (2 * x + 2 * y - 2 * z + 1); */
/* phi[5] = -x * z * (-1 + y) * (-2 * y + 2 * z - 3 + 2 * x); */
/* phi[6] = x * y * z * (2 * y + 2 * z - 5 + 2 * x); */
/* phi[7] = y * z * (-1 + x) * (2 * x - 2 * z + 3 - 2 * y); */
/* phi[8] = -4 * x * (-1 + z) * (-1 + y) * (-1 + x); */
/* phi[9] = -4 * y * (-1 + z) * (-1 + y) * (-1 + x); */
/* phi[10] = -4 * z * (-1 + z) * (-1 + y) * (-1 + x); */
/* phi[11] = 4 * x * y * (-1 + z) * (-1 + y); */
/* phi[12] = 4 * x * z * (-1 + z) * (-1 + y); */
/* phi[13] = 4 * x * y * (-1 + z) * (-1 + x); */
/* phi[14] = -4 * x * y * z * (-1 + z); */
/* phi[15] = 4 * y * z * (-1 + z) * (-1 + x); */
/* phi[16] = 4 * x * z * (-1 + y) * (-1 + x); */
/* phi[17] = 4 * y * z * (-1 + y) * (-1 + x); */
/* phi[18] = -4 * x * y * z * (-1 + y); */
/* phi[19] = -4 * x * y * z * (-1 + x); */

/* gradphi[0][0] = (-1 + z) * (-1 + y) * (2 * y + 2 * z + 4 * x - 3); */
/* gradphi[0][1] = (-1 + z) * (-1 + x) * (2 * x + 2 * z - 3 + 4 * y); */
/* gradphi[0][2] = (-1 + y) * (-1 + x) * (2 * x + 2 * y - 3 + 4 * z); */
/* gradphi[1][0] = (-1 + z) * (-1 + y) * (-2 * y - 2 * z - 1 + 4 * x); */
/* gradphi[1][1] = x * (-1 + z) * (-2 * z + 1 - 4 * y + 2 * x); */
/* gradphi[1][2] = x * (-1 + y) * (-2 * y + 1 - 4 * z + 2 * x); */
/* gradphi[2][0] = -y * (-1 + z) * (2 * y - 2 * z - 3 + 4 * x); */
/* gradphi[2][1] = -x * (-1 + z) * (-2 * z - 3 + 4 * y + 2 * x); */
/* gradphi[2][2] = -x * y * (2 * y - 4 * z - 1 + 2 * x); */
/* gradphi[3][0] = -y * (-1 + z) * (-2 * y + 2 * z - 1 + 4 * x); */
/* gradphi[3][1] = -(-1 + z) * (-1 + x) * (2 * x + 2 * z + 1 - 4 * y); */
/* gradphi[3][2] = -y * (-1 + x) * (2 * x - 1 + 4 * z - 2 * y); */
/* gradphi[4][0] = -z * (-1 + y) * (2 * y - 2 * z - 1 + 4 * x); */
/* gradphi[4][1] = -z * (-1 + x) * (2 * x - 2 * z - 1 + 4 * y); */
/* gradphi[4][2] = -(-1 + y) * (-1 + x) * (2 * x - 4 * z + 2 * y + 1); */
/* gradphi[5][0] = -z * (-1 + y) * (-2 * y + 2 * z + 4 * x - 3); */
/* gradphi[5][1] = -x * z * (2 * z - 4 * y - 1 + 2 * x); */
/* gradphi[5][2] = -x * (-1 + y) * (-2 * y - 3 + 4 * z + 2 * x); */
/* gradphi[6][0] = y * z * (2 * y + 2 * z + 4 * x - 5); */
/* gradphi[6][1] = x * z * (2 * z + 4 * y - 5 + 2 * x); */
/* gradphi[6][2] = x * y * (2 * y + 4 * z - 5 + 2 * x); */
/* gradphi[7][0] = y * z * (-2 * y - 2 * z + 1 + 4 * x); */
/* gradphi[7][1] = z * (-1 + x) * (2 * x - 2 * z + 3 - 4 * y); */
/* gradphi[7][2] = y * (-1 + x) * (2 * x - 2 * y + 3 - 4 * z); */
/* gradphi[8][0] = -4 * (-1 + z) * (-1 + y) * (2 * x - 1); */
/* gradphi[8][1] = -4 * x * (-1 + z) * (-1 + x); */
/* gradphi[8][2] = -4 * x * (-1 + y) * (-1 + x); */
/* gradphi[9][0] = -4 * y * (-1 + z) * (-1 + y); */
/* gradphi[9][1] = -4 * (-1 + z) * (2 * y - 1) * (-1 + x); */
/* gradphi[9][2] = -4 * y * (-1 + y) * (-1 + x); */
/* gradphi[10][0] = -4 * z * (-1 + z) * (-1 + y); */
/* gradphi[10][1] = -4 * z * (-1 + z) * (-1 + x); */
/* gradphi[10][2] = -4 * (2 * z - 1) * (-1 + y) * (-1 + x); */
/* gradphi[11][0] = 4 * y * (-1 + z) * (-1 + y); */
/* gradphi[11][1] = 4 * x * (-1 + z) * (2 * y - 1); */
/* gradphi[11][2] = 4 * x * y * (-1 + y); */
/* gradphi[12][0] = 4 * z * (-1 + z) * (-1 + y); */
/* gradphi[12][1] = 4 * x * z * (-1 + z); */
/* gradphi[12][2] = 4 * x * (2 * z - 1) * (-1 + y); */
/* gradphi[13][0] = 4 * y * (-1 + z) * (2 * x - 1); */
/* gradphi[13][1] = 4 * x * (-1 + z) * (-1 + x); */
/* gradphi[13][2] = 4 * x * y * (-1 + x); */
/* gradphi[14][0] = -4 * y * z * (-1 + z); */
/* gradphi[14][1] = -4 * x * z * (-1 + z); */
/* gradphi[14][2] = -4 * x * y * (2 * z - 1); */
/* gradphi[15][0] = 4 * y * z * (-1 + z); */
/* gradphi[15][1] = 4 * z * (-1 + z) * (-1 + x); */
/* gradphi[15][2] = 4 * y * (2 * z - 1) * (-1 + x); */
/* gradphi[16][0] = 4 * z * (-1 + y) * (2 * x - 1); */
/* gradphi[16][1] = 4 * x * z * (-1 + x); */
/* gradphi[16][2] = 4 * x * (-1 + y) * (-1 + x); */
/* gradphi[17][0] = 4 * y * z * (-1 + y); */
/* gradphi[17][1] = 4 * z * (2 * y - 1) * (-1 + x); */
/* gradphi[17][2] = 4 * y * (-1 + y) * (-1 + x); */
/* gradphi[18][0] = -4 * y * z * (-1 + y); */
/* gradphi[18][1] = -4 * x * z * (2 * y - 1); */
/* gradphi[18][2] = -4 * x * y * (-1 + y); */
/* gradphi[19][0] = -4 * y * z * (2 * x - 1); */
/* gradphi[19][1] = -4 * x * z * (-1 + x); */
/* gradphi[19][2] = -4 * x * y * (-1 + x); */

schnaps_real t1 = -1 + z;
schnaps_real t2 = -1 + y;
schnaps_real t3 = t1 * t2;
schnaps_real t4 = 2 * y;
schnaps_real t5 = 2 * z;
schnaps_real t6 = 4 * x;
schnaps_real t9 = -1 + x;
schnaps_real t10 = t1 * t9;
schnaps_real t11 = 2 * x;
schnaps_real t12 = 4 * y;
schnaps_real t15 = t2 * t9;
schnaps_real t16 = 4 * z;
schnaps_real t24 = x * t1;
schnaps_real t27 = x * t2;
schnaps_real t33 = y * t1;
schnaps_real t38 = x * y;
schnaps_real t48 = y * t9;
schnaps_real t54 = z * t2;
schnaps_real t57 = z * t9;
schnaps_real t67 = x * z;
schnaps_real t75 = y * z;
schnaps_real t94 = t11 - 1;
schnaps_real t98 = 4 * t24 * t9;
schnaps_real t100 = 4 * t27 * t9;
schnaps_real t104 = 4 * t33 * t2;
schnaps_real t105 = t4 - 1;
schnaps_real t111 = 4 * y * t2 * t9;
schnaps_real t114 = z * t1;
schnaps_real t116 = 4 * t114 * t2;
schnaps_real t118 = 4 * t114 * t9;
schnaps_real t119 = t5 - 1;
schnaps_real t128 = 4 * t38 * t2;
schnaps_real t132 = 4 * t67 * t1;
schnaps_real t141 = 4 * t38 * t9;
schnaps_real t145 = 4 * t75 * t1;
schnaps_real t158 = 4 * t67 * t9;
schnaps_real t162 = 4 * t75 * t2;

gradphi[0][0] = t3 * (t4 + t5 - 3 + t6);
gradphi[0][1] = t10 * (t11 + t5 - 3 + t12);
gradphi[0][2] = t15 * (t11 + t4 - 3 + t16);
gradphi[0][3] = t3 * t9 * (t11 + t4 + t5 - 1);
gradphi[1][0] = t3 * (-t4 - t5 - 1 + t6);
gradphi[1][1] = t24 * (-t5 + 1 + t11 - t12);
gradphi[1][2] = t27 * (-t4 + 1 + t11 - t16);
gradphi[1][3] = t24 * t2 * (-t4 - t5 + t11 - 1);
gradphi[2][0] = -t33 * (t4 - t5 - 3 + t6);
gradphi[2][1] = -t24 * (-t5 - 3 + t11 + t12);
gradphi[2][2] = -t38 * (t4 - 1 + t11 - t16);
gradphi[2][3] = -t38 * t1 * (t4 - t5 - 3 + t11);
gradphi[3][0] = -t33 * (-t4 + t5 - 1 + t6);
gradphi[3][1] = -t10 * (t11 + t5 + 1 - t12);
gradphi[3][2] = -t48 * (t11 - 1 - t4 + t16);
gradphi[3][3] = -t33 * t9 * (t11 + t5 - t4 + 1);
gradphi[4][0] = -t54 * (t4 - t5 - 1 + t6);
gradphi[4][1] = -t57 * (t11 - t5 - 1 + t12);
gradphi[4][2] = -t15 * (t11 + 1 + t4 - t16);
gradphi[4][3] = -t54 * t9 * (t11 + t4 - t5 + 1);
gradphi[5][0] = -t54 * (-t4 + t5 - 3 + t6);
gradphi[5][1] = -t67 * (t5 - 1 + t11 - t12);
gradphi[5][2] = -t27 * (-t4 - 3 + t11 + t16);
gradphi[5][3] = -t67 * t2 * (-t4 + t5 - 3 + t11);
gradphi[6][0] = t75 * (t4 + t5 - 5 + t6);
gradphi[6][1] = t67 * (t5 - 5 + t11 + t12);
gradphi[6][2] = t38 * (t4 - 5 + t11 + t16);
gradphi[6][3] = t38 * z * (t4 + t5 - 5 + t11);
gradphi[7][0] = t75 * (-t4 - t5 + 1 + t6);
gradphi[7][1] = t57 * (t11 - t5 + 3 - t12);
gradphi[7][2] = t48 * (t11 - t4 + 3 - t16);
gradphi[7][3] = t75 * t9 * (t11 - t5 - t4 + 3);
gradphi[8][0] = -4 * t3 * t94;
gradphi[8][1] = -t98;
gradphi[8][2] = -t100;
gradphi[8][3] = -4 * t24 * t15;
gradphi[9][0] = -t104;
gradphi[9][1] = -4 * t1 * t105 * t9;
gradphi[9][2] = -t111;
gradphi[9][3] = -4 * t33 * t15;
gradphi[10][0] = -t116;
gradphi[10][1] = -t118;
gradphi[10][2] = -4 * t119 * t2 * t9;
gradphi[10][3] = -4 * t114 * t15;
gradphi[11][0] = t104;
gradphi[11][1] = 4 * t24 * t105;
gradphi[11][2] = t128;
gradphi[11][3] = 4 * t38 * t3;
gradphi[12][0] = t116;
gradphi[12][1] = t132;
gradphi[12][2] = 4 * x * t119 * t2;
gradphi[12][3] = 4 * t67 * t3;
gradphi[13][0] = 4 * t33 * t94;
gradphi[13][1] = t98;
gradphi[13][2] = t141;
gradphi[13][3] = 4 * t38 * t10;
gradphi[14][0] = -t145;
gradphi[14][1] = -t132;
gradphi[14][2] = -4 * t38 * t119;
gradphi[14][3] = -4 * t38 * t114;
gradphi[15][0] = t145;
gradphi[15][1] = t118;
gradphi[15][2] = 4 * y * t119 * t9;
gradphi[15][3] = 4 * t75 * t10;
gradphi[16][0] = 4 * t54 * t94;
gradphi[16][1] = t158;
gradphi[16][2] = t100;
gradphi[16][3] = 4 * t67 * t15;
gradphi[17][0] = t162;
gradphi[17][1] = 4 * z * t105 * t9;
gradphi[17][2] = t111;
gradphi[17][3] = 4 * t75 * t15;
gradphi[18][0] = -t162;
gradphi[18][1] = -4 * t67 * t105;
gradphi[18][2] = -t128;
gradphi[18][3] = -4 * t38 * t54;
gradphi[19][0] = -4 * t75 * t94;
gradphi[19][1] = -t158;
gradphi[19][2] = -t141;
gradphi[19][3] = -4 * t38 * t57;
