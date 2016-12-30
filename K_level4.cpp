#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"Header.h"


void rollingB(double Bt[][6], double B[][24]){	//20160609debug

		for (int j = 0; j < 6; j++){
			for (int k = 0; k < 24; k++){
				Bt[k][j] = B[j][k];
			}
		}

}

void productBtD(double Bt[][6], double D[][6], double BD[][6]){	//20160609debug

	for (int j = 0; j < 24; j++){
		for (int k = 0; k < 6; k++){
			for (int l = 0; l < 6; l++){
				BD[j][k] += Bt[j][l] * D[l][k];
			}
		}
	}

}

void jacob(double J[][3], double dndr[][8], double nd[][8],double *r){		//20160616debug

	dndr[0][0] = -(1 - r[1])*(1 - r[2]) * 0.125;
	dndr[1][0] = -(1 - r[0])*(1 - r[2]) * 0.125;
	dndr[2][0] = -(1 - r[0])*(1 - r[1]) * 0.125;
	
	
	dndr[0][1] = (1 - r[1])*(1 - r[2]) * 0.125;
	dndr[1][1] = -(1 + r[0])*(1 - r[2]) * 0.125;
	dndr[2][1] = -(1 + r[0])*(1 - r[1]) * 0.125;
	
	dndr[0][2] = (1 + r[1])*(1 - r[2]) * 0.125;
	dndr[1][2] = (1 + r[0])*(1 - r[2]) * 0.125;
	dndr[2][2] = -(1 + r[0])*(1 + r[1]) * 0.125;
	
	dndr[0][3] = -(1 + r[1])*(1 - r[2]) * 0.125;
	dndr[1][3] = (1 - r[0])*(1 - r[2]) * 0.125;
	dndr[2][3] = -(1 - r[0])*(1 + r[1]) * 0.125;
	
	dndr[0][4] = -(1 - r[1])*(1 + r[2]) * 0.125;
	dndr[1][4] = -(1 - r[0])*(1 + r[2]) * 0.125;
	dndr[2][4] = (1 - r[0])*(1 - r[1]) * 0.125;

	dndr[0][5] = (1 - r[1])*(1 + r[2]) * 0.125;
	dndr[1][5] = -(1 + r[0])*(1 + r[2]) * 0.125;
	dndr[2][5] = (1 + r[0])*(1 - r[1]) * 0.125;

	dndr[0][6] = (1 + r[1])*(1 + r[2]) * 0.125;
	dndr[1][6] = (1 + r[0])*(1 + r[2]) * 0.125;
	dndr[2][6] = (1 + r[0])*(1 + r[1]) * 0.125;

	dndr[0][7] = -(1 + r[1])*(1 + r[2]) * 0.125;
	dndr[1][7] = (1 - r[0])*(1 + r[2]) * 0.125;
	dndr[2][7] = (1 - r[0])*(1 + r[1]) * 0.125;



			for (int l = 0; l < 8; l++){
				J[0][0] += dndr[0][l] * nd[0][l];
				J[0][1] += dndr[0][l] * nd[1][l];
				J[0][2] += dndr[0][l] * nd[2][l];

				J[1][0] += dndr[1][l] * nd[0][l];
				J[1][1] += dndr[1][l] * nd[1][l];
				J[1][2] += dndr[1][l] * nd[2][l];

				J[2][0] += dndr[2][l] * nd[0][l];
				J[2][1] += dndr[2][l] * nd[1][l];
				J[2][2] += dndr[2][l] * nd[2][l];
	
			}


			
	

}

void inv_jacob(double invJ[][3], double J[][3]){	//20160609debug

	//for (int j = 0; j < 3; j++){
	//	for (int k = 0; k < 3; k++){
	//		invJ[j][j] = 1;
	//	}
	//}

	//double inv_pivot, temp;

	//for (int ipv = 0; ipv < 3; ipv++){

	//	inv_pivot = 1.0 / J[ipv][ipv];

	//	for (int j = 0; j < 3; j++){
	//		J[ipv][j] *= inv_pivot;
	//		invJ[ipv][j] *= inv_pivot;
	//	}

	//	for (int j = 0; j < 3; j++){
	//		if (j != ipv){
	//			temp = J[j][ipv];
	//			for (int k = 0; k < 3; k++){
	//				J[j][k] -= temp*J[ipv][k];
	//				invJ[j][k] -= temp*invJ[ipv][k];
	//			}
	//		}
	//	}
	//}
	double detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
				- J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
				+ J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
	
	invJ[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) / detJ;
	invJ[0][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) / detJ;
	invJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / detJ;

	invJ[1][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) / detJ;
	invJ[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / detJ;
	invJ[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) / detJ;

	invJ[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) / detJ;
	invJ[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) / detJ;
	invJ[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / detJ;
}

void DnDx(double dndx[][8], double invJ[][3],double dndr[][8]){	//20160609debug

	for (int j = 0; j < 3; j++){
		for (int k = 0; k < 8; k++){
			for (int l = 0; l < 3; l++){
				dndx[j][k] += invJ[j][l] * dndr[l][k];
			}
		}
	}

}


void productBDB(double BD[][6], double B[][24], double Ke[][24]){	//20160609debug

	for (int j = 0; j < 24; j++){
		for (int k = 0; k < 24; k++){
			for (int l = 0; l < 6; l++){
				Ke[j][k] += BD[j][l] * B[l][k];
			}
		}
	}

}