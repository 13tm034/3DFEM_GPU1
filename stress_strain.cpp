#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"Header.h"


void strain_stress_calc(node *no, element *el, material *m, int E, int N, int M){

	for (int i = 0; i < E; i++){			//要素iについて
		if (el[i].rm == false){
			//point
			double p[5] = {};
			p[0] = -0.9061798459;
			p[1] = -0.5384693101;
			p[2] = 0;
			p[3] = 0.5384693101;
			p[4] = 0.9061798459;

			//weight
			double weight[5];
			weight[0] = 0.2369268851;
			weight[1] = 0.4786286704;
			weight[2] = 0.5688888889;
			weight[3] = 0.4786286704;
			weight[4] = 0.2369268851;

			double B[6][24] = {};
			double DB[6][24] = {};
			//double D[6][6] = {};
			//Dmatrix(D, m, i, el);

			for (int j = 0; j < 5; j++){
				for (int k = 0; k < 5; k++){
					for (int l = 0; l < 5; l++){

						double r[3] = {};
						r[0] = p[l];
						r[1] = p[k];
						r[2] = p[j];


						double w[3] = {};
						w[0] = weight[l];
						w[1] = weight[k];
						w[2] = weight[j];




						double b[6][24] = {};			//Bマトリクス	
						double J[3][3] = {};			//jacobマトリクス
						Bmatrix(no, el, b, i, r, J);

						double detJ = 0;				//jacobのdetの値
						double *det;
						det = &detJ;
						det_jacob(J, det);
						double D[6][6] = {};
						Dmatrix(D, m, i, el);

						addB(B, b, w, detJ);
						addDB(DB, b, D, w, detJ);

					}
				}
			}

			//for (int j = 0; j < 6; j++){
			//	for (int k = 0; k < 24; k++){
			//		printf("%lf,", B[j][k]);
			//	}
			//	printf("\n");
			//}
			//for (int j = 0; j < 6; j++){
			//	for (int k = 0; k < 24; k++){
			//		printf("%lf,", DB[j][k]);
			//	}
			//	printf("\n");
			//}

			//system("PAUSE");
			calc_strain(B, el, no, i);
			calc_stress(DB, el, no, i);
			mises_strain(el, i, m);

		}
	}
}




void mises_strain(element *el, int i,material *m){


	double s[6] = {};

	for (int j = 0; j < 6; j++){
		s[j] = el[i].stress[j];
	}

	double sig1 = (s[0] - s[1])*(s[0] - s[1]);
	double sig2 = (s[1] - s[2])*(s[1] - s[2]);
	double sig3 = (s[2] - s[0])*(s[2] - s[0]);
	double sig4 = 2 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
	el[i].mises = sqrt((sig1 + sig2 + sig3 + 3 * sig4)/ 2);
	//printf("mises=%lf\n", el[i].mises);

}

void addB(double B[][24], double b[][24], double *w,double detJ){


	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 24; j++){
			B[i][j] += detJ*b[i][j] * w[0] * w[1] * w[2];
		}
	}



}


void addDB(double DB[][24], double b[][24], double D[][6], double *w, double detJ){

	double db[6][24] = {};


	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 24; j++){
			for (int k = 0; k < 6; k++){
				db[i][j] += D[i][k] * b[k][j];
			}
		}
	}


	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 24; j++){
			DB[i][j] += detJ*db[i][j] * w[0] * w[1] * w[2];
		
		}

	}



}



void calc_strain(double B[][24], element *el, node *no, int i){


	int n[8] = {};					//全体節点番号

	for (int j = 0; j < 8; j++){
		n[j] = el[i].node[j];
	}

	double nd[24] = {};			//一要素の各節点の座標値

	for (int j = 0; j < 8; j++){
		nd[j+j+j] = no[n[j]].xd[0];
		nd[j+j+j+1] = no[n[j]].xd[1];
		nd[j+j+j+2] = no[n[j]].xd[2];
	}
	double strain[6] = {};

	for (int j = 0; j < 6; j++){
		for (int k = 0; k < 24; k++){
			strain[j] += B[j][k] * nd[k];
		}
	}

	for (int j = 0; j < 6; j++){
		el[i].strain[j] = strain[j];
	}
}

void calc_stress(double DB[][24], element *el, node *no, int i){

	double strain[6] = {};
	int n[8] = {};					//全体節点番号

	for (int j = 0; j < 8; j++){
		n[j] = el[i].node[j];
	}

	double nd[24] = {};			//一要素の各節点の座標値

	for (int j = 0; j < 8; j++){
		nd[j + j + j] = no[n[j]].xd[0];
		nd[j + j + j + 1] = no[n[j]].xd[1];
		nd[j + j + j + 2] = no[n[j]].xd[2];
	}
	/******************debug*******************/
	//for (int j = 0; j < 24; j++){
	//	printf("%1.10f\n", nd[j]);
	//}
	/******************debug*******************/
	double stress[6] = {};

	/******************debug*******************/
	for (int j = 0; j < 6; j++){
		for (int k = 0; k < 24; k++){
				stress[j] += DB[j][k]*nd[k];
				//printf("%lf,",DB[j][k]);
		}
		//printf("\n");
	}
	/******************debug*******************/

	for (int j = 0; j < 6; j++){
		el[i].stress[j] = stress[j];
	}


}

void calc_strain1(double B[][24], element *el, node *no, int i){


	int n[8] = {};					//全体節点番号

	for (int j = 0; j < 8; j++){
		n[j] = el[i].node[j];
	}

	double nd[24] = {};			//一要素の各節点の座標値

	for (int j = 0; j < 8; j++){
		nd[j + j + j] = no[n[j]].xd[0];
		nd[j + j + j + 1] = no[n[j]].xd[1];
		nd[j + j + j + 2] = no[n[j]].xd[2];
	}

	double strain[6] = {};

	for (int j = 0; j < 6; j++){
		for (int k = 0; k < 24; k++){
			strain[j] += B[j][k] * nd[k];
		}
	}

	for (int j = 0; j < 6; j++){
		el[i].strain[j] = strain[j];
	}

}