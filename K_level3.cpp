#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"Header.h"

void Bmatrix(node *no, element *el, double B[][24], int i,double *r,double J[][3]){


	/*要素接点から座標を出力*/
	int n[8] = {};					//全体節点番号
	


	for (int j = 0; j < 8; j++){
		n[j] = el[i].node[j];
	}


	double nd[3][8] = {};			//一要素の各節点の座標値

	for (int j = 0; j < 8; j++){
		nd[0][j] = no[n[j]].x[0];
		nd[1][j] = no[n[j]].x[1];
		nd[2][j] = no[n[j]].x[2];
		//printf("node %d = (%lf,%lf,%lf)\n", j, nd[0][j], nd[1][j], nd[2][j]);
	}

	

	/*jacobマトリクスの作成*/
	double dndr[3][8] = {};		//dn/dr
	double kariJ[3][3] = {};
	jacob(J, dndr, nd, r);


	for (int j = 0; j < 3; j++){
		for (int k = 0; k < 3; k++){
			kariJ[j][k] = J[j][k];
		}
	}



	/*jacobマトリクスの逆行列*/
	double invJ[3][3] = {};		//jacobの逆行列
	inv_jacob(invJ, J);
	for (int j = 0; j < 3; j++){
		for (int k = 0; k < 3; k++){
			J[j][k] = kariJ[j][k];
		}
	}


	/*dndxの計算*/
	double dndx[3][8] = {};		//dn/dx
	DnDx(dndx, invJ,dndr);
	

	/*Bマトリクスへ展開*/
	for (int j = 0; j < 8; j++){
		B[0][j + j + j] = dndx[0][j];
		B[1][j + j + j + 1] = dndx[1][j];
		B[2][j + j + j + 2] = dndx[2][j];
		
		B[3][j + j + j] = dndx[1][j];
		B[3][j + j + j + 1] = dndx[0][j];
		
		B[5][j + j + j] = dndx[2][j];
		B[5][j + j + j + 2] = dndx[0][j];
		
		B[4][j + j + j + 1] = dndx[2][j];
		B[4][j + j + j + 2] = dndx[1][j];
	}


}

void Dmatrix(double D[][6], material *m, int i, element *el){	//20160609debug

	double e = m[el[i].m].e;		// ポアソン比
	double v = m[el[i].m].v;		//縦弾性係数

	double a = e*(1 - v) / ((1 + v)*(1 - 2 * v));


	for (int j = 0; j < 3; j++){			//対角成分
		for (int k = 0; k < 3; k++){
			D[j][k] = 1;
		if (j == k){
			D[j + 3][k + 3] = (1 - 2 * v) / (2 * (1 - v));//if
		}
		}
	}

	D[0][1] = v / (1 - v);					//対角成分以外
	D[0][2] = v / (1 - v);
	D[1][2] = v / (1 - v);
	D[1][0] = v / (1 - v);
	D[2][0] = v / (1 - v);
	D[2][1] = v / (1 - v);


	for (int j = 0; j < 6; j++){
		for (int k = 0; k < 6; k++){
			D[j][k] = a*D[j][k];
		}
	}

}


void BDBmatrix(double B[][24], double D[][6], double Ke[][24]){

	double Bt[24][6] = {};					//Bの転地
	rollingB(Bt, B);



	double BD[24][6] = {};					//BとDの積
	productBtD(Bt, D, BD);


	productBDB(BD, B, Ke);					//BDとBの積


}

void set_eRC(double **K, node *no, double *S, double *temp, int i, int j,int N){

	int count = 1;
	if (no[i].xrc[j] == 1){									//変位拘束が入っていれば

		printf("Restraint Node:\t%d\n", i);
		for (int k = 0; k < N + N + N; k++){
			temp[k] = K[k][i + i + i + j] * no[i].xd[j];	//tempに節点変位と剛性行列の3i+j行の一列を変えたものを代入
		}
		temp[i + i + i + j] = 0;							//変位拘束がかかった節点の拘束方向に0を代入
		S[i + i + i + j] = no[i].xd[j];						//変位拘束のかかった節点の変位方向の変位を入力情報にする。


		for (int k = 0; k < N + N + N; k++){
			S[k] = S[k] - temp[k];							//初期の応力ベクトルからtempを引く
			K[k][i + i + i + j] = 0;						//3j+i行目の一列、3j+i列目の一行に0を代入
			K[i + i + i + j][k] = 0;
		}

		K[i + i + i + j][i + i + i + j] = 1;				//3j+i行3j+i列目に1を代入
		count++;
	}

}

void det_jacob(double J[][3], double *detJ){


	*detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
		- J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
		+ J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
	
	/*printf("detJ=%lf", *detJ);

	if (*detJ <= 0){
		system("PAUSE");


	}*/

}

void product_detJ_weight(double pKe[][24], double detJ, double *w){

	for (int i = 0; i < 24; i++){
		for (int j = 0; j < 24; j++){

			pKe[i][j] = pKe[i][j] * detJ*w[0] * w[1] * w[2];

		}
	}


}

void addKe(double pKe[][24], double Ke[][24]){

	for (int i = 0; i < 24; i++){
		for (int j = 0; j < 24; j++){
			Ke[i][j] = Ke[i][j] + pKe[i][j];
		}
	}


}


void QsortPreCOO(double *COO_Kval, int *COO_col, int *COO_row, int left, int right){
	int pivot_row, pivot_col, l_hold, r_hold;
	double pivot_Kval;
	l_hold = left;
	r_hold = right;
	pivot_row = COO_row[left];
	pivot_col = COO_col[left];
	pivot_Kval = COO_Kval[left];

	while (left < right)
	{
		while ((COO_row[right] >= pivot_row) && (left < right))right--;
		if (left != right)
		{
			COO_row[left] = COO_row[right];
			COO_col[left] = COO_col[right];
			COO_Kval[left] = COO_Kval[right];
			left++;
		}
		while ((COO_row[left] <= pivot_row) && (left < right))left++;
		if (left != right)
		{
			COO_row[right] = COO_row[left];
			COO_col[right] = COO_col[left];
			COO_Kval[right] = COO_Kval[left];
			right--;
		}
	}
	COO_row[left] = pivot_row;
	COO_col[left] = pivot_col;
	COO_Kval[left] = pivot_Kval;
	pivot_row = left;
	left = l_hold;
	right = r_hold;
	if (left < pivot_row)QsortPreCOO(COO_Kval, COO_col, COO_row, left, pivot_row - 1);
	if (right > pivot_row)QsortPreCOO(COO_Kval, COO_col, COO_row, pivot_row + 1, right);
}

void set_eRC_CSR(node *no,double *S,double *TempVector,double *CSR_Kval, int *CSR_col, int *CSR_row, int NodeNumber, int XYZ, int N){

	if (no[NodeNumber].xrc[XYZ] == 1){
		printf("Restrain Node:\t%d\n", NodeNumber);
		for (int i = 0; i < N + N + N; i++){
			for (int j = CSR_row[i]; j < CSR_row[i + 1]; j++){
				if (CSR_col[j] == 3 * NodeNumber + XYZ){
					TempVector[i] = CSR_Kval[j] * no[NodeNumber].xd[XYZ];
					break;
				}
			}
		}
		TempVector[3 * NodeNumber + XYZ] = 0;
		S[3 * NodeNumber + XYZ] = no[NodeNumber].xd[XYZ];

		for (int i = 0; i < N + N + N; i++){
			S[i] = S[i] - TempVector[i];
			for (int j = CSR_row[i]; j < CSR_row[i + 1]; j++){
				if (CSR_col[j] == 3 * NodeNumber + XYZ){
					//printf("j=%d\t", j);
					//printf("CSR_col=%d\n", CSR_col[j]);
					CSR_Kval[j] = 0;
				}
			}
			for (int j = CSR_row[3 * NodeNumber + XYZ]; j < CSR_row[3 * NodeNumber + XYZ + 1]; j++){
				CSR_Kval[j] = 0;
			}
		}

		for (int i = CSR_row[3 * NodeNumber + XYZ]; i < CSR_row[3 * NodeNumber + XYZ+1]; i++){
			if (CSR_col[i] == 3 * NodeNumber + XYZ){
				CSR_Kval[i] = 1;
				break;
			}
		}
	}
	//int count = 1;
	//if (no[i].xrc[j] == 1){									//変位拘束が入っていれば

	//	printf("Restraint Node:\t%d\n", i);
	//	for (int k = 0; k < N + N + N; k++){
	//		temp[k] = K[k][i + i + i + j] * no[i].xd[j];	//tempに節点変位と剛性行列の3i+j行の一列を変えたものを代入
	//	}
	//	temp[i + i + i + j] = 0;							//変位拘束がかかった節点の拘束方向に0を代入
	//	S[i + i + i + j] = no[i].xd[j];						//変位拘束のかかった節点の変位方向の変位を入力情報にする。


	//	for (int k = 0; k < N + N + N; k++){
	//		S[k] = S[k] - temp[k];							//初期の応力ベクトルからtempを引く
	//		K[k][i + i + i + j] = 0;						//3j+i行目の一列、3j+i列目の一行に0を代入
	//		K[i + i + i + j][k] = 0;
	//	}

	//	K[i + i + i + j][i + i + i + j] = 1;				//3j+i行3j+i列目に1を代入
	//	count++;
	//}

}
