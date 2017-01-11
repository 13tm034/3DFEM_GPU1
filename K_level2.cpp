#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"Header.h"

void addK(element *el, double Ke[][24], double **K, int ElementNumber,FILE *fp){

	/*要素接点番号から全体節点番号を出力*/
	int n[8] = {};					//全体節点番号

	for (int j = 0; j < 8; j++){
		n[j] = el[ElementNumber].node[j];
	}


	//n[0] = el[i].node[4];
	//n[1] = el[i].node[7];
	//n[2] = el[i].node[6];
	//n[3] = el[i].node[5];
	//n[4] = el[i].node[0];
	//n[5] = el[i].node[3];
	//n[6] = el[i].node[2];
	//n[7] = el[i].node[1];
	

	for (int i = 0; i < 24; i++){
		for (int j = 0; j < 24; j++){
			fprintf_s(fp, "%lf,", Ke[i][j]);
		}
	}
	fprintf_s(fp, "\n");



	if (ElementNumber == 314){
		/*****************************debug*********************************/
		FILE *fp_rmelement;
		errno_t errors_rm;
		char file_name_rm[] = "RemoveElement";
		if (errors_rm = fopen_s(&fp_rmelement, file_name_rm, "w") != 0){
			printf("\n file open failed \n");
		}
		for (int i = 0; i < 24; i++){
			for (int j = 0; j < 24; j++){
				fprintf(fp_rmelement, "%lf,", Ke[i][j]);
			}
		}
		for (int j = 0; j < 8; j++){
			fprintf(fp_rmelement, "%d\n", n[j]);
		}
		fclose(fp_rmelement);
		/*****************************debug*********************************/
	}



	/*全体節点番号に従い要素剛性行列を全体剛性行列に加算*/
	for (int j = 0; j < 8; j++){
		for (int k = 0; k < 8; k++){
			//printf("%d,%d\n", n[j] + n[j] + n[j], i);
			K[n[j] + n[j] + n[j]][n[k] + n[k] + n[k]] += Ke[j + j + j][k + k + k];
			
			K[n[j] + n[j] + n[j] + 1][n[k] + n[k] + n[k]] += Ke[j + j + j + 1][k + k + k];
			K[n[j] + n[j] + n[j] + 2][n[k] + n[k] + n[k]] += Ke[j + j + j + 2][k + k + k];
			
			K[n[j] + n[j] + n[j]][n[k] + n[k] + n[k] + 1] += Ke[j + j + j][k + k + k + 1];
			K[n[j] + n[j] + n[j]][n[k] + n[k] + n[k] + 2] += Ke[j + j + j][k + k + k + 2];
			
			K[n[j] + n[j] + n[j] + 1][n[k] + n[k] + n[k] + 1] += Ke[j + j + j + 1][k + k + k + 1];
			K[n[j] + n[j] + n[j] + 2][n[k] + n[k] + n[k] + 2] += Ke[j + j + j + 2][k + k + k + 2];
			
			K[n[j] + n[j] + n[j] + 1][n[k] + n[k] + n[k] + 2] += Ke[j + j + j + 1][k + k + k + 2];
			K[n[j] + n[j] + n[j] + 2][n[k] + n[k] + n[k] + 1] += Ke[j + j + j + 2][k + k + k + 1];

		}
	}


}

void AddKe2preCOO(element *el, double Ke[][24], double *COO_Kval, int *COO_col, int *COO_row, int MaximumNumberOfValues, int ElementNumber){
	
	/*要素接点番号から全体節点番号を出力*/
	int n[8] = {};					//全体節点番号

	for (int j = 0; j < 8; j++){
		n[j] = el[ElementNumber].node[j];
	}

	int Index = ElementNumber*24*24;

	for (int ElementNode_row = 0; ElementNode_row < 8; ElementNode_row++){
		for (int ElementNode_col = 0; ElementNode_col < 8; ElementNode_col++){
			/****************************************************************************/
			//bool ThereIsInTheElements = false;
			for (int RowXYZ = 0; RowXYZ < 3; RowXYZ++){
				for (int ColXYZ = 0; ColXYZ < 3; ColXYZ++){

							COO_Kval[Index] = (double)Ke[ElementNode_row * 3 + RowXYZ][ElementNode_col * 3 + ColXYZ];
							COO_col[Index] = n[ElementNode_col] * 3 + ColXYZ;
							COO_row[Index] = n[ElementNode_row] * 3 + RowXYZ;
							if (ElementNumber==350) printf("%d,%d,%lf\n", n[ElementNode_row] * 3 + RowXYZ, n[ElementNode_col] * 3 + ColXYZ,COO_Kval[Index]);
							Index++;
				}
			}
			/****************************************************************************/	
		}
	}
	
}

void Sort_preCOO(double *COO_Kval, int *COO_col, int *COO_row, int MaximumNumberOfValues){
	QsortPreCOO(COO_Kval, COO_col,COO_row, 0, MaximumNumberOfValues - 1);
}





void preCOO2COO(double *COO_Kval, int *COO_col, int *COO_row, int MaximumNumberOfValues, int *RealNumberOfValues){
	
	int NonZeroCount = 0;
	for (int i = 0; i < MaximumNumberOfValues; i++){		
		if (COO_Kval[i] != 0){
			int j = i+1;
			while (COO_row[i]==COO_row[j]){
				if (COO_col[i] == COO_col[j] ){
					COO_Kval[i] = COO_Kval[i]+ COO_Kval[j];
					COO_Kval[j] = 0;
					COO_col[j] = 0;
				}
				j++;
			}
		}
	}
	printf("-------------------------\n");
	for (int i = 0; i < 100; i++){
		printf("%d\t%d\t%lf\n", COO_col[i], COO_row[i], COO_Kval[i]);
	}
	for (int i = 0; i < MaximumNumberOfValues; i++){
		if (COO_Kval[i] != 0&& COO_col[i]!=0)NonZeroCount++;
	}
	printf("NonZeroCount=%d\n", NonZeroCount);

	int temp = NonZeroCount;
	double *tempCOO_Kval;
	tempCOO_Kval = (double *)malloc(sizeof(double)*(NonZeroCount));
	int *tempCOO_col;
	tempCOO_col = (int *)malloc(sizeof(int)*(NonZeroCount));
	int *tempCOO_row;
	tempCOO_row = (int *)malloc(sizeof(int)*(NonZeroCount));
	for (int i = 0; i < NonZeroCount; i++){
		tempCOO_Kval[i] = 0;
		tempCOO_col[i] = 0;
		tempCOO_row[i] = 0;
	}


	printf("REMOVE ZERO ELEMENTS\n");
	int NonZeroIndex = 0;
	for (int i = *RealNumberOfValues; i < MaximumNumberOfValues; i++){
		if (COO_Kval[i] != 0 && COO_col[i]!=0){
			tempCOO_Kval[NonZeroIndex] = COO_Kval[i];
			tempCOO_col[NonZeroIndex] = COO_col[i];
			tempCOO_row[NonZeroIndex] = COO_row[i];
			NonZeroIndex++;
		}
	}

	for (int i = 0; i < MaximumNumberOfValues; i++){
		COO_Kval[i] = 0;
		COO_col[i] = 0;
		COO_row[i] = 0;
	}
	for (int i = 0; i < NonZeroIndex; i++){
		COO_Kval[i] = tempCOO_Kval[i];
		COO_col[i] = tempCOO_col[i];
		COO_row[i] = tempCOO_row[i];
	}

	*RealNumberOfValues = NonZeroIndex;
}



void Sort_preCOO_col(double *preCOO_Kval, int *preCOO_col, int *preCOO_row, int StartIndex,int EndIndex){
	int SameRowIndex = EndIndex - StartIndex + 1;

	double *Kval;
	Kval = (double *)malloc(sizeof(double)*(SameRowIndex));
	int *Col;
	Col = (int *)malloc(sizeof(int)*(SameRowIndex));
	int *Row;
	Row = (int *)malloc(sizeof(int)*(SameRowIndex));


	for (int i = 0; i < SameRowIndex; i++){
		Kval[i] = preCOO_Kval[i + StartIndex];
		Col[i] = preCOO_col[i + StartIndex];
		Row[i] = preCOO_row[i + StartIndex];
	}

	QsortPreCOO(Kval, Row, Col, 0, SameRowIndex - 1);

	for (int i = 0; i < SameRowIndex; i++){
		preCOO_Kval[i + StartIndex]=Kval[i];
		preCOO_col[i + StartIndex]=Col[i] ;
		preCOO_row[i + StartIndex] = Row[i];
	}

}

void Kematrix(node *no, element *el, material *m, double Ke[][24],int i){

	//order5 point5
	/*double a = sqrt(1120.0);
	double b = (70 + a) / 126;
	double c = (70 - a) / 126;
	double d = 1 / (15 * (c - b));*/

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

				double pKe[24][24] = {};


				double B[6][24] = {};			//Bマトリクス	
				double J[3][3] = {};			//jacobマトリクス
				Bmatrix(no, el, B, i,r,J);



				double detJ = 0;				//jacobのdetの値
				double *det;
				det = &detJ;
				det_jacob(J, det);
				double D[6][6] = {};			//Dマトリクス
				Dmatrix(D, m, i, el);

				BDBmatrix(B, D, pKe);			//BとD行列の積

				product_detJ_weight(pKe, detJ, w);		//detJとweightをpKeにかける





				addKe(pKe, Ke);							//pKeをKeに足していく

				}
			}
		}
}


	





void load_stress(node *no, double *S,int N){

	for (int i = 0; i < N; i++){
		S[i + i + i] = no[i].xf[0];
		S[i + i + i + 1] = no[i].xf[1];
		S[i + i + i + 2] = no[i].xf[2];
	}

}

void set_RC(node *no, double *S, double **K, int E, int N, int M){

	double *temp;
	temp = (double *)malloc(sizeof(double)*(N + N + N));
	for (int i = 0; i < N + N + N; i++){
		temp[i] = 0;
	}

	for (int i = 0; i < N; i++){			//節点iの
		for (int j = 0; j < 3; j++){		//j軸方向について
			set_eRC(K, no, S, temp, i, j, N);
			
		}
	}

	//free(temp);
}
void set_RC_CSR(node *no, double *S, double *CSR_Kval, int *CSR_col, int *CSR_row,int E,int N,int M)
{	
	double *TempVector;
	TempVector = (double *)malloc(sizeof(double)*(N + N + N));
	for (int i = 0; i < N + N + N; i++){
		TempVector[i] = 0;
	}

	for (int i = 0; i < N; i++){
		for (int j = 0; j < 3; j++){
			set_eRC_CSR(no,S,TempVector,CSR_Kval, CSR_col, CSR_row, i, j, N);
		}
	}
}
