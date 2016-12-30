#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"Header.h"

void Kmatrix(node *no, element *el, material *m, double *preCOO_Kval, int *preCOO_col, int *preCOO_row, int E, int N, int M, int *RealNumberOfValues){

	FILE *fp_element;
	errno_t errors;
	char file_name[] = "ElementStiffnes1";
	if (errors = fopen_s(&fp_element, file_name, "w") != 0){
		printf("\n file open failed1 \n");
	}
	fclose(fp_element);
	
	int MaximumNumberOfValues = 24 * 24 * E;	//値がある要素の理論上最大値	

	for (int i = 0; i < E; i++){			//要素iについて
		double Ke[24][24] = {};
		//Kematrix_log(i);
		Kematrix(no, el, m, Ke, i);			//Keマトリクスの作成
		//Kematrix_log1(Ke, i);
		AddKe2preCOO(el, Ke, preCOO_Kval, preCOO_col, preCOO_row, MaximumNumberOfValues, i);
		//addK(el, Ke,K,i,fp_element);
	}


	//並び変えて0要素を圧縮
	printf("SORTING PROCESS\n");
	Sort_preCOO(preCOO_Kval, preCOO_col, preCOO_row, MaximumNumberOfValues);
	preCOO2COO(preCOO_Kval, preCOO_col, preCOO_row, MaximumNumberOfValues,RealNumberOfValues);


	for (int i = 0; i < *RealNumberOfValues; i++){
		int StartIndex;
		int EndIndex;
		if (i == 0){
			StartIndex = 0;
		}
		else if (preCOO_row[i - 1] != preCOO_row[i]){
			StartIndex = i;
		}
		if (preCOO_row[i] != preCOO_row[i+1]){
			EndIndex = i;
			Sort_preCOO_col(preCOO_Kval, preCOO_col, preCOO_row, StartIndex, EndIndex);
		}
	}

	//for (int i = 900; i < 1000; i++){
	//	printf("%d,%d\t%lf\t%lf\n", preCOO_row[i], preCOO_col[i], preCOO_Kval[i], K[preCOO_row[i]][preCOO_col[i]]);
	//}


	printf("RealNumberOfValues=%d\n", *RealNumberOfValues);
	system("PAUSE");
}

void setBC(node *no, double *S, double **K, int E, int N, int M){
	load_stress(no, S,N);

	set_RC(no, S, K, E, N, M);

}


void COO2CSR(int *preCSR_row, int *COO_row,int *CountRow,int RealNumberOfValues){

	int count = 0;
	for (int i = 0; i < RealNumberOfValues; i++){
		if (i == 0){	
			preCSR_row[count] = 0;
			count++;
		}
		else if (i > 0){
			if (COO_row[i - 1] != COO_row[i]){
				preCSR_row[count] = i;
				count++;
			}
		}
	}
	preCSR_row[count] = RealNumberOfValues;
	*CountRow = count;
}
void setBC_CSR(node *no, double *S, double *CSR_Kval, int *CSR_col, int *CSR_row,int E,int N,int M){
	load_stress(no, S, N);
	set_RC_CSR(no, S, CSR_Kval, CSR_col, CSR_row, E, N, M);
}

void CSR2BC_CSR(double *CSR_Kval, int *CSR_col, int *CSR_row, double *BC_CSR_Kval, int *BC_CSR_col, int *BC_CSR_row, int N){
	int count = 0;
	for (int i = 0; i < N + N + N; i++){
		BC_CSR_row[i] = count;
		for (int j = CSR_row[i]; j < CSR_row[i + 1]; j++){
			if (CSR_Kval[j] != 0){
				BC_CSR_Kval[count] = CSR_Kval[j];
				BC_CSR_col[count] = CSR_col[j];
				count++;
			}
		}
	}

	BC_CSR_row[N + N + N] = count;
}
