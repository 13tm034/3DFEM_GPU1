#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"Header.h"

#define WORD 16384

void remove_matrix(int ElementToBeRemoved[][100], int num, double **K, element *el, node *no,material *m){

	int lim;
	lim = ElementToBeRemoved[num][0];
	int point = 1;
	while( point != lim+1){
		int ETR = ElementToBeRemoved[num][point];
		printf("ETR;\t%d", ETR);
		el[ETR].rm = true;
		
		double RKe[24][24] = {};
		/*テキストから要素剛性行列を読み込み*/
		FILE *fp_element;
		errno_t errors;
		char file_name[] = "ElementStiffnes";
		char readline[WORD] = { '\0' };
		if (errors = fopen_s(&fp_element, file_name, "r") != 0){
			printf("\n file open failed \n");
		}
		int count = 0;
		while (fgets(readline, WORD, fp_element) != NULL){
			if (ETR==count){
				convert(readline, RKe);
			}
			count++;
		}

		/*****************************debug*********************************/
		FILE *fp_rmelement;
		errno_t errors_rm;
		char file_name_rm[] = "RemovedElement";
		if (errors_rm = fopen_s(&fp_rmelement, file_name_rm, "w") != 0){
			printf("\n file open failed \n");
		}
		for (int i = 0; i < 24; i++){
			for (int j = 0; j < 24; j++){
				fprintf(fp_rmelement, "%lf,", RKe[i][j]);
			}
		}
		for (int j = 0; j < 8; j++){
			fprintf(fp_rmelement, "%d\n", el[ETR].node[j]);
		}
		fclose(fp_rmelement);
		/*****************************debug*********************************/

		/*行列を剛性行列から削除*/
		removeK(el, RKe, K, ETR,no);
		point++;
	}

}

void RemoveMatrixCSR(int ElementToBeRemoved[][100], int num, double *CSR_Kval, int *CSR_col, int *CSR_row, element *el, node *no, material *m){
	
	int lim = ElementToBeRemoved[num][0]+1;
	for (int Index = 1; Index < lim; Index++){
		int ENTR = ElementToBeRemoved[num][Index];
		el[ENTR].rm = true;
		double RKe[24][24] = {};
		Kematrix(no, el, m, RKe, ENTR);			//Keマトリクスの作成
		int n[8] = {};
		for (int i = 0; i < 8; i++){
			n[i] = el[ENTR].node[i];
			printf("node%d\n", n[i]);
			no[n[i]].point--;
		}
		int count = 0;
		for (int Row = 0; Row < 8; Row++){
		for (int Col = 0; Col < 8; Col++){
			
			for (int RowXYZ = 0; RowXYZ < 3; RowXYZ++){
			for (int ColXYZ = 0; ColXYZ < 3; ColXYZ++){
						
					int TargetRow = 3 * n[Row] + RowXYZ;
					int TargetCol = 3 * n[Col] + ColXYZ;
					double DividingVal = RKe[3 * Row + RowXYZ][3 * Col + ColXYZ];
					printf("%d,%d,%lf", TargetCol,TargetRow,DividingVal);
					int check = 0;
					for (int ValIndex = CSR_row[TargetRow]; ValIndex < CSR_row[TargetRow + 1]; ValIndex++){
						if (CSR_col[ValIndex] == TargetCol){
							CSR_Kval[ValIndex] = CSR_Kval[ValIndex] - DividingVal;
							printf("\t%d,%d,%lf\n",TargetRow, TargetCol,CSR_Kval[ValIndex]);
							count++;
							check = 1;
						}
					}
					if (check == 0)printf("\n");
					
			}				
			}
		}
		}
		printf("count=%d\n", count);

	}
	system("PAUSE");
}


void removeK(element *el, double Ke[][24], double **K, int i,node *no){


	/*要素接点番号から全体節点番号を出力*/
	int n[8] = {};					//全体節点番号

	for (int j = 0; j < 8; j++){
		n[j] = el[i].node[j];
		printf("%d\n", n[j]);
		no[n[j]].point--;
		printf("no.point\t%d", no[n[j]].point);
	}
	system("PAUSE");




	/*全体節点番号に従い要素剛性行列を全体剛性行列から除去*/
	for (int j = 0; j < 8; j++){
		for (int k = 0; k < 8; k++){
			//printf("%d,%d\n", n[j] + n[j] + n[j], i);
			K[n[j] + n[j] + n[j]][n[k] + n[k] + n[k]] -= Ke[j + j + j][k + k + k];

			K[n[j] + n[j] + n[j] + 1][n[k] + n[k] + n[k]] -= Ke[j + j + j + 1][k + k + k];
			K[n[j] + n[j] + n[j] + 2][n[k] + n[k] + n[k]] -= Ke[j + j + j + 2][k + k + k];

			K[n[j] + n[j] + n[j]][n[k] + n[k] + n[k] + 1] -= Ke[j + j + j][k + k + k + 1];
			K[n[j] + n[j] + n[j]][n[k] + n[k] + n[k] + 2] -= Ke[j + j + j][k + k + k + 2];

			K[n[j] + n[j] + n[j] + 1][n[k] + n[k] + n[k] + 1] -= Ke[j + j + j + 1][k + k + k + 1];
			K[n[j] + n[j] + n[j] + 2][n[k] + n[k] + n[k] + 2] -= Ke[j + j + j + 2][k + k + k + 2];

			K[n[j] + n[j] + n[j] + 1][n[k] + n[k] + n[k] + 2] -= Ke[j + j + j + 1][k + k + k + 2];
			K[n[j] + n[j] + n[j] + 2][n[k] + n[k] + n[k] + 1] -= Ke[j + j + j + 2][k + k + k + 1];

		}
	}


}

void convert(char *readline, double Ke[][24]){
	char *ctx;
	int count = 0;
	Ke[0][0] = atof(strtok_s(readline, ",",&ctx));
	count++;
	for (int i = 1; i < 24; i++){
		Ke[0][i] = atof(strtok_s(NULL, ",", &ctx));
		count++;
	}
	for (int i = 1; i < 24; i++){
		for (int j = 0; j < 24; j++){
			Ke[i][j] = atof(strtok_s(NULL, ",", &ctx));
			count++;
			printf("\n%d\n",count);
		}
	}

}