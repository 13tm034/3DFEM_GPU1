#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"Header.h"

void RemoveNodesCSR(double *BC_CSR_Kval, int *BC_CSR_col, int *BC_CSR_row,double *S, node *no,int N,int BC_RealNum){

	for (int num = 0; num < N; num++){
		if (no[num].point == 0){
			for (int axis = 0; axis < 3; axis++){

				int DenseIndex = num * 3 + axis;
				S[DenseIndex] = 0;
				for (int TargetVal = 0; TargetVal < BC_RealNum; TargetVal++){
					if (BC_CSR_col[TargetVal] == DenseIndex)BC_CSR_Kval[TargetVal] = 0;
				}
				for (int TargetVal = BC_CSR_row[DenseIndex]; TargetVal < BC_CSR_row[DenseIndex]; TargetVal++){
					if (BC_CSR_col[TargetVal] == DenseIndex)BC_CSR_Kval[TargetVal] = 1;
					else BC_CSR_Kval[TargetVal] = 0;
				}

			}

		}
	}



}
