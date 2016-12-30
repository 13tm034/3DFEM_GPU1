// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include<math.h>

/* Using updated (v2) interfaces to cublas */
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>

#include <helper_functions.h>  // helper for shared functions common to CUDA Samples
#include <helper_cuda.h>       // helper function CUDA error checking and initialization
#include"Header.h"

void solve_matrix(double **A, double *b, int N, node *no,int *PreCount){

	double error[1000000] = {};

	long double sum = 0;
	double *p0, *p1;
	p0 = (double *)malloc(3 * N*sizeof(double));
	p1 = (double *)malloc(3 * N*sizeof(double));



	double *x0, *x1;
	x0 = (double *)malloc(3 * N*sizeof(double));
	x1 = (double *)malloc(3 * N*sizeof(double));

	double *r0, *r1;
	r0 = (double *)malloc(3 * N*sizeof(double));
	r1 = (double *)malloc(3 * N*sizeof(double));

	double num = 0, *denov, deno = 0;
	denov = (double *)malloc(3 * N*sizeof(double));

	double *numv;
	numv = (double *)malloc(3 * N*sizeof(double));

	double alpha = 0;
	double beta = 0;

	double *Ax0;
	Ax0 = (double *)malloc(3 * N*sizeof(double));

	double *TempDiagonal;
	TempDiagonal = (double *)malloc(3 * N*sizeof(double));


	for (int i = 0; i < N + N + N; i++){
		p0[i] = 0;
		p1[i] = 0;
		x0[i] = 0;
		x1[i] = 0;
		r0[i] = 0;
		r1[i] = 0;
		denov[i] = 0;
		numv[i] = 0;
		Ax0[i] = 0;
		TempDiagonal[i] = 0;
	}


	for (int i = 0; i < N + N + N; i++){
		TempDiagonal[i] = A[i][i];
	}


	for (int i = 0; i < N + N + N; i++){
		for (int j = 0; j < N + N + N; j++){
			A[i][j] = A[i][j] / sqrt(TempDiagonal[i]);
			A[i][j] = A[i][j] / sqrt(TempDiagonal[j]);
		}
		b[i] = b[i] / sqrt(TempDiagonal[i]);
	}




	for (int i = 0; i<N+N+N; i++){
		for (int j = 0; j < N + N+N; j++){
			Ax0[i] += A[i][j] * x0[j];
		}
		r0[i] = b[i] - Ax0[i];
	}



	for (int i = 0; i<N + N + N; i++){
		p0[i] = r0[i];
	}




	for (int i = 0; i < N + N + N; i++){
		sum += r0[i] * r0[i];
	}
	int count = 1;
	error[0] = 1;


	printf("Trial:\t%d\tResidual Error:\t%1.20f\n", 0, sum);
	clock_t start, end;
	start = clock();
	int check = 0;
	while (1){
		/*ˆÀ’è«‚Ì”»’f‚ðŠÜ‚ß‚½Žc·Žû‘©”»’è*/
		//if (sum < EPS*EPS){
		//	double  stability = fabs(error[count] - error[count - 100]);
		//	if (stability < EPS / 10){
		//		break;		//EPS‚Ì”»’f
		//	}
		//}
	
		check = 0;
		num = 0;
		deno = 0;
		
		for (int i = 0; i < N + N + N; i++){
			denov[i] = 0;
		}

		for (int i = 0; i < N + N + N; i++){
			num += p0[i] * r0[i];					//a‚Ì•ªŽq
			for (int j = 0; j < N + N + N; j++){
				denov[i] += p0[j] * A[j][i];		//a‚Ì•ª•ê	
			}
		}
		
		for (int i = 0; i < N + N + N; i++){
			deno += denov[i] * p0[i];				//a‚Ì•ª•ê
		}
		alpha = num / deno;							//a‚ÌŒvŽZ

		for (int i = 0; i < N + N + N; i++){				//xk+1‚ÌŒvŽZ
			x1[i] = x0[i] + alpha*p0[i];
		}

		for (int i = 0; i < N + N + N; i++){				//rk+1‚ÌŒvŽZ
			double aAp = 0;
			for (int j = 0; j < N + N + N; j++){
				aAp += alpha*A[i][j] * p0[j];
			}
			r1[i] = r0[i] - aAp;
			/*printf("r1[%d]=%lf\n", i, r1[i]);*/
		}

		sum = 0;
		for (int i = 0; i < N + N + N; i++){
			sum += r1[i] * r1[i];
		}

		error[count] = sum;
		printf("Trial:\t%d\tResidual Error:\t%1.20f\n", count, sum);
		/*ˆÀ’è«‚Ì”»’f‚ðŠÜ‚ß‚½Žc·Žû‘©”»’è*/
		//if (sum < EPS*EPS){
		//	double  stability= fabs(error[count] - error[count - 100]);
		//	if (stability < EPS*EPS / 10){
		//		break;		//EPS‚Ì”»’f
		//	}
		//}
		if (sum <= 0)break;

		check = 1;

		for (int i = 0; i < N + N + N; i++){
			denov[i] = 0;
		}

		num = 0;
		for (int i = 0; i < N + N + N; i++){
			numv[i] = 0;
		}
		for (int i = 0; i < N + N + N; i++){
			for (int j = 0; j < N + N + N; j++){
				numv[i] += r1[j] * A[i][j];
			}
			num += numv[i] * p0[i];
		}

		deno = 0;
		for (int i = 0; i < N + N + N; i++){
			for (int j = 0; j < N + N + N; j++){
				denov[i] += p0[j] * A[i][j];
			}
			deno += denov[i] * p0[i];
		}

		beta = num / deno;
		//printf("beta=%lf\n", beta);
		//printf("num=%lf\n", num);
		//printf("deno=%lf\n", deno);

		for (int i = 0; i < N + N + N; i++){
			p0[i] = r1[i] - beta*p0[i];
			x0[i] = x1[i];
			r0[i] = r1[i];
		}

		//printf("\n\n");


		count++;
		//if(count == *PreCount)break;
	}
	end = clock();
	printf("ˆ—ŽžŠÔ-CG@%d\n", end - start);
	if (check == 1){
		for (int i = 0; i < N + N + N; i++){
			b[i] = x0[i] / sqrt(TempDiagonal[i]);
		}
	}
	else{
		for (int i = 0; i < N + N + N; i++){
			b[i] = x1[i] / sqrt(TempDiagonal[i]);
		}
	}



	printf("SOLVING PROCESS END");
	*PreCount = count;
	free(p0);


}