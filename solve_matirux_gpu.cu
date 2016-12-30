
#include<math.h>

/* using updated (v2) interfaces to cublas */
#include <cuda_runtime.h>
#include<cusparse.h>
#include<cusparse_v2.h>
#include<cublas_v2.h>

#include <helper_functions.h>  // helper for shared functions common to cuda samples
#include <helper_cuda.h>       // helper function cuda error checking and initialization
#include"header.h"


void solve_matrix_gpu(double **A, double *b, int N, node *no){
	
	//PreConditoned
	double *TempDiagonal;
	TempDiagonal = (double *)malloc(3 * N*sizeof(double));
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
	
	
	
	//ホスト側のメモリを確保
	double *hA,*hx,*hb,*hr,*hp,*HostTempVector;
	hA = (double *)malloc(sizeof(double)*(9 * N*N));
	hb = (double *)malloc(sizeof(double)*(N + N + N));
	hx = (double *)malloc(sizeof(double)*(N + N + N));
	hr = (double *)malloc(sizeof(double)*(N + N + N));
	hp = (double *)malloc(sizeof(double)*(N + N + N));
	HostTempVector = (double *)malloc(sizeof(double)*(N + N + N));

	for (int i = 0; i < N + N + N; i++){
		hp[i] = 0;
		hr[i] = 0;
		hp[i] = 0;
		hx[i] = 0;
		hb[i] = 0;
		HostTempVector[i] = 0;
	}


	//データをコピー
	for (int i = 0; i < N + N+ N; i++){
		for (int j = 0; j < N + N + N; j++){
			hA[i * 3 * N + j] = A[i][j];
		}
	}

	for (int i = 0; i < N; i++){
		hx[i + i + i] = no[i].xd[0];
		hx[i + i + i + 1] = no[i].xd[1];
		hx[i + i + i + 2] = no[i].xd[2];
	}

	for (int i = 0; i < N+N+N; i++){
		hb[i] = b[i];
	}
	

	cusparseHandle_t cusparseHandle = 0;
	cusparseMatDescr_t descr = 0;

	cusparseCreate(&cusparseHandle);


	//デバイス側のメモリを確保
	double *dA, *db, *dx, *dr, *dp, *TempVector;
	double CorrectionCoefficientA = 0;	//修正係数1
	double CorrectionCoefficientB = 0;	//修正係数2
	int *npr;
	cudaMalloc((void**)&dA, sizeof(double)* 9 * N*N);
	cudaMalloc((void**)&db, sizeof(double)* 3 * N);
	cudaMalloc((void**)&dx, sizeof(double)* 3 * N);
	cudaMalloc((void**)&dr, sizeof(double)* 3 * N);
	cudaMalloc((void**)&dp, sizeof(double)* 3 * N);


	cudaMalloc((void**)&npr, sizeof(int)* 9 * N*N);
	cudaMalloc((void**)&TempVector, sizeof(double)* 3 * N);
	/*Create CUSPARSE context*/
	cusparseCreateMatDescr(&descr);
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

	/*Create CUBLAS context*/
	cublasHandle_t cublasHandle = 0;
	cublasStatus_t cublasStatus;
	cublasStatus = cublasCreate(&cublasHandle);


	//初期値を転送
	cudaMemcpy(dA, hA, sizeof(double)*(9 * N*N), cudaMemcpyHostToDevice);
	cudaMemcpy(db, hb, sizeof(double)*(N + N + N), cudaMemcpyHostToDevice);
	cudaMemcpy(dx, hx, sizeof(double)*(N + N + N), cudaMemcpyHostToDevice);
	cudaMemcpy(dp, hp, sizeof(double)*(N + N + N), cudaMemcpyHostToDevice);
	cudaMemcpy(dr, hr, sizeof(double)*(N + N + N), cudaMemcpyHostToDevice);
	cudaMemcpy(TempVector, HostTempVector, sizeof(double)*(N + N + N), cudaMemcpyHostToDevice);


	int total;

	cusparseDnnz(cusparseHandle, CUSPARSE_DIRECTION_ROW, 3 * N, 3 * N, descr, dA, 3 * N, npr, &total);

	double *csrV_A;
	int *csrC_A;
	int *csrR_A;
	cudaMalloc((void**)&csrV_A, sizeof(double)*total);
	cudaMalloc((void**)&csrR_A, sizeof(int)* 3 * N + 1);
	cudaMalloc((void**)&csrC_A, sizeof(int)*total);

	cusparseDdense2csr(cusparseHandle, 3 * N, 3 * N, descr, dA, 3 * N, npr, csrV_A, csrR_A, csrC_A);
	double *h_csrV_A;
	h_csrV_A = (double *)malloc(sizeof(double)*total);
	int *h_csrC_A;
	h_csrC_A = (int *)malloc(sizeof(int)* total);
	int *h_csrR_A;
	h_csrR_A = (int *)malloc(sizeof(int)* 3 * N + 1);

	cudaMemcpy(h_csrV_A, csrV_A, sizeof(double)*total, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_csrC_A, csrC_A, sizeof(int)*total, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_csrR_A, csrR_A, sizeof(int)*(N + N + N), cudaMemcpyDeviceToHost);
	printf("total=%d\n", total);
	FILE *fp_V,*fp_C,*fp_R;
	errno_t errorV, errorC, errorR;
	char fnameV[] = "gpuV";
	char fnameC[] = "gpuC";
	char fnameR[] = "gpuR";
	if (errorV = fopen_s(&fp_V, fnameV, "w") != 0){
		printf("\n file open failed \n");
	}

	if (errorC = fopen_s(&fp_C, fnameC, "w") != 0){
		printf("\n file open failed \n");
	}
	
	if (errorR = fopen_s(&fp_R, fnameR, "w") != 0){
		printf("\n file open failed \n");
	}
	for (int i = 0; i < total; i++){
		fprintf(fp_V, "%19.15f\n", h_csrV_A[i]);
		fprintf(fp_C, "%d\n", h_csrC_A[i]);
	}
	for (int i = 0; i < N + N + N + 1; i++){
		fprintf(fp_R, "%d\n", h_csrR_A[i]);
	}
	fclose(fp_V);
	fclose(fp_C);
	fclose(fp_R);

	
	double alpha = -1.0;
	double beta = 1.0;
	cublasDcopy(cublasHandle, N + N + N, db, 1, TempVector, 1);
	cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, 3 * N, 3 * N, total, &alpha, descr, csrV_A, csrR_A, csrC_A, dx, &beta, db);
	cublasDcopy(cublasHandle, N + N + N, db, 1, dr, 1);
	cublasDcopy(cublasHandle, N + N + N, TempVector, 1, db, 1);
	cublasDcopy(cublasHandle, N + N + N, dr, 1, dp, 1);

	double ResidualError = 0;
	cublasDdot(cublasHandle, N + N + N, dr, 1, dr, 1, &ResidualError);
	printf("iteration:\t%d\tresidual error:\t%1.20f\n", 0, ResidualError);
	

	int Iteration = 1;
	double error[1000000] = {};
	error[0] = ResidualError;
	clock_t start, end;
	start = clock();

	while (1){
		double Denominator = 0.0f;	//分母
		double Numerator = 0.0f;		//分子

		cublasDdot(cublasHandle, N + N + N, dp, 1, dr, 1, &Numerator);
	
		alpha = 1.0;
		beta = 0.0f;
		cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, 3 * N, 3 * N, total, &alpha, descr, csrV_A, csrR_A, csrC_A, dp, &beta, TempVector);
		cublasDdot(cublasHandle, N + N + N, dp, 1, TempVector, 1, &Denominator);
		
		CorrectionCoefficientA = Numerator / Denominator;
		
		cublasDaxpy(cublasHandle, N + N + N, &CorrectionCoefficientA, dp, 1, dx, 1);
	
		alpha = -CorrectionCoefficientA;
		beta = 1;
		cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, 3 * N, 3 * N, total, &alpha, descr, csrV_A, csrR_A, csrC_A, dp, &beta, dr);
		
		cublasDdot(cublasHandle, N + N + N, dr, 1, dr, 1, &ResidualError);
			
		Iteration++;
		error[Iteration] = ResidualError;
		//if (ResidualError < EPS*EPS){
		//	double  stability = fabs(error[Iteration] - error[Iteration - 1000]);
		//	if (stability < EPS*EPS / 10){
		//		break;		//EPSの判断
		//	}
		//}
		if (ResidualError == 0)break;


		alpha = 1.0f;
		beta = 0.0f;
		cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, 3 * N, 3 * N, total, &alpha, descr, csrV_A, csrR_A, csrC_A, dp, &beta, TempVector);
		cublasDdot(cublasHandle, N + N + N, dr, 1, TempVector, 1, &Numerator);
		cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, 3 * N, 3 * N, total, &alpha, descr, csrV_A, csrR_A, csrC_A, dp, &beta, TempVector);
		cublasDdot(cublasHandle, N + N + N, dp, 1, TempVector, 1, &Denominator);
		CorrectionCoefficientB = -Numerator / Denominator;
		cublasDscal(cublasHandle, N + N + N, &CorrectionCoefficientB, dp, 1);
		CorrectionCoefficientB = 1;
		cublasDaxpy(cublasHandle, N + N + N,&CorrectionCoefficientB , dr, 1, dp, 1);
	}
	end = clock();
	printf("処理時間-CG　%d\n", end - start);
	cudaMemcpy(hx, dx, sizeof(double)*(N + N + N), cudaMemcpyDeviceToHost);
	for (int i = 0; i < N + N + N; i++){
		b[i] = hx[i] / sqrt(TempDiagonal[i]);
	}
}




void solve_matrix_gpu_CSR(double *CSR_Kval,int *CSR_col,int *CSR_row, double *b, int N, node *no,int RealNumberOfValues){

	//Precondtioned
	double *TempDiagonal;
	TempDiagonal = (double *)malloc(3 * N*sizeof(double));
	for (int i = 0; i < N + N + N; i++){
		for (int j = CSR_row[i]; j < CSR_row[i+1]; j++){
			if (CSR_col[j] == i)TempDiagonal[i] = CSR_Kval[j];
		}
	}





	//PreConditoned
	//double *TempDiagonal;
	//TempDiagonal = (double *)malloc(3 * N*sizeof(double));
	//for (int i = 0; i < N + N + N; i++){
	//	TempDiagonal[i] = A[i][i];
	//}
	//for (int i = 0; i < N + N + N; i++){
	//	for (int j = 0; j < N + N + N; j++){
	//		A[i][j] = A[i][j] / sqrt(TempDiagonal[i]);
	//		A[i][j] = A[i][j] / sqrt(TempDiagonal[j]);
	//	}
	//	b[i] = b[i] / sqrt(TempDiagonal[i]);
	//}
	//
	
	
	//ホスト側のメモリを確保
	double *h_CSR_Kval,*hx,*hb,*hr,*hp,*HostTempVector;
	int *h_CSR_col, *h_CSR_row;
	h_CSR_Kval = (double *)malloc(sizeof(double)*(RealNumberOfValues));
	h_CSR_col = (int *)malloc(sizeof(int)*(RealNumberOfValues));
	h_CSR_row = (int *)malloc(sizeof(int)*(N + N + N));
	hb = (double *)malloc(sizeof(double)*(N + N + N));
	for (int i = 0; i < RealNumberOfValues; i++){
		h_CSR_Kval[i] = CSR_Kval[i];
	}

	//対角スケーリング
	for (int i = 0; i < N + N + N; i++){
		for (int j = CSR_row[i]; j < CSR_row[i + 1]; j++){
			h_CSR_Kval[j] = h_CSR_Kval[j] / sqrt(TempDiagonal[i]);
			h_CSR_Kval[j] = h_CSR_Kval[j] / sqrt(TempDiagonal[CSR_col[j]]);
		}
	}

	for (int i = 0; i < RealNumberOfValues; i++){
		h_CSR_col[i] = CSR_col[i];
	}
	for (int i = 0; i < N + N + N ; i++){
		h_CSR_row[i] = CSR_row[i];
	}
	for (int i = 0; i < N + N + N; i++){
		hb[i] = b[i] / sqrt(TempDiagonal[i]);
	}

	FILE *fp_V, *fp_C, *fp_R;
	errno_t errorV, errorC, errorR;
	char fnameV[] = "cpuV";
	char fnameC[] = "cpuC";
	char fnameR[] = "cpuR";
	if (errorV = fopen_s(&fp_V, fnameV, "w") != 0){
		printf("\n file open failed \n");
	}

	if (errorC = fopen_s(&fp_C, fnameC, "w") != 0){
		printf("\n file open failed \n");
	}

	if (errorR = fopen_s(&fp_R, fnameR, "w") != 0){
		printf("\n file open failed \n");
	}
	for (int i = 0; i < RealNumberOfValues; i++){
		fprintf(fp_V, "%19.15f\n", h_CSR_Kval[i]);
		fprintf(fp_C, "%d\n", h_CSR_col[i]);
	}
	for (int i = 0; i < N + N + N; i++){
		fprintf(fp_R, "%d\n", h_CSR_row[i]);
	}
	fclose(fp_V);
	fclose(fp_C);
	fclose(fp_R);


	hx = (double *)malloc(sizeof(double)*(N + N + N));
	hr = (double *)malloc(sizeof(double)*(N + N + N));
	hp = (double *)malloc(sizeof(double)*(N + N + N));
	HostTempVector = (double *)malloc(sizeof(double)*(N + N + N));

	for (int i = 0; i < N + N + N; i++){
		hp[i] = 0;
		hr[i] = 0;
		hx[i] = 0;
		HostTempVector[i] = 0;
	}


	//データをコピー
	//for (int i = 0; i < N; i++){
	//	hx[i + i + i] = no[i].xd[0];
	//	hx[i + i + i + 1] = no[i].xd[1];
	//	hx[i + i + i + 2] = no[i].xd[2];
	//}


	

	cusparseHandle_t cusparseHandle = 0;
	cusparseMatDescr_t descr = 0;

	cusparseCreate(&cusparseHandle);


	//デバイス側のメモリを確保
	double *d_CSR_Kval,*db, *dx, *dr, *dp, *TempVector;
	int *d_CSR_col, *d_CSR_row;
	double CorrectionCoefficientA = 0;	//修正係数1
	double CorrectionCoefficientB = 0;	//修正係数2
	int *npr;

	cudaMalloc((void**)&d_CSR_Kval, sizeof(double)*RealNumberOfValues);
	cudaMalloc((void**)&d_CSR_col, sizeof(int)*RealNumberOfValues);
	cudaMalloc((void**)&d_CSR_row, sizeof(int)*(N + N + N )+1);
	cudaMalloc((void**)&db, sizeof(double)* 3 * N);
	cudaMalloc((void**)&dx, sizeof(double)* 3 * N);
	cudaMalloc((void**)&dr, sizeof(double)* 3 * N);
	cudaMalloc((void**)&dp, sizeof(double)* 3 * N);


	cudaMalloc((void**)&npr, sizeof(int)* 9 * N*N);
	cudaMalloc((void**)&TempVector, sizeof(double)* 3 * N);
	/*Create CUSPARSE context*/
	cusparseCreateMatDescr(&descr);
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

	/*Create CUBLAS context*/
	cublasHandle_t cublasHandle = 0;
	cublasStatus_t cublasStatus;
	cublasStatus = cublasCreate(&cublasHandle);


	//初期値を転送
	cudaMemcpy(d_CSR_Kval, h_CSR_Kval, sizeof(double)*RealNumberOfValues, cudaMemcpyHostToDevice);
	cudaMemcpy(d_CSR_col, h_CSR_col, sizeof(int)*RealNumberOfValues, cudaMemcpyHostToDevice);
	cudaMemcpy(d_CSR_row, h_CSR_row, sizeof(int)*(N + N + N ), cudaMemcpyHostToDevice);
	cudaMemcpy(db, hb, sizeof(double)*(N + N + N), cudaMemcpyHostToDevice);
	cudaMemcpy(dx, hx, sizeof(double)*(N + N + N), cudaMemcpyHostToDevice);
	cudaMemcpy(dp, hp, sizeof(double)*(N + N + N), cudaMemcpyHostToDevice);
	cudaMemcpy(dr, hr, sizeof(double)*(N + N + N), cudaMemcpyHostToDevice);
	cudaMemcpy(TempVector, HostTempVector, sizeof(double)*(N + N + N), cudaMemcpyHostToDevice);


	int total=RealNumberOfValues;


	//double *csrV_A;
	//int *csrC_A;
	//int *csrR_A;
	//cudaMalloc((void**)&csrV_A, sizeof(double)*total);
	//cudaMalloc((void**)&csrR_A, sizeof(int)* 3 * N + 1);
	//cudaMalloc((void**)&csrC_A, sizeof(int)*total);

	
	double alpha = -1.0;
	double beta = 1.0;
	cublasDcopy(cublasHandle, N + N + N, db, 1, TempVector, 1);
	cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, 3 * N, 3 * N, total, &alpha, descr, d_CSR_Kval, d_CSR_row, d_CSR_col, dx, &beta, db);
	cublasDcopy(cublasHandle, N + N + N, db, 1, dr, 1);
	cublasDcopy(cublasHandle, N + N + N, TempVector, 1, db, 1);
	cublasDcopy(cublasHandle, N + N + N, dr, 1, dp, 1);

	double ResidualError = 0;
	cublasDdot(cublasHandle, N + N + N, dr, 1, dr, 1, &ResidualError);
	//printf("iteration:\t%d\tresidual error:\t%1.20f\n", 0, ResidualError);
	

	int Iteration = 1;
	double error[1000000] = {};
	error[0] = ResidualError;
	clock_t start, end;
	start = clock();

	while (1){
		double Denominator = 0.0f;	//分母
		double Numerator = 0.0f;		//分子

		cublasDdot(cublasHandle, N + N + N, dp, 1, dr, 1, &Numerator);
	
		alpha = 1.0;
		beta = 0.0f;
		cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, 3 * N, 3 * N, total, &alpha, descr, d_CSR_Kval, d_CSR_row, d_CSR_col, dp, &beta, TempVector);
		cublasDdot(cublasHandle, N + N + N, dp, 1, TempVector, 1, &Denominator);
		
		CorrectionCoefficientA = Numerator / Denominator;
		
		cublasDaxpy(cublasHandle, N + N + N, &CorrectionCoefficientA, dp, 1, dx, 1);
	
		alpha = -CorrectionCoefficientA;
		beta = 1;	
		cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, 3 * N, 3 * N, total, &alpha, descr, d_CSR_Kval, d_CSR_row, d_CSR_col, dp, &beta, dr);
		cublasDdot(cublasHandle, N + N + N, dr, 1, dr, 1, &ResidualError);
			
		Iteration++;
		error[Iteration] = ResidualError;
		//printf("iteration:\t%d\tresidual error:\t%1.20f\n", Iteration, error[Iteration]);
		//if (ResidualError < EPS*EPS){
		//	double  stability = fabs(error[Iteration] - error[Iteration - 1000]);
		//	if (stability < EPS*EPS / 10){
		//		break;		//EPSの判断
		//	}
		//}
		if (ResidualError <= 1e-15f*1e-15f)break;


		alpha = 1.0f;
		beta = 0.0f;
		cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, 3 * N, 3 * N, total, &alpha, descr, d_CSR_Kval, d_CSR_row, d_CSR_col, dp, &beta, TempVector);
		cublasDdot(cublasHandle, N + N + N, dr, 1, TempVector, 1, &Numerator);
		cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, 3 * N, 3 * N, total, &alpha, descr, d_CSR_Kval, d_CSR_row, d_CSR_col, dp, &beta, TempVector);
		cublasDdot(cublasHandle, N + N + N, dp, 1, TempVector, 1, &Denominator);
		CorrectionCoefficientB = -Numerator / Denominator;
		cublasDscal(cublasHandle, N + N + N, &CorrectionCoefficientB, dp, 1);
		CorrectionCoefficientB = 1;
		cublasDaxpy(cublasHandle, N + N + N,&CorrectionCoefficientB , dr, 1, dp, 1);
	}
	end = clock();
	printf("処理時間-CG　%d\n", end - start);
	cudaMemcpy(hx, dx, sizeof(double)*(N + N + N), cudaMemcpyDeviceToHost);
	for (int i = 0; i < N + N + N; i++){
		b[i] = hx[i] / sqrt(TempDiagonal[i]);
	}

	FILE *fp_errors;
	errno_t errors;
	char file_name[100] = {};
	printf("FILE NAME:");
	scanf("%s", file_name);
	if (errors = fopen_s(&fp_errors, file_name, "w") != 0){
		printf("\n file open failed \n");
	}
	for (int i = 0; i < Iteration; i++){
		fprintf(fp_errors, "%d,%31.30f\n", i,error[i]);
	}
	fclose(fp_errors);
}

