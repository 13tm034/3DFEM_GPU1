#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"Header.h"

void declare_check(void){
	FILE *fp;
	errno_t error;
	char file_name[] = "log.txt";

	if (error = fopen_s(&fp, file_name, "a") != 0){
		printf("\n file open failed \n");
	}

	fprintf(fp, "----------------------DECLARE CHECK START----------------------\n");
	fclose(fp);
}


void declare_check1(int E,int N,int M){
	FILE *fp;
	errno_t error;
	char file_name[] = "log.txt";

	if (error = fopen_s(&fp, file_name, "a") != 0){
		printf("\n file open failed \n");
	}

	fprintf(fp, "elemental number E %d\n",E);
	fprintf(fp, "node number N %d\n", N);
	fprintf(fp, "elemental number M %d\n", M);

	fclose(fp);
}


void declare_check2(int E, int N, int M){
	FILE *fp;
	errno_t error;
	char file_name[] = "log.txt";

	if (error = fopen_s(&fp, file_name, "a") != 0){
		printf("\n file open failed \n");
	}

	fprintf(fp, "%dÅ~%d stiffness matrix\n", N,N);
	fprintf(fp, "%d vector for strain\n", N);

	fclose(fp);
}

void Kmatrix_log(void){

	FILE *fp;
	errno_t error;
	char file_name[] = "log.txt";

	if (error = fopen_s(&fp, file_name, "a") != 0){
		printf("\n file open failed \n");
	}

	fprintf(fp, "----------------------KMATRIX CHECK START----------------------\n");

	fclose(fp);


}


void Kematrix_log(int i){

	FILE *fp;
	errno_t error;
	char file_name[] = "log.txt";

	if (error = fopen_s(&fp, file_name, "a") != 0){
		printf("\n file open failed \n");
	}

	fprintf(fp, "ELEMENTAL %d\n",i);

	fclose(fp);


}


void Kematrix_log1(double Ke[][24],int i){

	FILE *fp;
	errno_t error;
	char file_name[] = "log.txt";

	if (error = fopen_s(&fp, file_name, "a") != 0){
		printf("\n file open failed \n");
	}


	for (int j = 0; j < 24; j++){
		for (int k = 0; k < 24; k++){
			fprintf(fp, "%lf,", Ke[j][k]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	fclose(fp);


}
void Kmatrix_BC(double **K ,int N){

	FILE *fp;
	errno_t error;
	char file_name[] = "log.txt";

	if (error = fopen_s(&fp, file_name, "a") != 0){
		printf("\n file open failed \n");
	}

	fprintf(fp, "\nKmatrix_BC\n");

	for (int i = 0; i < N + N + N; i++){
		for (int j = 0; j < N + N + N; j++){
			fprintf(fp,"%lf,", K[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

}

void S_BC(double *S, int N){

	FILE *fp;
	errno_t error;
	char file_name[] = "log.txt";

	if (error = fopen_s(&fp, file_name, "a") != 0){
		printf("\n file open failed \n");
	}

	fprintf(fp, "S_BC");

	for (int i = 0; i < N + N + N; i++){
		fprintf(fp, "%lf\n", S[i]);
	}
	fclose(fp);

}