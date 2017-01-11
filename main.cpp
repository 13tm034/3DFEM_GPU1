#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"Header.h"


int main(void){


	printf("    _/      _/                      _/        _/                            _/                          _/      _/                  _/       \n");
	printf("   _/_/  _/_/    _/_/_/    _/_/_/  _/_/_/        _/_/_/      _/_/        _/_/_/_/    _/_/      _/_/    _/      _/          _/_/_/  _/_/_/    \n");
	printf("  _/  _/  _/  _/    _/  _/        _/    _/  _/  _/    _/  _/_/_/_/        _/      _/    _/  _/    _/  _/      _/        _/    _/  _/    _/   \n");
	printf(" _/      _/  _/    _/  _/        _/    _/  _/  _/    _/  _/              _/      _/    _/  _/    _/  _/      _/        _/    _/  _/    _/    \n");
	printf("_/      _/    _/_/_/    _/_/_/  _/    _/  _/  _/    _/    _/_/_/          _/_/    _/_/      _/_/    _/      _/_/_/_/    _/_/_/  _/_/_/	     \n");

	printf("ooooooooooo o       oooooooo8 ooooooooooo      ooooooooooo ooooooooooo oooo     oooo       oooooooo8 ooooo  oooo oooooooo8 ooooooooooo ooooooooooo oooo     oooo \n");
	printf(" 888    88 888     888        88  888  88       888    88   888    88   8888o   888       888          888  88  888        88  888  88  888    88   8888o   888  \n");
	printf(" 888ooo8  8  88     888oooooo     888           888ooo8     888ooo8     88 888o8 88        888oooooo     888     888oooooo     888      888ooo8     88 888o8 88  \n");
	printf(" 888     8oooo88           888    888           888         888    oo   88  888  88               888    888            888    888      888    oo   88  888  88  \n");
	printf("o888o  o88o  o888o o88oooo888    o888o         o888o       o888ooo8888 o88o  8  o88o      o88oooo888    o888o   o88oooo888    o888o    o888ooo8888 o88o  8  o88o \n");

	/*load number*/
	DISPLAY("LOAD NUMBER");
	int N, *pN, E, *pE, M, *pM;
	N = 0;
	E = 0;
	M = 0;
	pN = &N;		//節点数
	pE = &E;		//要素数
	pM = &M;		//材料種数
	
	load_number(pN, pE, pM);
	printf("NODE NUMBER:\tN=%d\n", N);
	printf("ELEMENTAL NUMBER:\tE=%d\n", E);
	printf("MATERIAL NUMBER:\tM=%d\n", M);


	/*declaret*/
	DISPLAY("DACLARE");
	declare_check();
	node *no;		//節点情報
	element *el;	//要素情報
	material *m;	//材料種
	no = (node *)malloc(N*sizeof(node));
	el = (element *)malloc(E*sizeof(element));
	m = (material *)malloc(M*sizeof(material));
	declare_check1(E, N, M);
	printf("no[%d]\n", N);
	printf("el[%d]\n", E);
	printf("m[%d]\n", M);


	double *S;
	S = (double *)malloc(sizeof(double)*(N + N + N));
	double *tempS;
	tempS = (double *)malloc(sizeof(double)*(N + N + N));
	declare_check2(E, N, M);
	printf("%d*%d STIFFNES MATRIX\n", N+N+N,N+N+N);
	printf("%d STRAIN VECTOR\n", N + N + N);
	
	/*initialization*/
	initial(S, no, el, m, N, E, M);
	
	/*load_node_element*/
	DISPLAY("LOAD");
	info_N(no,N);
	int *NASTRAN_no;
	NASTRAN_no = (int *)malloc(sizeof(int)*(no[N-1].num));
	for (int i = 0; i < no[N - 1].num; i++){
		NASTRAN_no[i] = 0;
	}
	for (int i = 0; i < N; i++){
		NASTRAN_no[no[i].num] = i;
	}
	

	info_E(el,E);
	

	for (int i = 0; i < E; i++){
		for (int j = 0; j < 8; j++){
			el[i].node[j] = NASTRAN_no[el[i].node[j]];
		}
	}
	//要素節点番号をNASTRANnodeから本来のノードに変換+ある節点を共有する要素の数

	int element[8] = {};

	for (int i = 0; i < E; i++){
		for (int j = 0; j < 8; j++){
			element[j] = el[i].node[j];
			no[element[j]].point++;
		}
		el[i].node[0] = element[4];
		el[i].node[1] = element[7];
		el[i].node[2] = element[6];
		el[i].node[3] = element[5];
		el[i].node[4] = element[0];
		el[i].node[5] = element[3];
		el[i].node[6] = element[2];
		el[i].node[7] = element[1];//16

		//el[i].node[0] = element[3];
		//el[i].node[1] = element[7];
		//el[i].node[2] = element[4];
		//el[i].node[3] = element[0];
		//el[i].node[4] = element[2];
		//el[i].node[5] = element[6];
		//el[i].node[6] = element[5];
		//el[i].node[7] = element[1]; //900


		//el[i].node[0] = element[0];
		//el[i].node[1] = element[1];
		//el[i].node[2] = element[2];
		//el[i].node[3] = element[3];
		//el[i].node[4] = element[4];
		//el[i].node[5] = element[5];
		//el[i].node[6] = element[6];
		//el[i].node[7] = element[7];
		
		
	}


	for (int i = 0; i < E; i++){
		for (int j = 0; j < 8; j++){
			int nd = el[i].node[j];
			el[i].point += no[nd].point;
			for (int xyz = 0; xyz < 3; xyz++){
				if (no[nd].xrc[xyz] == 1)el[i].Restraint = true;
			}
		}
	}

	info_M(m);
	printf("N=%d,E=%d,M=%d\n", N, E, M);
	for (int i = 0; i < N; i++){
		printf("node=%d, x=%10f,%10f,%10f, xf=%10f,%10f,%10f, xrc=%d,%d,%d, point=%d\n"
			, i, no[i].x[0], no[i].x[1], no[i].x[2], no[i].xf[0], no[i].xf[1], no[i].xf[2]
			, no[i].xrc[0], no[i].xrc[1], no[i].xrc[2],no[i].point);
	}

	for (int i = 0; i < E; i++){
		printf("element:%d material: %d point:%d Restraint:%d (",i,el[i].m,el[i].point,el[i].Restraint);

		for (int j = 0; j < 8; j++){

			printf("%d,", el[i].node[j]);

		}
		printf(")\n");

	}

	printf("e=%lf,v=%lf\n",m[0].e,m[0].v);
	original_output(no, N);

	/******************************Kmatrix start********************************/
	DISPLAY("STIFFNESS MATRIX");
	/*preCOOの宣言*/
	int RealNumberOfValues = 0;
	int MaximumNumberOfValues = 24 * 24 * E;	//値がある要素の理論上最大値					 
	double *preCOO_Kval;
	preCOO_Kval = (double *)malloc(sizeof(double)*(MaximumNumberOfValues));
	int *preCOO_col;
	preCOO_col = (int *)malloc(sizeof(int)*(MaximumNumberOfValues));
	int *preCOO_row;
	preCOO_row = (int *)malloc(sizeof(int)*(MaximumNumberOfValues));
	for (int i = 0; i < MaximumNumberOfValues; i++){
		preCOO_Kval[i] = 0;
		preCOO_col[i] = 0;
		preCOO_row[i] = 0;
	}
	/*剛性行列の作成*/
	Kmatrix(no, el, m,preCOO_Kval,preCOO_col,preCOO_row, E, N, M,&RealNumberOfValues);

	/*COO形式に変換*/
	double *COO_Kval;
	COO_Kval = (double *)malloc(sizeof(double)*(RealNumberOfValues));
	int *COO_col;
	COO_col = (int *)malloc(sizeof(int)*(RealNumberOfValues));
	int *COO_row;
	COO_row = (int *)malloc(sizeof(int)*(RealNumberOfValues));
	for (int i = 0; i < RealNumberOfValues; i++){
		COO_Kval[i] = preCOO_Kval[i];
		COO_col[i] = preCOO_col[i];
		COO_row[i] = preCOO_row[i];
	}
	free(preCOO_Kval);
	free(preCOO_col);
	free(preCOO_row);


	/*preCSR*/
	printf("RealNum=%d\n", RealNumberOfValues);
	int *preCSR_row;
	preCSR_row = (int *)malloc(sizeof(int)*(RealNumberOfValues));
	if (preCSR_row == NULL){
		printf("malloc failed\n");
		system("PAUSE");
	}
	for (int i = 0; i < RealNumberOfValues; i++){
		preCSR_row[i] = 0;
	}

	int CountRow = 0;
	COO2CSR(preCSR_row, COO_row,&CountRow, RealNumberOfValues);

	double *CSR_Kval;
	CSR_Kval = (double *)malloc(sizeof(double)*(RealNumberOfValues));
	int *CSR_col;
	CSR_col = (int *)malloc(sizeof(int)*(RealNumberOfValues));
	int *CSR_row;
	CSR_row = (int *)malloc(sizeof(int)*(CountRow));
	for (int i = 0; i < RealNumberOfValues; i++){
		CSR_Kval[i] = COO_Kval[i];
		CSR_col[i] = COO_col[i];
	}
	for (int i = 0; i <= CountRow; i++){
		CSR_row[i] = preCSR_row[i];
	}
	printf("\n");
	for (int i = CSR_row[3099]; i <CSR_row[3100]; i++){
		if (CSR_col[i] == 1703)printf("CSRK=%lf\n", CSR_Kval[i]);
	}
	system("PAUSE");
	free(COO_Kval);
	free(COO_col);
	free(COO_row);
	printf("count=%d,RealNum=%d\n", CSR_row[CountRow], RealNumberOfValues);
	DISPLAY("SET PERMANENT BOUNDARY CONDITIONS")
	//setBC(no, S, K, E, N, M);
	setBC_CSR(no, S, CSR_Kval, CSR_col, CSR_row, E, N, M);

	int BC_RealNumberOfValues = 0;

	for (int i = 0; i < RealNumberOfValues; i++){
		if (CSR_Kval[i] != 0 && CSR_col[i]!=0 )BC_RealNumberOfValues++;
	}


	double *BC_CSR_Kval;
	BC_CSR_Kval = (double *)malloc(sizeof(double)*(BC_RealNumberOfValues));
	int *BC_CSR_col;
	BC_CSR_col = (int *)malloc(sizeof(int)*(BC_RealNumberOfValues)); 
	int *BC_CSR_row;
	BC_CSR_row = (int *)malloc(sizeof(int)*(N+N+N+1));

	CSR2BC_CSR(CSR_Kval, CSR_col, CSR_row, BC_CSR_Kval, BC_CSR_col, BC_CSR_row, N);
	for (int i = BC_CSR_row[3099]; i < BC_CSR_row[3100]; i++){
		if (BC_CSR_col[i] == 1703)printf("BC_CSRK=%lf", BC_CSR_Kval[i]);
	}

	S_BC(S, N);

	/******************************solve matrix********************************/


	int ElementToBeRemoved[2][100] = {}; //第一列→消す要素数,第二列→実際に消す要素番号
	ElementToBeRemoved[0][0] =0;
	//ElementToBeRemoved[0][1] = 350;
	ElementToBeRemoved[1][0] = 1;
	ElementToBeRemoved[1][1] = 350;
	int ETRcount = 2;				//ElementToBeRemovedの行数
	int count = 0;					//行の位置
	int NumberOfRemoeved = 0;		//削除した要素の数
	int PreCount = 0;
	

	for (int i = 0; i < N + N + N; i++){
		for (int j = BC_CSR_row[i]; j < BC_CSR_row[i + 1]; j++){
			if(i==BC_CSR_col[j])printf("%d,%d\t%lf\n", i, BC_CSR_col[j], BC_CSR_Kval[j]);
		}
	}

	while (count != ETRcount){
		DISPLAY("SOLVING MATRIX");
		NumberOfRemoeved += ElementToBeRemoved[count][0];
		RemoveMatrixCSR(ElementToBeRemoved, count, BC_CSR_Kval, BC_CSR_col, BC_CSR_row, el, no, m);
		RemoveNodesCSR(BC_CSR_Kval, BC_CSR_col, BC_CSR_row,S, no,N,BC_RealNumberOfValues);
		clock_t start, end;
		start = clock();
		
		//GPUsolver
		DISPLAY("RUN GPU SOLVER");
		start = clock();
		solve_matrix_gpu_CSR(BC_CSR_Kval, BC_CSR_col, BC_CSR_row, S, N, no, BC_RealNumberOfValues);
		end = clock();
		printf("処理時間　%d\n", end - start);

		strain_stress_calc(no, el, m, E, N, M);
		result_output(el, no, S, E, N, count, NumberOfRemoeved);
		system("PAUSE");


		count++;
	}
}


