#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"Header.h"

void  initial( double *S, node *no, element *el, material *ml,int N,int E,int M){

	

	for (int i = 0; i < N+N+N; i++){
		S[i] = 0;

	}

	//for (int i = 0; i < N + N + N; i++){
	//	for (int j = 0; j < N + N + N; j++){
	//		K[i][j] = 0;
	//	}
	//}

	for (int i = 0; i < N; i++){
		for (int j = 0; j < 3; j++){
			no[i].x[j] = 0;
			no[i].xd[j] = 0;
			no[i].xf[j] = 0;
			no[i].xrc[j] = 0;
		}
		no[i].point = 0;
	}


	for (int i = 0; i < E; i++){
		el[i].m = 0;
		el[i].mises = 0;
		el[i].rm = false;
		for (int j = 0; j < 8; j++){
			el[i].node[j] = 0;
		}

		for (int j = 0; j < 6; j++){
			el[i].strain[j] = 0;
			el[i].stress[j] = 0;
		}
		el[i].rm = false;
	}

	for (int i = 0; i < M; i++){
		ml[i].e = 0;
		ml[i].v = 0;
	}







}
