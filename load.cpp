#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"Header.h"


#define word 256

/****************************CLASS_NODE*******************************/
class Node{
	int number;
public:
	Node() :number(0){};
	void count();
	int finish();
	void convert_x(char *s,node *no,int count);
	void convert_f(char *s, node *no,int N);
	void convert_rc(char *s, node *no,int N);
	void convert_num(char *s, node *no,int count);
};


void Node::count()
{
	number++;
}

int Node::finish()
{
	return number;
}


void Node::convert_x(char *s, node *no, int count){

	char x1[8] = { '\0' };
	char x2[8] = { '\0' };
	char x3[8] = { '\0' };

	for (int i = 0; i < 8; i++){
		x1[i] = s[24 + i];
		x2[i] = s[32 + i];
		x3[i] = s[40 + i];
	}

	no[count].x[0] = atof(x1);
	no[count].x[1] = atof(x2);
	no[count].x[2] = atof(x3);


}


void Node::convert_num(char*s, node *no, int count){
	char num[8] = { '\0' };

	for (int i = 0; i < 8; i++){
		num[i] = s[8 + i];
	}

	no[count].num = atoi(num);
}

void Node::convert_f(char *s, node *no,int N){

	char num[8] = { '\0' };
	char force1[8] = { '\0' };
	char force2[8] = { '\0' };
	char force3[8] = { '\0' };


	for (int i = 0; i < 8; i++){
		num[i] = s[16 + i];
		force1[i] = s[40 + i];
		force2[i] = s[48 + i];
		force3[i] = s[56 + i];
	}

	for (int i = 0; i < N; i++){
		if (atoi(num) == no[i].num){
			no[i].xf[0] = atof(force1);
			no[i].xf[1] = atof(force2);
			no[i].xf[2] = atof(force3);
		}
	}

}


void Node::convert_rc(char *s, node *no,int N){

	char num[8] = { '\0' };
	char rc1[8] = { '\0' };
	char rc2[8] = { '\0' };
	char rc3[8] = { '\0' };



	for (int i = 0; i < 8; i++){
		num[i] = s[24 + i];
		}
	

	for (int i = 0; i < N; i++){
		if (atoi(num) == no[i].num){
			no[i].xrc[0] = 1;
			no[i].xrc[1] = 1;
			no[i].xrc[2] = 1;
		}
	}
}

/****************************CLASS_ELEMENT*******************************/

class Element{
	int number;
public:
	Element() :number(0){};
	void count();
	int finish();
	void convert_e(char *s,element *el, int count);
	void convert_e2(char *s, element *el, int count);
};


void Element::count()
{
	number++;
}

int Element::finish()
{
	return number;
}


void Element::convert_e(char *s,element *el, int count){

	char mat[8] = { '\0' };
	char nd[6][8] = {'\0'};


	for (int i = 0; i < 8; i++){
		mat[i] = s[16 + i] ;
	}

	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 8; j++){
			nd[i][j] = s[24 + j + 8 * i];
		}
	}


	el[count].m = atoi(mat)-1;

	for (int i = 0; i < 6; i++){
		el[count].node[i] = atoi(nd[i]);
	}

}

void Element::convert_e2(char *s, element *el, int count){

	char nd[2][8] = { '\0' };


	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 8; j++){

			nd[i][j] = s[8 + j + 8 * i];

		}
	}

	for (int i = 0; i < 2; i++){

		el[count].node[6 + i] = atoi(nd[i]);

	}



}




/****************************CLASS_MATERIAL*******************************/
class Material{
	int number;
public:
	Material() :number(0){};
	void count();
	int finish();
	void convert_m(char *s, material *m);
};


void Material::count()
{
	number++;
}

int Material::finish()
{
	return number;
}



void Material::convert_m(char *s, material *m){

	char e[8] = { '\0' };
	char v[8] = { '\0' };


	for (int i = 0; i < 8; i++){
		e[i] = s[16 + i];
		v[i] = s[32 + i];
	}

	char en[8] = { '\0' };
	char es[8] = { '\0' };
	int i = 0;
	while (e[i] != '+'){
		en[i] = e[i];
		i++;
	}

	int lim = i;
	int j = 0;
	for (i = lim-1; i < 8; i++){
		es[i] = e[j];
		j++;
	}
	printf("<<<<<<<%c\n", es[0]);

	m[0].e = atof(en)*pow(10.0, 11);
	m[0].v = atof(v);
}

/*****************************************************************************/



void load_number(int *N, int *E,int *M){


	FILE *fp;
	errno_t errors;
	char file_name[] = "mesh.dat";
	char readline[word] = { '\0' };
	Node No;
	Element El;
	Material Ml;
	int countN = 0;

	if (errors = fopen_s(&fp, file_name, "r") != 0){
		printf("\n file open failed \n");
	}

	printf("%d",*N);
	while (fgets(readline, word, fp) != NULL){


		if (readline[0] == 'G'&&readline[3] == 'D')No.count();
		else if (readline[0] == 'C'&&readline[4] == 'A')El.count();
		else if (readline[0] == '$'&&readline[24] == 'M')Ml.count();

	}


	*N = No.finish();
	*E = El.finish();
	*M = Ml.finish();


	fclose(fp);
}


void info_N(node *no,int N){

	FILE *fp;
	errno_t errors;
	char file_name[] = "mesh.dat";
	char readline[word] = { '\0' };

	Node No;

	if (errors = fopen_s(&fp, file_name, "r") != 0){
		printf("\n file open failed \n");
	}

	int count_x=0;


	while (fgets(readline, word, fp) != NULL){

		if (readline[0] == 'G'&&readline[3] == 'D'){
			No.convert_x(readline, no, count_x);
			No.convert_num(readline, no, count_x);
			count_x++;
		}
	}

	rewind(fp);

	while (fgets(readline, word, fp) != NULL){
		if (readline[0] == 'F'&&readline[4] == 'E'){
			No.convert_f(readline, no, N);
		}

		else if (readline[0] == 'S'&&readline[3] == '1'){
			No.convert_rc(readline, no, N);
		}
		
	}

	fclose(fp);


}


void info_E(element *el,int E){


	FILE *fp;
	errno_t errors;
	char file_name[] = "mesh.dat";
	char readline[word] = { '\0' };

	Element El;

	if (errors = fopen_s(&fp, file_name, "r") != 0){
		printf("\n file open failed \n");
	}


	unsigned int count = 0;

	while (fgets(readline, word, fp) != NULL){
		if (readline[0] == 'C'&&readline[4] == 'A'){
			El.convert_e(readline, el, count);
		}

		else if (readline[0] == '+'&&readline[4] == ' '&&readline[9] == ' '){
			El.convert_e2(readline, el, count);
			count++;
		}
	
	}


	//for (int i = 0; i < E; i++){
	//	printf("element:%d material: %d (", i, el[i].m);

	//	for (int j = 0; j < 8; j++){

	//		printf("%d,", el[i].node[j]);

	//	}
	//	printf(")\n");

	//}


	int element[8] = {};

	//for (int i = 0; i < E; i++){
	//	for (int j = 0; j < 8; j++){
	//		element[j] = el[i].node[j];
	//	}
	//	el[i].node[0] = element[4];
	//	el[i].node[1] = element[7];
	//	el[i].node[2] = element[6];
	//	el[i].node[3] = element[5];
	//	el[i].node[4] = element[0];
	//	el[i].node[5] = element[3];
	//	el[i].node[6] = element[2];
	//	el[i].node[7] = element[1];


	//	//el[i].node[0] = element[0];
	//	//el[i].node[1] = element[1];
	//	//el[i].node[2] = element[2];
	//	//el[i].node[3] = element[3];
	//	//el[i].node[4] = element[4];
	//	//el[i].node[5] = element[5];
	//	//el[i].node[6] = element[6];
	//	//el[i].node[7] = element[7];
	//	
	//}



}


void info_M(material *m){


	FILE *fp;
	errno_t errors;
	char file_name[] = "mesh.dat";
	char readline[word] = { '\0' };

	Material Ml;
	if (errors = fopen_s(&fp, file_name, "r") != 0){
		printf("\n file open failed \n");
	}


	while (fgets(readline, word, fp) != NULL){
		if (readline[0] == 'M'&&readline[3] == '1'){
			Ml.convert_m(readline, m);

		}
	}
}