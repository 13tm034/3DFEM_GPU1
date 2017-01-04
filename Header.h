//20160516�@�f�[�^�\���̌���
//20160519	B�}�g���N�X�쐬��dndx�̌v�Z�̑O�܂Ł�dndx�̌v�Z���Ƃ肠����B�}�g���N�X����
//20160520	�֐��̑��K�w����D�}�g���N�X��BDB�̎�O�܂ŁiB�}�g���N�X�����Ŋ����Ă��Ȃ����Ƃɒ��Ӂj
//20160523	�\�[�X�t�@�C���̐������֐�Bmatrix�̍ו�����B�}�g���N�X�����Ŋ���B�i�K�v�Ȃ��H�Ǝv����j���v�f�����s���S�̍����s��ɉ��Z�����E�����̓K�p
//20160524	���E�����̓K�p���A�����������������v�f�����s��Ɋe�v�f�̑̐ς�������
//20160609	�ϕ����@�����������̂Ńf�o�b�O�J�n�A�����w������s���ꕔ�������Č��݂ł���Ƃ���̓f�o�b�O�������K�E�X�ϕ��̎������@���l����B
//20160610	5���̃K�E�X�ϕ��������������������f�o�b�O���\���o����(�K�E�X�W�������_��)
//20160613	malloc�֐��̎��������[�h����������
//20160614	�S�v���O�����̎��s��jacob�̕����Ɍ�肪����
//20160615  jacob���f�o�b�O�J�n
//20160616	jacob�̃o�O�𔭌���Kematrix�Ƀo�O�����肻��
//20160617	Kematrix���Ӄf�o�b�O
//20160620	Bmatrix�̍Ō�ŃX�^�b�N�A����͖���������ꂽ�ψʂ͂������������E�����K�����K�͖��Ȃ��������ꉞ���ʂ͏o����paraview�Ő}��
//20160622	�V���ȏ����ōČv�Z
//20160624	�v�Z���ʂɓ��Asolve�̉���MATLAB�Ɠ����ł͂Ȃ��B
//20160629	�����������̌��ʂ���������
//20160630	�f�[�^�\�������R����ɓn�����߂ɍ\���̕������m�F
//20160711	���v���z�@�\���o�[������
//20160712	�O�����̎����A�t�@�C������̃��[�h����������
//20160719	�i�X�g�����t�@�C���ǂݍ��ݕ����̍쐬
//20160720  ����
//20160721	�ߓ_�ԍ��̂����Ƀi�X�g�����ƃv���O�����ō�������dn/dr�̕�����ς���B
//20160721	�ߓ_�ԍ��̂�����dn/dr�̕��������ɖ߂��A���͂̎��_�ł�����
//20160726	�ߓ_�ԍ��̂����̎��s
//20160811	�f�o�b�O���i�X�g������v�f�ł͔��U�����R����Ɋm�F
//20160907  �ψʂ͐������o����悤�ɂȂ��������ɂ͉��́A�Ђ���
//20160914	NASTRAN�̃m�[�h���K���������Ԓʂ�łȂ����Ƃ����o��NASTRANnode����{���̃m�[�h���Q�Ƃł���悤�ɂ���
//20161003  �؍�ɉ������v�f�����̃V�X�e���̎����J�n
//20161014	���b�V���ύX���ɋC�����邱�Ɓ��v�f�ߓ_�ԍ��A�����O���̏搔
//20161025	�\���o�[�����O����������
//20161026	�v�f�����s��t�@�C���ۑ��̎���
//20161107	CSR�`���ϊ��̎���
//20161122	���v���z�@GPU�����ACOO�@���@CSR�̎����ACSR�����s��ɋ��E������K�p����
//20161201  CSR�`�����v���z�@�̎������v�f����
//20161206  CSR�ł̗v�f�����̋�̉�
#include<time.h>
#include<string.h>

#ifndef _HEADER_H_
#define _HEADER_H_

/*macros*/
//#define E 1		//�v�f��
//#define N E+7	//�ߓ_��
//#define M 1		//�ޗ���
#define EPS 0.0000000001
#define BLOCK 16

#define DISPLAY(string) printf("---------------------------------------------------------%s---------------------------------------------------------\n",string);


/*structure*/
typedef struct nodes{

	double x[3];		//���W
	int num;
	int xrc[3];			//�S������
	double xd[3];		//�ψ�
	double xf[3];		//�ߓ_��
	int point;			//������v�f�̐�

}node;

typedef struct elements{

	int node[8];		//�v�f�ߓ_���W
	double stress[6];	//����
	double strain[6];   //�c
	//double stiff[24][24];
	double mises;		//�~�[�[�X����
	bool rm;			//false 
	int m;				//�ޗ��ԍ�
	int point;
	bool Restraint;
}element;

typedef struct materials{
	double v;		//�|�A�\����
	double e;		//�c�e���W��
}material;




/*function*/

/*level 1*/
void Kmatrix(node *no, element *el,material *m,double*COO_Kval,int *COO_col,int *COO_row,int E,int N,int M,int *RealNumOfValues);	//�S�̍����v�f�}�g���N�X
void COO2CSR(int *preCSR_row, int *COO_row,int *CountRow,int RealNumberOfValues);
void setBC(node *no, double *S, double **K,int E,int N,int M);						//���E�����K�p
void setBC_CSR(node *no, double *S, double *CSR_Kval, int *CSR_col, int *CSR_row,int E,int N,int M);
void CSR2BC_CSR(double *CSR_Kval, int *CSR_col, int *CSR_row, double *BC_CSR_Kval, int *BC_CSR_col, int *BC_CSR_row,int N);
/*level 2*/
void Kematrix(node *no, element *el, material *m,double Ke[][24],int i);			//�v�f�����}�g���N�X
void addK(element *el, double Ke[][24], double **K, int i,FILE *fp);							//�v�f�����s���S�̍����s��ɉ��Z
void load_stress(node *no, double *S,int N);												//���̓��[�h
void set_RC(node *no, double *S, double **K,int E,int N,int M);						//�ψʍS���ݒ�
void AddKe2preCOO(element *el, double Ke[][24], double *COO_Kval, int *COO_col, int *COO_row, int MaximumNumberOfValues, int ElementNumber);
void Sort_preCOO(double *COO_Kval, int *COO_col, int *COO_row,int MaximumNumberOfValues);		//�\�[�g�̌Ăяo���֐�
void preCOO2COO(double *COO_Kval, int *COO_col, int *COO_row,int MaximumNumberOfValues ,int *RealNumberOfValues);	//preCOO ����COO
void Sort_preCOO_col(double *COO_Kval, int *COO_col, int *COO_row,int StartIndex,int EndIndex);
void set_RC_CSR(node *no, double *S, double *CSR_Kval, int *CSR_col, int *CSR_row, int E, int N, int M);
/*level 3*/
void Bmatrix(node *no, element *el, double B[][24],int i,double *r,double J[][3]);	//B�}�g���N�X
void Dmatrix(double D[][6], material *m,int i,element *el);							//D�}�g���N�X
void BDBmatrix(double B[][24], double D[][6],double Ke[][24]);						//BDB�̌v�Z�i�����IKe�̌v�Z�j
void set_eRC(double **K, node *no, double *S, double *temp,int i,int j,int N);		//���E�����̓K�p
void set_eRC_CSR(node *no,double *S,double *TempVector,double *CSR_Kval, int *CSR_col, int *CSR_row, int i, int j, int N);
void det_jacob(double J[][3], double *detJ);										//jacob�s���det
void product_detJ_weight(double pKe[][24], double detJ, double *w);					//detJ��weight��pKe�ɂ�����
void addKe(double pKe[][24], double Ke[][24]);										//pKe��Ke�̗v�f�ɑ���
void QsortPreCOO(double *COO_Kval, int *COO_col, int *COO_row, int left, int right); //quick sort �̊֐�row�C���f�b�N�X
/*level 4*/
void rollingB(double Bt[][6], double B[][24]);										//B�̓]�n
void productBtD(double Bt[][6], double D[][6],double BD[][6]);						//BtD�̐�
void productBDB(double BD[][6], double B[][24], double Ke[][24]);					//BDB�̐�
void jacob(double J[][3], double dndr[][8], double nd[][8],double *r);				//jacob�s��
void inv_jacob(double invJ[][3], double J[][3]);									//jacob�s��̋t�s��
void DnDx(double dndx[][8], double invJ[][3],double dndr[][8]);						//DnDx�̍쐬

//solve
void solve_matrix(double **K, double *S, int N,node *no,int *PreCount);									//�}�g���N�X�̃\���o
void solve_matrix_gpu(double **A, double *b, int N, node *no);				//gpu�ɂ�����(���́@dense)
void solve_matrix_gpu_CSR(double *CSR_Kval, int *CSR_col, int *CSR_row, double *b, int N, node *no, int RealNumberOfValues);	//gpu�ɂ�����(���́@CSR)
//input_outputload
void difference_output(double *S, int N);
void input(node *no, element *el, material *m);										//�e���̓���(�����)
void original_output(node *no, int N);


//load
void load_number(int *N, int *E,int *M);											//�v�f���A�ߓ_���̓���
void info_N(node *no,int N);																//�ߓ_���̓ǂݍ���
void info_E(element *el,int E);															//�v�f���̓ǂݍ���
void info_M(material *m);

//initial
void initial( double *S, node *no, element *el, material *ml,int N,int E,int M); //���I�ȕϐ��̏�����


//stress_strain
void strain_stress_calc(node *no, element *el, material *m, int E, int N, int M);
void addB(double B[][24], double b[][24], double *w, double detJ);
void addDB(double DB[][24], double b[][24], double D[][6], double *w, double detJ);
void calc_strain(double B[][24], element *el, node *no, int i);
void mises_strain(element *el, int i,material *m);
void calc_stress(double DB[][24], element *el, node *no, int i);

//log
void declare_check(void);
void declare_check1(int E, int N, int M);
void declare_check2(int E, int N, int M);
void Kmatrix_log(void);
void Kematrix_log(int i);
void Kematrix_log1(double Ke[][24], int i);
void Kmatrix_BC(double **K, int N);
void S_BC(double *S, int N);


//timer
void time_display(char *s);


//kari
void inv_K(double **K);


//RemoveMatrix
void remove_matrix(int ElementToBeRemoved[][100], int num, double **K, element *el, node *no, material *m);
void removeK(element *el, double Ke[][24], double **K, int i,node *no);
void convert(char *readline, double Ke[][24]);
void RemoveMatrixCSR(int ElementToBeRemoved[][100], int num, double *CSR_Kval, int *CSR_col, int *CSR_row, element *el, node *no, material *m);
void RemoveNodesCSR(double *BC_CSR_Kval, int *BC_CSR_col, int *BC_CSR_row,double *S, node *no,int N,int BC_RealNum);
//output
void coordinate_output(node *no, int N);
void element_node_output(element *el, int E);
void result_output(element *el, node *no, double *S, int E, int N,int count,int NumberOfRemoved);

#endif