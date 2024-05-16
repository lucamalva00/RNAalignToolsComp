#ifndef RNAalign_H_H
#define RNAalign_H_H
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <cmath>
#include <sstream>
#include <unistd.h>
#include <fcntl.h>
#include "Kabsch.h"
#include "function.h"
#include "debug_function.h"
//#include"global_variable.h"
using namespace std;

//extern varibale
extern const int MAXLEN ;
extern double D0_MIN;                             //for d0
extern double Lnorm;                              //normalization length
extern double score_d8, d0, d0_search, dcu0;      //for TMscore search
extern double **score;                            //Input score table for dynamic programming
extern bool   **path;                             //for dynamic programming  
extern double **val;                              //for dynamic programming  
extern int    xlen, ylen, minlen;                 //length of proteins
extern double **xa, **ya;                         //for input vectors xa[0...xlen-1][0..2], ya[0...ylen-1][0..2]
extern int    *xresno, *yresno;                   //residue numbers, used in fragment gapless threading 
extern double **xtm, **ytm;                       //for TMscore search engine
extern double **xt;                               //for saving the superposed version of r_1 or xtm
extern char   *seqx, *seqy;                       //for the protein sequence 
extern int    *secx, *secy;                       //for the secondary structure 
extern double **r1, **r2;                         //for Kabsch rotation 
extern double t[3], u[3][3];                      //Kabsch translation vector and rotation matrix
extern char out_reg[10000];
extern double Lnorm_ass, Lnorm_d0, d0_scale, d0A, d0B, d0u, d0a;
extern bool o_opt, a_opt, u_opt, d_opt, v_opt;
extern double TM3, TM4, TM5;

extern char   *chidx, *chidy;

//  extern variable
extern string aseq;
extern string bseq;
extern string item;

extern string achain,bchain;


const string version = "20150227";


void parameter_set4search(int xlen, int ylen);
void parameter_set4final(double len);
void parameter_set4final(double len,double );
void parameter_set4scale(int len, double d_s);

void get_xyz(string line,double& x,double& y,double& z,int& no,char& seq);
//void read_PDB(ifstream& infile, double** coor,string chainid,int *resno,string &seq,string atom, bool c_opt);
void read_PDB(ifstream& infile, double** &coor,string chainid,int *&resno,string & seq,char* &chid,string item,bool c_opt);
void print_help();
//void load_PDB(bool c_opt1, bool c_opt2); 
//void load_PDB(string aname,string bname,bool Ac_opt,bool Bc_opt,double ** acoor,
//		double ** bcoor,int* aresno,int *bresno,string& aseq,string& bseq,
//		string& item,string & achain, string &bchain);
void load_PDB(string& aname,string& bname,bool Ac_opt,bool Bc_opt,double ** &acoor,
		double ** &bcoor,int* &aresno,int * &bresno,string& aseq,string& bseq,
		string &item,string &achain, string &bchain);

int get_APDB_length(ifstream& name,bool c_opt,string& item, string& achian);
int get_BPDB_length(ifstream& name,bool c_opt,string& item, string& bchain);


template <class A> void NewArray(A*** array, int Narray1,int Narray2){
	*array = new A* [Narray1];
	for (int i=0;i<Narray1;i++)
		*(*array + i) = new A [Narray2];
};

template <class A> void DeleteArray(A ** array, int Narray){
	for(int i=0; i<Narray; i++)
		if(*(*array+i)) delete [] *(*array+i);
	if(Narray) delete [] (*array);
	(*array)=NULL;
};


void print_2Array(double** &array2,int n1,int n2) ;

int getmin(int alen,int blen);

double detailed_search( double **x,
                        double **y, 
                        int x_len, 
                        int y_len, 
                        int invmap0[],
                        double t[3],
                        double u[3][3],
                        int simplify_step,
                        int score_sum_method                        
                       );

int score_fun8( double **xa, 
                double **ya, 
                int n_ali,
                double d,
                int i_ali[], 
                double *score1,
                int score_sum_method
              );

double TMscore8_search( double **xtm, 
                        double **ytm,
                        int Lali, // aligned length of residue
                        double t0[3],
                        double u0[3][3],
                        int simplify_step,
                        int score_sum_method,
                        double *Rcomm
                       );

double get_score_fast(double **x, double **y, int x_len, int y_len, int invmap[]);

double get_initial( double **x, 
                    double **y, 
                    int x_len,
                    int y_len, 
                    int *y2x
                   );

bool get_initial_local(  double **x, 
						 double **y, 
						 int x_len,
						 int y_len, 
						 int *y2x
						 );

void score_matrix_rmsd(  double **x, 
						 double **y, 
						 int x_len,
						 int y_len,
						 int *y2x
						 );

void score_matrix_rmsd_sec(  double **x, 
							 double **y, 
							 int x_len,
							 int y_len,
							 int *y2x
							 );

void get_initial_ssplus( double **x, 
						 double **y, 
						 int x_len,
						 int y_len,
						 int *y2x0,
						 int *y2x						
						 );

void find_max_frag(double **x, int *resno, int len, int *start_max, int *end_max);

double get_initial_fgt( double **x, 
						double **y, 
						int x_len,
						int y_len, 
						int *xresno,
						int *yresno,
						int *y2x
						);

double DP_iter( double **x,
                double **y, 
                int x_len, 
                int y_len, 
                double t[3],
                double u[3][3],
				int invmap0[],
				int g1,
				int g2,
				int iteration_max                                   
				);

void NWDP_TM(int len1, int len2, double gap_open, int j2i[]);

void NWDP_TM(double **x, double **y, int len1, int len2, double t[3], double u[3][3], double d02, double gap_open, int j2i[]);

void NWDP_TM(int *secx, int *secy, int len1, int len2, double gap_open, int j2i[]);

void PrintErrorAndQuit(string sErrorString);

bool locate_memory_temp(int alen, int blen);

void print_1Array(int *array, int size);
void print_1Array(char *array, int size);

double TMSCORE(double t[3], double u[3][3], int n_ali, double* &distance);

void get_initial_splus(double **x, double **y, int x_len, int y_len, int *y2x0 ,int *y2x) ;
//void score_matrix_rmsd_s(x, y, x_len, y_len, y2x0);


void output_alignment(string xname,string yname,int x_len,int y_len,double t[3],double u[3][3],
		double TM1,double TM2,double rmsd,double d0_out,int m1[], int m2[],int n_ali8,
		int n_ali,double TM_0,double Lnorm_0,double d0_0, int* & ano,int* & bno,bool a_opt,
		bool u_opt, bool _opt,bool o_opt,char superposed_file[100]);

int read_dssra(string name, string chain_id);

int read_dssrb(string name, string chain_id);

int read_b_second_structual(string name,string chainid);

int read_a_second_structual(string name,string chainid);
template<class T> void print_1Array(T * arrat, int size);
void free_memory();
void get_xyz(string line,double& x,double& y,double& z);
void parameter_set4rowsearch(int xlen, int ylen);
#endif
