#ifndef global_variable_H_H
#define global_variable_H_H

#include<string>
#include<vector>

const int MAXLEN = 10000;
// globale variable
int alen(3),blen;  //  for locate the length of RNA length
std::string aseq,bseq; //  for RNA sequence 
std::string aname,bname,achain(""),bchain(""); // for A B PDB file and chain ids
std::string item("C3'");
double  u0;
char* A0;
// argument variable  add with the pro

double D0_MIN;                             //for d0
double Lnorm;                              //normalization length
double score_d8, d0, d0_search, dcu0;      //for TMscore search
double **score;            			       //Input score table for dynamic programming
bool   **path;                             //for dynamic programming  
double **val;                              //for dynamic programming  
int    xlen, ylen, minlen;                 //length of proteins
double **xa, **ya;                         //for input vectors xa[0...xlen-1][0..2], ya[0...ylen-1][0..2]
                                           //in general, ya is regarded as native structure --> superpose xa onto ya
int    *xresno, *yresno;                   //residue numbers, used in fragment gapless threading 
double **xtm, **ytm;                       //for TMscore search engine
double **xt;                               //for saving the superposed version of r_1 or xtm

char   *chidx, *chidy;
char   *seqx, *seqy;                       //for the RNA sequence 
int    *secx, *secy;                       //for the secondary structure 
// no useful for this program

double **r1, **r2;                         //for Kabsch rotation 
double t[3], u[3][3];                      //Kabsch translation vector and rotation matrix

char out_reg[MAXLEN];
double Lnorm_ass, Lnorm_d0, d0_scale, d0A, d0B, d0u, d0a;
bool o_opt, a_opt, u_opt, d_opt, v_opt;
double TM3, TM4, TM5;



#endif
