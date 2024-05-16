#ifndef SARATMscore_H_HH
#define SARATMscore_H_HH
string SARA_align(string,string,string,string);
float read_SARA_align(int *,string,string,string,int *,int*,float &,float &);
float read_SARA_align(int *,string,string,string,float &,float &);
void split_whitespace(const string &str, vector<string> &result);
#endif
