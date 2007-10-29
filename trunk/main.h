// for primer3 input
//#include <LEDA/core/string.h>

#define _CRT_SECURE_NO_DEPRECATE 1
#define _CRT_NONSTDC_NO_DEPRECATE 1
#define NUM_MAX_GENE	7000
#define BLAST_THRESHOLD 1.00E-3
#define CLUSTER_SIZE 96  //96	24
#define SUFFIX_MAX 1 //8
#define BIG_DELETION_INSERTION 6 // min big is 6

#include "string.h"
#include "stdio.h"
#include "time.h"
#include "math.h"
#include <vector>
#include <string>
#include <iostream>
#include <streambuf>
#include <fstream>
#include <algorithm>
#include <hash_map>

#define min(a, b)  (((a) < (b)) ? (a) : (b)) 
#define max(a, b)  (((a) > (b)) ? (a) : (b)) 
#define MAX_JMP 1000

using namespace std;
using namespace stdext;

class ORF{
public:
	string sequence, upstream, downstream, left_primer, right_primer;
	string sGeneName, sTitle;
	string db_sequence, db_sequence_reverse;
	int nSuffix, nLeftExtension, nRightExtension;
	bool isActive, isPCR;
	bool isForward;

	//for blast
	string sIdentities;
	char cOrientation; // 'F' or 'R'
	double e_value;
	string no_change, normal_change, stop_codon_change, no_start_codon_change;
	string normal_insertion, frame_shift_insertion;
	string normal_deletion, frame_shift_deletion;
	int nHSP_length, num_normal_change;

public:
	ORF();
	~ORF();
	void getSequence(); // form sequence form cds, utr and primer, output to db_sequence
	string replaceStartStopCodon(string str, int flag); // replace TGA, TAG or TAA to TTA
};

// global variables
string sCDS, s5UTR, s3UTR, sPrimer, sDB;
string sSequencing, sQuery, sParsing;
basic_string <char>::size_type index;
ifstream in_file;
ofstream out_file;
hash_map <string, ORF> hORF;
hash_map <string, ORF> hCLONE;
hash_map<string, ORF> ::iterator ORF_Iter;

string sFix5 = "AAAAAAGCAGGCTCCGAATTCGCCCTT";
string sFix3 = "AAGGGCGAATTCGACCCAGCTTTCTTGTACAAA";

string sFix5_R = "AAGGGCGAATTCGGAGCCTGCTTTTTT";
string sFix3_R = "TTTGTACAAGAAAGCTGGGTCGAATTCGCCCTT";

// global functions
void build_database(); 
void build_query();
void call_blast();
void analyze_blast();

template< typename StrT >
int SplitString(const char* str, const char* delim, vector<StrT>& results);
char * _strupr_new (char * string);
void str_trim(string& str);
int parsing(string sID, string sFileName);
string getProteinName(string sDNA);