//	sequencing application

#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
#define _CRT_SECURE_NO_DEPRECATE 1

#include "main.h"
#include "class.h"
#include "file_operation.h"

using namespace std;
using namespace stdext;

int main(int argc, char* argv[])
{

	string sCategory;

	// file names
	sCDS = "TAIR7_cds_20070425.txt";
	s5UTR = "TAIR7_5_utr_20070226.txt"; 
	s3UTR = "TAIR7_3_utr_20070226.txt";
	sPrimer = "primers_all_info.txt";
	sDB = "sequencing_db.txt";
	sQuery = "sequencing_query.txt";
	sSequencing = "SEQ0012.txt";
	sParsing = "blast_result.txt";
   
	////////////////////////////////////////////
	// deal with parameters
	////////////////////////////////////////////
	if(argc<2){

		cout<<"Syntax Error: Parameter is missing!"<<endl;
		cout<<"Syntax Help: $sequencing	<parameter>"<<endl;
		cout<<"	Available parameters are:"<<endl;
		cout<<"		build_database"<<endl;
		cout<<"		build_query"<<endl;
		cout<<"		call_blast"<<endl;
		cout<<"		analyze_blast"<<endl;
		
		// for easy debug
		cout<<endl<<"...for debug..."<<endl;
		sCategory = "analyze_blast";

	}else{

		sCategory = argv[1];
		
	}// end for parameters

	if(sCategory=="build_database"){
		build_database();  
	}

	if(sCategory=="build_query"){
		build_query();  
	}

	if(sCategory=="call_blast"){
		call_blast();  
	}

	if(sCategory=="analyze_blast"){
		analyze_blast();  
	}

	return 0;
}

// read cds, upstream, downstream, primers to build blast database
void build_database(){

	string s1, s2,s3,s4, sName;
	vector<string> vLine;
	basic_string <char>::size_type index;
	int i; 
	size_t j;
	char buffer[200];

	////////////////////////////////////////////////////////////////////////////////
	// read CDS file
	////////////////////////////////////////////////////////////////////////////////
	in_file.open(sCDS.c_str());

	if(!in_file.good()){
		cout<<"ERROR: CDS file can't open!"<<endl;
		return;
	}

	s1.clear();
	s2.clear();
	sName.clear();
	i = 0;

	while(in_file.good()){

		getline(in_file, s1, '\n');
		if(s1[0]=='!' || s1[0]=='#') continue;

		vLine.clear();
		SplitString(s1.c_str(), " ", vLine); // use space to saperate columns

		if(vLine.size()>1){  // title line

			if(i%MAX_JMP==0) cout<<"processing gene CDS #"<<i<<"..."<<endl;
			
			if(!s2.empty()) {
				hORF[sName].sequence = s2;
				s2.clear();
			}

			sName = vLine[0];
			sName.erase(sName.begin());
			index = sName.find_first_of ( "." );
			s3 = sName.substr(0, index);
			s4 = sName.substr(index+1, 1);
			hORF[sName].sGeneName = s3;
			hORF[sName].nSuffix = atoi(s4.c_str());
			hORF[sName].isActive = true;
			hORF[sName].sTitle = s1;
			j = vLine.size()-1;
			if(vLine[j]=="FORWARD") hORF[sName].isForward = true;
			else hORF[sName].isForward = false;
			i++;

		}else{
			s2 = s2 + s1;
		}
	}// end while file

	// last line
	if(!s2.empty()) {
		hORF[sName].sequence = s2;
		s2.clear();
	}

	in_file.close();
	in_file.clear();

	////////////////////////////////////////////////////////////////////////////////
	// read UTR files
	////////////////////////////////////////////////////////////////////////////////
	in_file.open(s5UTR.c_str());

	if(!in_file.good()){
		cout<<"ERROR: 5' UTR file can't open!"<<endl;
		return;
	}

	s1.clear();
	s2.clear();
	sName.clear();
	i = 0;
	while(in_file.good()){

		getline(in_file, s1, '\n');
		if(s1[0]=='!' || s1[0]=='#') continue;

		vLine.clear();
		SplitString(s1.c_str(), " ", vLine); // use space to saperate columns

		if(vLine.size()>1){  // title line

			if(i%MAX_JMP==0) cout<<"processing gene 5' UTR #"<<i<<"..."<<endl;
			
			if(!s2.empty() && hORF[sName].isActive == true) {
				hORF[sName].upstream = s2;
				s2.clear();
			}

			sName = vLine[0];
			sName.erase(sName.begin());
			i++;

		}else{
			s2 = s2 + s1;
		}
	}// end while file

	// last line
	if(!s2.empty() && hORF[sName].isActive == true) {
		hORF[sName].upstream = s2;
		s2.clear();
	}

	in_file.close();
	in_file.clear();

	in_file.open(s3UTR.c_str()); 

	if(!in_file.good()){
		cout<<"ERROR: 3' UTR file can't open!"<<endl;
		return;
	}

	s1.clear();
	s2.clear();
	sName.clear();
	i = 0;

	while(in_file.good()){

		getline(in_file, s1, '\n');
		if(s1[0]=='!' || s1[0]=='#') continue;

		vLine.clear();
		SplitString(s1.c_str(), " ", vLine); // use space to saperate columns

		if(vLine.size()>1){  // title line

			if(i%MAX_JMP==0) cout<<"processing gene 3' UTR #"<<i<<"..."<<endl;
			
			if(!s2.empty() && hORF[sName].isActive == true) {
				hORF[sName].downstream = s2;
				s2.clear();
			}

			sName = vLine[0];
			sName.erase(sName.begin());
			i++;

		}else{
			s2 = s2 + s1;
		}
	}// end while file

	// last line
	if(!s2.empty() && hORF[sName].isActive == true) {
		hORF[sName].downstream = s2;
		s2.clear();
	}

	in_file.close();
	in_file.clear();

	////////////////////////////////////////////////////////////////////////////////
	// read Primer files
	////////////////////////////////////////////////////////////////////////////////
	// problem 1. primer only has locus name, no gene models
	// solution 1. try suffix from 1 to 8

	in_file.open(sPrimer.c_str());

	if(!in_file.good()){
		cout<<"ERROR: primer file can't open!"<<endl;
		return;
	}

	s1.clear();
	s2.clear();
	sName.clear();
	i = 0;

	while(in_file.good()){

		getline(in_file, s1, '\n');
		if(s1[0]=='!' || s1[0]=='#') continue;

		vLine.clear();
		SplitString(s1.c_str(), "\t", vLine); // use \t to saperate columns

		if(vLine.size()==0) continue;
		
		//cout<<"processing primer #"<<i<<"..."<<endl;

		s2 = vLine[0];

		for(int k=1; k<=SUFFIX_MAX; k++){
			_itoa( k, buffer, 10 );
			s3 = s2 + "."+buffer;
			if(hORF[s3].isActive == true){
				hORF[s3].isPCR = true;
				hORF[s3].left_primer = vLine[15];
				hORF[s3].right_primer = vLine[16];
				hORF[s3].nLeftExtension = atoi(vLine[8].c_str());
				hORF[s3].nRightExtension = atoi(vLine[9].c_str());
			}
		}
			
		i++;

	}// end while file

	in_file.close();
	in_file.clear();

	////////////////////////////////////////////////////////////////////////////////
	// setup output format 
	////////////////////////////////////////////////////////////////////////////////
	// problem 1. replace stop codons in CDS and UTR - done!
	// problem 2. add fixed sequences on both end - done!
	// problem 3. right primer need to be reversed - seems BLAST can handle it...

	cout<<"Output database..."<<endl;

	out_file.open(sDB.c_str());

	i = 0;

	for ( ORF_Iter = hORF.begin(); ORF_Iter != hORF.end(); ORF_Iter++ ){
		if(ORF_Iter->second.isActive==true && ORF_Iter->second.isPCR == true) {
			
			//output the sequence in FASTA format
			
			sName = ORF_Iter->first;
			cout<<"output #"<<i++<<" "<<sName<<endl;

			hORF[sName].getSequence();

			if(!hORF[sName].db_sequence.empty()){
				out_file<<hORF[sName].sTitle<<"_Forward"<<endl;
				out_file<<hORF[sName].db_sequence<<endl;
			}

			if(!hORF[sName].db_sequence_reverse.empty()){
				out_file<<hORF[sName].sTitle<<"_Reverse"<<endl;
				out_file<<hORF[sName].db_sequence_reverse<<endl;
			}

		}
	}

	out_file.close();
	out_file.clear();

	cout<<"Build sequencing BLAST database DONE!"<<endl;

	return;
}


// read seq001- in folder trimmed.seq
void build_query(){

	string sFileName, sFolderName;
	string SEQ, POS, sGeneName, sID;
	string s1, s2, s3, s4;
	vector<string>vLine;
	ifstream seq_file; 
	int i;

	sFolderName = "trimmed.seq\\";

	in_file.open(sSequencing.c_str());
	if(!in_file.good()){
		cout<<"ERROR: sequencing index file can't open!"<<endl;
		return;
	}

	while(in_file.good()){

		in_file>>SEQ>>POS>>sID;

		// forward seq
		sGeneName = sID+"_F";
		sFileName = sFolderName+SEQ+"_"+POS+"_F.trimmed.seq"; //"SEQ001_A01_F.trimmed.seq"
		
		seq_file.open(sFileName.c_str());
		
		if(!seq_file.good()){
			cout<<"ERROR: sequencing file "<<sFileName<<" can't open!"<<endl;
		}else{

			// read sequence
			s2.clear();

			while(seq_file.good()){

				getline(seq_file, s1, '\n');
				if(s1[0]=='!' || s1[0]=='#') continue;

				vLine.clear();
				SplitString(s1.c_str(), " ", vLine); // use space to saperate columns

				if(vLine.size()>1){  // title line

					if(!s2.empty()) {
						hCLONE[sGeneName].sequence = s2;
						s2.clear();
					}

					index = sGeneName.find_first_of ( "." );
					s3 = sGeneName.substr(0, index);
					s4 = sGeneName.substr(index+1, 1);
					hCLONE[sGeneName].sGeneName = s3;
					hCLONE[sGeneName].nSuffix = atoi(s4.c_str());
					hCLONE[sGeneName].isActive = true;
					hCLONE[sGeneName].sTitle = ">"+sGeneName+" "+s1;
				}else{
					s2 = s2 + s1;
				}
			}// end while file

			// last line
			if(!s2.empty() && hCLONE[sGeneName].isActive == true) {
				hCLONE[sGeneName].sequence = s2;
				s2.clear();
			}

			// find AAA AAA
			index = hCLONE[sGeneName].sequence.find( "AAAAAA" );
			
			if(index==string::npos) {
				cout<<"ERROR: cannot find AAA AAA in sequence "<<sGeneName<<" in file "<<sFileName<<endl;
			}else{
				hCLONE[sGeneName].sequence.erase(0, index); // remove any seq before AAA AAA
			}
		} // end else

		seq_file.close();
		seq_file.clear();

		// backward seq
		sGeneName = sID+"_R";
		sFileName = sFolderName+SEQ+"_"+POS+"_R.trimmed.seq"; //"SEQ001_A01_F.trimmed.seq"
		
		seq_file.open(sFileName.c_str());
		if(!seq_file.good()){
			cout<<"ERROR: sequencing file "<<sFileName<<" can't open!"<<endl;
			seq_file.close();
			seq_file.clear();	
			continue;
		}

		// read sequence
		s2.clear();

		while(seq_file.good()){

			getline(seq_file, s1, '\n');
			if(s1[0]=='!' || s1[0]=='#') continue;

			vLine.clear();
			SplitString(s1.c_str(), " ", vLine); // use space to saperate columns

			if(vLine.size()>1){  // title line

				if(!s2.empty()) {
					hCLONE[sGeneName].sequence = s2;
					s2.clear();
				}

				index = sGeneName.find_first_of ( "." );
				s3 = sGeneName.substr(0, index);
				s4 = sGeneName.substr(index+1, 1);
				hCLONE[sGeneName].sGeneName = s3;
				hCLONE[sGeneName].nSuffix = atoi(s4.c_str());
				hCLONE[sGeneName].isActive = true;
				hCLONE[sGeneName].sTitle = ">"+sGeneName+" "+s1;
			}else{
				s2 = s2 + s1;
			}
		}// end while file

		// last line
		if(!s2.empty() && hCLONE[sGeneName].isActive == true) {
			hCLONE[sGeneName].sequence = s2;
			s2.clear();
		}

		// // find TTT GTA
		index = hCLONE[sGeneName].sequence.find( "TTTGTA" );
		
		if(index==string::npos) {
			cout<<"ERROR: cannot find TTT GTA in sequence "<<sGeneName<<" in file "<<sFileName<<endl;
		}else{
			hCLONE[sGeneName].sequence.erase(0, index); // remove any seq before AAA AAA
		}
		
		seq_file.close();
		seq_file.clear();	

	} // end while	

	in_file.close();
	in_file.clear();

	// save to output file
	out_file.open(sQuery.c_str());

	i = 0;

	for ( ORF_Iter = hCLONE.begin(); ORF_Iter != hCLONE.end(); ORF_Iter++ ){
		if(ORF_Iter->second.isActive==true) {
			
			//output the sequences in FASTA format
			
			sGeneName = ORF_Iter->first;
			//cout<<"output #"<<i++<<" "<<sGeneName<<endl;

			out_file<<hCLONE[sGeneName].sTitle<<endl;
			out_file<<hCLONE[sGeneName].sequence<<endl;
		}
	}

	out_file.close();
	out_file.clear();

	return;
}

void call_blast(){

	string sGeneName;
	string s1, s2, s3;
	char buffer[20];
	ofstream bat_file;
	string sDIR;

	build_database();  
	
	build_query();  

	sDIR = "blast_result\\";
	s1 = sDIR+"run_blast.bat";
	bat_file.open(s1.c_str());

	for ( ORF_Iter = hCLONE.begin(); ORF_Iter != hCLONE.end(); ORF_Iter++ ){

		// build query file
		if(ORF_Iter->second.isActive==true) {
			sGeneName = ORF_Iter->first;
			out_file.open((sDIR+sGeneName).c_str());
			out_file<<hCLONE[sGeneName].sTitle<<endl;
			out_file<<hCLONE[sGeneName].sequence<<endl;
			out_file.close();
			out_file.clear();
		}

		// build db file
		s1 = sGeneName; // e.g. "AT3G13050_R"
		index =s1.find( "_" );
		
		if(index==string::npos) {
			cout<<"ERROR:  format error in name "<<sGeneName<<endl;
		}else{
			s1.erase(index, s1.size()-1); 
		}

		for(int k=1; k<=SUFFIX_MAX; k++){

			_itoa( k, buffer, 10 );
			s2 = s1 + "."+buffer;

			if(hORF[s2].isActive == true){
				s3 = s2 + "_Forward";
				out_file.open((sDIR+s3).c_str());
				out_file<<hORF[s2].sTitle<<"_Forward"<<endl;
				out_file<<hORF[s2].db_sequence<<endl;
				out_file.close();
				out_file.clear();

				// write a bat file
				s3 = "bl2seq -i "+sGeneName +" -j "+s3+" -p blastn -F F > " + sGeneName +"_"+ s3+".result";
				bat_file<<s3<<endl;

				s3 = s2 + "_Reverse";
				out_file.open((sDIR+s3).c_str());
				out_file<<hORF[s2].sTitle<<"_Reverse"<<endl;
				out_file<<hORF[s2].db_sequence_reverse<<endl;
				out_file.close();
				out_file.clear();

				// write a bat file
				s3 = "bl2seq -i "+sGeneName +" -j "+s3+" -p blastn -F F> " + sGeneName +"_"+ s3+".result";
				bat_file<<s3<<endl;

			}
		
		} // end for

	}// end for each query

	bat_file.close();
	bat_file.clear();

	return;
}

void analyze_blast(){

	int i, j;
	string sDIR = "parsing\\";
	string sID, sHit;
	string s1, s2, s3, s4, sFileName, sIndexName;
	vector<string> vLine;

	// build hCLONE
	build_query();

	//read index.txt, which has 6 columns 
	//#Query	#Sbjct	#Orientation 1	#File 1	#Orientation 2	#File 2
	// blast result saved in file <FileName>
	sIndexName = "index.txt";

	in_file.open((sDIR+sIndexName).c_str());
	out_file.open(sParsing.c_str());
	
	//title 
	out_file<<"ID	Orientation	Length	Identities	e_value	big_insertion/deletion	no_aa_change	substitution_aa_change	num_of_substitution_aa_change	stop_codon_change	no_start_codon_change	";
	out_file<<"normal_insertion	frame_shift_insertion	normal_deletion	frame_shift_deletion"<<endl;

	if(!in_file.good()){
		cout<<"ERROR: Cannot find the index file!"<<endl;
		return;
	}

	i = 1;
	while(in_file.good()){
		//hHIT
		getline(in_file, s1, '\n');
		if(s1[0]=='!' || s1[0]=='#') continue;

		vLine.clear();
		SplitString(s1.c_str(), "\t", vLine); // use space to saperate columns
		
		if(vLine.size()<6) continue; //non-effective line			
		
		//if(i%MAX_JMP==0) 
			cout<<"parsing gene #"<<i<<"..."<<endl;
		
		sID = vLine[0];  // e.g. AT1G06080_R
		sHit = vLine[1];
		if(hCLONE[sID].isActive!=true) {
			continue;
		}

		//parsing forward file
		i++;
		sFileName = vLine[3];
		j = parsing(sID, sDIR+sFileName);

		if(j!=1){
			//parsing reverse file
			sFileName = vLine[5];
			j = parsing(sID, sDIR+sFileName);
			if(j!=1) {
				cout<<"ERROR: Cannot find matched sequences for "<<sID<<"!"<<endl;
				out_file<<sID<<"	NO_MATCH"<<endl;
				continue;
			}			
		}

		//save to out_file
			out_file<<sID<<"	"<<hCLONE[sID].cOrientation<<"	"<<hCLONE[sID].sequence.size()<<"	";
			out_file<<hCLONE[sID].sIdentities<<"	"<<hCLONE[sID].e_value<<"	";
			// big_insertion/deletion
			if(abs(hCLONE[sID].sequence.size()-hCLONE[sID].nHSP_length<BIG_DELETION_INSERTION)) out_file<<"NO	";
			else out_file<<"YES	";
			//replacment
			out_file<<hCLONE[sID].no_change<<"	"<<hCLONE[sID].normal_change<<"	"<<hCLONE[sID].num_normal_change<<"	";
			out_file<<hCLONE[sID].stop_codon_change<<"	"<<hCLONE[sID].no_start_codon_change<<"	";
			//inserstion
			out_file<<hCLONE[sID].normal_insertion<<"	"<<hCLONE[sID].frame_shift_insertion<<"	";
			//deletion
			out_file<<hCLONE[sID].normal_deletion<<"	"<<hCLONE[sID].frame_shift_deletion<<endl;

	}

	in_file.close();
	in_file.clear();

	out_file.close();
	out_file.clear();

	return;
}

///////////////////////////////////////////////////////////////////////////
// kernal function: parsing blast results
// return 1 if hCLONE[sID] matches sFileName, otherwise 0
//////////////////////////////////////////////////////////////////////////
int parsing(string sID, string sFileName){

	ifstream ifParsing;
	int nScore, nReturn, nStartPos, nHSP_length;
	size_t j, k;
	char cOrientation;
	string sIdentities, sQueryChange, sSbjctChange;
	string sQueryProtein, sSbjctProtein, s1, s2, s3, s_k;
	string sSbjct;
	bool bExit;
	basic_string <char>::size_type index1, index2;
	vector<string> vParsing;
	char buffer[200];
	double e_value;

	//init
	nReturn = 0;
	nStartPos = -1;
	nScore = -1;
	bExit = false;

	// file operation
	ifParsing.open(sFileName.c_str());
	if(!ifParsing.good()){
		cout<<"ERROR: File "<<sFileName<<" can't open!"<<endl;
		return nReturn;
	}

	while(ifParsing.good() && bExit==false ){
		getline(ifParsing, s1, '\n');
		str_trim(s1);
		if(s1.empty()) continue;  //skip empty lines
		vParsing.clear();
		SplitString(s1.c_str(), " ", vParsing);

		if(vParsing[0]=="Score") {
			if(nScore!=-1) bExit = true;  // skip the 2nd, 3rd... matches
			else {
				nScore = atoi(vParsing[2].c_str());
				e_value = atof(vParsing[vParsing.size()-1].c_str());
			}
		}

		if(vParsing[0]=="Identities") {
			sIdentities = vParsing[2];
			
			index1 = sIdentities.find_first_of ( "/" );
			s3 = sIdentities.substr(0, index1); //sIdentities.substr(index1+1,sIdentities.size()-index1);
			nHSP_length = atoi(s3.c_str());
		}
		if(vParsing[0]=="Strand" ){
			if(vParsing[2] == vParsing[4]) cOrientation = 'F';
			else cOrientation = 'R';
		}

		//start details...
		if(vParsing[0]=="Query:"){
			
			if(nStartPos == -1){  // first line
				nStartPos = atoi(vParsing[1].c_str());

				if(nStartPos!=1) { // not start from AAA AAA or TTT GTA
					bExit = true; // skip this file
					continue;
				}else{
					nReturn = 1;  // this is the maching file
					hCLONE[sID].cOrientation = cOrientation;
					hCLONE[sID].e_value = e_value;
					hCLONE[sID].sIdentities = sIdentities;
					hCLONE[sID].nHSP_length = nHSP_length;
				}
			}else{ 
				nStartPos = atoi(vParsing[1].c_str());
			}

			sQuery = vParsing[2];

			// read the following 2 lines
			getline(ifParsing, s2, '\n');  //||||||||||||||  ||||
			str_trim(s2);
			getline(ifParsing, s3, '\n'); // sbjct
			vParsing.clear();
			SplitString(s3.c_str(), " ", vParsing);
			sSbjct = vParsing[2];

			//is there any mismatch?
			if(s2.find(" ")==string::npos) continue;

			//find mismatches
			for(j=0; j<s2.size(); j=j+3){

				if(j+1>=s2.size() || j+2>=s2.size()) break; //not in pair
			
				if(s2[j]==' ' || s2[j+1]==' ' || s2[j+2]==' '){
					
					//p = (int)(j*1.0/3.0)*3; q = j%3;

					sQueryChange.clear();
					sSbjctChange.clear();
					for(k=j; k<j+3; k++){ // 3 nuclear acids = 1 protein
						sQueryChange = sQueryChange+sQuery[k]; 
						sSbjctChange = sSbjctChange+sSbjct[k];
					}
				
					sQueryProtein = getProteinName(sQueryChange);
					sSbjctProtein = getProteinName(sSbjctChange);

					index1 = sQueryChange.find("-");  // find blank
					index2 = sSbjctChange.find("-");

					k = j+nStartPos;  // position
					_itoa( (int)k, buffer, 10 );
					s_k = buffer;

					if(index1==string::npos && index2==string::npos) { //replacement
								
						if(sQueryProtein == sSbjctProtein) {
							//cout<<"Replaced without amino acid change"<<endl;
							hCLONE[sID].no_change = hCLONE[sID].no_change + \
								"[R "+s_k +", " +sSbjctChange + " -> " +sQueryChange + "] ";
						}else{
							if(sQueryProtein == "stop") {
								//cout<<"Replaced with a stop codon"<<endl;
								hCLONE[sID].stop_codon_change = hCLONE[sID].stop_codon_change + \
									"[R "+s_k +", " +sSbjctChange + " -> " +sQueryChange + ", " \
									+ sSbjctProtein + " => " + sQueryProtein+ "] ";
							}else if (sSbjctProtein == "start" && k<60) { // 27 vector+15 max utr < 60
								//cout<<"Replaced a start codon"<<endl;
								hCLONE[sID].no_start_codon_change = hCLONE[sID].no_start_codon_change + \
									"[R "+s_k +", " +sSbjctChange + " -> " +sQueryChange + ", " \
									+ sSbjctProtein + " => " + sQueryProtein+ "] ";
							}else {
								//cout<<"Replaced with normal amino acid change"<<endl;
								hCLONE[sID].normal_change = hCLONE[sID].normal_change + \
									"[R "+s_k +", " +sSbjctChange + " -> " +sQueryChange + ", " \
									+ sSbjctProtein + " => " + sQueryProtein+ "] ";
								hCLONE[sID].num_normal_change++;
							}
						}
					}

					if(index1!=string::npos && index2==string::npos) { 	//deletion	
						if(sQueryChange=="---"){
							//cout<<"Deletion without shift"<<endl;
							hCLONE[sID].normal_deletion = hCLONE[sID].normal_deletion + \
								"[D "+s_k+ ", 3] ";
						}else {
							int i = 0; 
							for(int n=0; n<3; n++) if(sQueryChange[n]=='-') i++;
							_itoa( i, buffer, 10 );
							//cout<<"Deletion with shift"<<endl;
							hCLONE[sID].frame_shift_deletion = hCLONE[sID].frame_shift_deletion + \
								"[D "+s_k+ ", "+buffer + "] ";
						}
					}
					
					if(index1==string::npos && index2!=string::npos) { 	//insertion
						if(sSbjctChange=="---"){
							//cout<<"Insertion without shift"<<endl;
							hCLONE[sID].normal_insertion = hCLONE[sID].normal_insertion + \
								"[I "+s_k+ "3] ";
						}else {
							int i = 0; 
							for(int n=0; n<3; n++) if(sSbjctChange[n]=='-') i++;
							_itoa( i, buffer, 10 );
							//cout<<"Insertion with shift"<<endl;
							hCLONE[sID].frame_shift_insertion = hCLONE[sID].frame_shift_insertion +\
								"[I "+s_k+ ", "+buffer + "] ";
						}
					}
					if(index1!=string::npos && index2!=string::npos){ // insertion and deletion
						cout<<"Insersion and deletion"; //not happened in BIONEER192
					}
				}
			}// end for mismatches j
		}// end if details
	} // end while file

	ifParsing.close();
	ifParsing.clear();

	return nReturn;
}

string getProteinName(string sDNA){

	string s;

	s = "NO_MATCH";
	std::transform( sDNA.begin(), sDNA.end(), sDNA.begin(), ::toupper );

	if(sDNA=="TTT") s= "Phe";
	if(sDNA=="TTC") s= "Phe";
	if(sDNA=="TTA") s= "Leu";
	if(sDNA=="TTG") s= "Leu";
	if(sDNA=="CTT") s= "Leu";
	if(sDNA=="CTC") s= "Leu";
	if(sDNA=="CTA") s= "Leu";
	if(sDNA=="CTG") s= "Leu";
	if(sDNA=="ATT") s= "Ile";
	if(sDNA=="ATC") s= "Ile";
	if(sDNA=="ATA") s= "Ile";
	if(sDNA=="ATG") s= "start"; //Met
	if(sDNA=="GTT") s= "Val";
	if(sDNA=="GTC") s= "Val";
	if(sDNA=="GTA") s= "Val";
	if(sDNA=="GTG") s= "Val";
	if(sDNA=="TCT") s= "Ser";
	if(sDNA=="TCC") s= "Ser";
	if(sDNA=="TCA") s= "Ser";
	if(sDNA=="TCG") s= "Ser";
	if(sDNA=="CCT") s= "Pro";
	if(sDNA=="CCC") s= "Pro";
	if(sDNA=="CCA") s= "Pro";
	if(sDNA=="CCG") s= "Pro";
	if(sDNA=="ACT") s= "Thr";
	if(sDNA=="ACC") s= "Thr";
	if(sDNA=="ACG") s= "Thr";
	if(sDNA=="ACA") s= "Thr";
	if(sDNA=="GCT") s= "Ala";
	if(sDNA=="GCC") s= "Ala";
	if(sDNA=="GCA") s= "Ala";
	if(sDNA=="GCG") s= "Ala";
	if(sDNA=="TAT") s= "Tyr";
	if(sDNA=="TAC") s= "Tyr";
	if(sDNA=="TAA") s= "stop";
	if(sDNA=="TAG") s= "stop";
	if(sDNA=="CAT") s= "His";
	if(sDNA=="CAC") s= "His";
	if(sDNA=="CAA") s= "Gln";
	if(sDNA=="CAG") s= "Gln";
	if(sDNA=="AAT") s= "Asn";
	if(sDNA=="AAC") s= "Asn";
	if(sDNA=="AAA") s= "Lys";
	if(sDNA=="AAG") s= "Lys";
	if(sDNA=="GAT") s= "Asp";
	if(sDNA=="GAC") s= "Asp";
	if(sDNA=="GAA") s= "Glu";
	if(sDNA=="GAG") s= "Glu";
	if(sDNA=="TGT") s= "Cys";
	if(sDNA=="TGC") s= "Cys";
	if(sDNA=="TGA") s= "stop";
	if(sDNA=="TGG") s= "Trp";
	if(sDNA=="CGT") s= "Arg";
	if(sDNA=="CGC") s= "Arg";
	if(sDNA=="CGA") s= "Arg";
	if(sDNA=="CGG") s= "Arg";
	if(sDNA=="AGT") s= "Ser";
	if(sDNA=="AGC") s= "Ser";
	if(sDNA=="AGA") s= "Arg";
	if(sDNA=="AGG") s= "Arg";
	if(sDNA=="GGT") s= "Gly";
	if(sDNA=="GGC") s= "Gly";
	if(sDNA=="GGA") s= "Gly";
	if(sDNA=="GGG") s= "Gly";

	return s;
}

char * _strupr_new (char * string){
    char *cp;       // traverses string for C locale conversion
    
	for (cp = string; *cp; ++cp){
            if ('a' <= *cp && *cp <= 'z')
                *cp += 'A' - 'a';
    }

    return(string);
}