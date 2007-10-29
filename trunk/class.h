

ORF::ORF(){

	this->sequence.clear();
	this->upstream.clear();
	this->downstream.clear();
	this->left_primer.clear();
	this->right_primer.clear();
	this->db_sequence.clear();
	this->db_sequence_reverse.clear();

	this->isActive = false;
	this->nSuffix = 0;
	this->nLeftExtension = 0;
	this->nRightExtension = 0;
	this->isPCR = false;

	//for blast
	this->no_change.clear(); this->normal_change.clear();
	this->stop_codon_change.clear(); this->no_start_codon_change.clear();
	this->normal_insertion.clear(); this->frame_shift_insertion.clear();
	this->normal_deletion.clear(); this->frame_shift_deletion.clear();
	this->num_normal_change = 0;

	return;
}

ORF::~ORF(){

	this->sequence.clear();
	this->upstream.clear();
	this->downstream.clear();
	this->left_primer.clear();
	this->right_primer.clear();
	this->db_sequence.clear();
	this->db_sequence_reverse.clear();

	this->isActive = false;
	this->nSuffix = 0;
	this->nLeftExtension = 0;
	this->nRightExtension = 0;
	this->isPCR = false;

	//for blast
	this->no_change.clear(); this->normal_change.clear();
	this->stop_codon_change.clear(); this->no_start_codon_change.clear();
	this->normal_insertion.clear(); this->frame_shift_insertion.clear();
	this->normal_deletion.clear(); this->frame_shift_deletion.clear();

	return;
}


// form sequence form cds, utr and primer, output to db_sequence
void ORF::getSequence(){

	string str, s1, s2, s3, sName;
	size_t i, j, len;
	int start, end;
	char up[20];

	start = this->nLeftExtension;
	end = this->nRightExtension;
	sName = this->sGeneName;

	s1.clear();
	s2.clear();

	// for upstream
	str = this->upstream;

	if(start<0 && str.size()<(0-start)) {  // abnormal return
		this->db_sequence.clear();
		return;
	}

	if(start<0){
		len = str.size(); 
		j = 0;
		for(i = len + start;i<=len;i++)	up[j++] = str[i];
		s1 = up;
		s1 = ORF::replaceStartStopCodon(s1, 5); 
	}

	// for downstream
	str = this->downstream;

	if(end>0 && str.size()<end) {  // abnormal return
		this->db_sequence.clear();
		return;
	}

	if(end>0)	{
		s2 = str.substr(0, end);
		s2 = ORF::replaceStartStopCodon(s2, 3); 
	}

	// for sequence
	str = this->sequence;
	if(end==-3){
		s3 = str.substr(0, str.size()-3);
	}else{//end = 3,6,9,12
		s3 = this->sequence;
		s3 = ORF::replaceStartStopCodon(s3, 0);
	}

	// prepare both sequences
//	if (this->isForward == true){
//		this->db_sequence = sFix5+ "GGG" + s1 + s3 + s2 + "CCC" + sFix3;
//	}else{
//		this->db_sequence = sFix3_R+ "GGG" + s1 + s3 + s2 + "CCC" + sFix5_R;
//	}

	this->db_sequence = sFix5+ "GGG" + s1 + s3 + s2 + "CCC" + sFix3;
	this->db_sequence_reverse = sFix3_R+ "GGG" + s1 + s3 + s2 + "CCC" + sFix5_R;

	return;
}

// if flag ==1, find stop codon in the full string, otherwise only find stop codon at the last 3 amino acides
string ORF::replaceStartStopCodon(string seq, int flag){	
	size_t i, n;
	string str1, str2; // use STL here
	basic_string <char>::size_type index;
	string cs; 
	bool stop;
	char *str, *p;

	str = (char*) malloc(sizeof(char)*seq.size());

	n = 0;
	i = seq.size();
	strcpy(str, seq.c_str());

	if(strlen(str)<3) return seq;

	if(flag ==0){ //only find stop codon at the last 3 amino acides
	
		p = &str[i-3];

		if(strcmp(p, "TGA")==0) {str[i-3]='T';str[i-2]='C';str[i-1]='A';}
		else if(strcmp(p, "TAA")==0) {str[i-3]='T';str[i-2]='T';str[i-1]='A';}
		else if(strcmp(p, "TAG")==0) {str[i-3]='T';str[i-2]='T';str[i-1]='G';}

		n=1;

	}else if(flag ==5){// find stop codons in the full string, must in-frame
		stop = false;
		str1 = str;
		while(stop==false){

			stop = true;

			str2 = "TGA";
			cs="TCA";
			index = str1.find (str2,0);
			if (index != string::npos && index%3==0 ){
				str1.replace(index, 3, cs);
				stop = false;
				n++;
			}
				
			str2 = "TAG";
			cs="TTG";
			index = str1.find (str2,0);
			if (index != string::npos && index%3==0 ){
				str1.replace(index, 3, cs);
				stop = false;
				n++;
			}

			str2 = "TAA";
			cs="TTA";
			index = str1.find (str2,0);
			if (index != string::npos && index%3==0 ){
				str1.replace(index, 3, cs);
				stop = false;
				n++;
			}
		}// end while

		str2 = "ATG"; //start codon; any position in 5' UTR
		cs="TTG";
		index = str1.find (str2,0);
		if (index != string::npos ){
			str1.replace(index, 3, cs);
			n=n+1; //1000
		}

		for(i = 0; i<str1.length(); i++)	str[i]= str1[i];

	}else  if(flag ==3){
		stop = false;
		str1 = str;
		while(stop==false){

			stop = true;

			str2 = "TGA";
			cs="TCA";
			index = str1.find (str2,0);
			if (index != string::npos && index%3==0  ){
				str1.replace(index, 3, cs);
				stop = false;
				n++;
			}
				
			str2 = "TAG";
			cs="TTG";
			index = str1.find (str2,0);
			if (index != string::npos && index%3==0 ){
				str1.replace(index, 3, cs);
				stop = false;
				n++;
			}

			str2 = "TAA";
			cs="TTA";
			index = str1.find (str2,0);
			if (index != string::npos && index%3==0 ){
				str1.replace(index, 3, cs);
				stop = false;
				n++;
			}
		}// end while

		for(i = 0; i<str1.length(); i++)	str[i]= str1[i];
	} //end else if

	cs = str;

	return cs;
}

