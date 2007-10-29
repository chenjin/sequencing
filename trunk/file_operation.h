template< typename StrT >
int SplitString(const char* str, const char* delim, vector<StrT>& results){
  char* pstr = const_cast<char*>(str);
  char* r = NULL;
  r = strstr(pstr, delim);
  int dlen =(int) strlen(delim);
  bool empties = false;  // no empty entity
  
  while( r != NULL ) {
    
	char* cp = new char[(r-pstr)+1];
    memcpy(cp, pstr, (r-pstr));
    cp[(r-pstr)] = '\0';

	//trim cp

    if( strlen(cp) > 0 || empties ) {
      StrT s(cp);
      results.push_back(s);
    }

    delete[] cp;
    pstr = r + dlen;
    r = strstr(pstr, delim);
  } //end while

  if( strlen(pstr) > 0 || empties )   results.push_back(StrT(pstr));

  return (int) results.size();

}

void str_trim(string& str){

  string::size_type pos = str.find_last_not_of(' ');

  if(pos != string::npos) {

    str.erase(pos + 1);
    pos = str.find_first_not_of(' ');
    if(pos != string::npos) str.erase(0, pos);

  }else {
	  str.erase(str.begin(), str.end());
  }

  return;
}