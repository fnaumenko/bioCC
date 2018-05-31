#include "common.h"
#include <sstream>
#ifdef OS_Windows
	#include <algorithm>
	#define SLASH '\\'		// standard Windows path separator
	#define REAL_SLASH '/'	// is permitted in Windows too
#else
	#define SLASH '/'	// standard Linux path separator
#endif


/************************ common Functions ************************/

// Gets number of digist in a integral value
//	@val: integral value
//	@isLocale: if true then adds number of '1000' separators
//	return: number of digist without minus symbol or 0 if value is 0
BYTE DigitsCount (LLONG val, bool isLocale)
{
	BYTE res = 0;
	for(; val; val/=10, res++);
	if(isLocale)	res += (res-1)/3;
	return res;
}

// Returns string represents the percent of part relatively total
//	@percent: value of percent
//	@precision: count of mapped digits; 
//	if count of value's mapped digits is more then that (too little percent), printed "<n%"
//	or exactly by default
//	@fieldWith: the width of the display field insine parentheses or exactly by default
//	@parentheses: if true parenthesize the value (default)
string	sPercent(float percent, BYTE precision, BYTE fieldWith, bool parentheses)
{
	float threshold = (float)pow(10.0, -precision);
	stringstream sout;
	if(parentheses)		sout << " (";
	if( percent && percent < threshold )
		sout << '<' << threshold;
	else {
		if(precision) {
			if(percent >= 100)	precision = 3;
			sout << setprecision(precision);
		}
		sout << setw(fieldWith) <<  percent;
	}
	sout << PERS;
	if(parentheses)		sout << ')';
	return sout.str();
}

// Prints horizontal line
//	@w: width of line
void PrintHorLine(int w)
{
#ifdef OS_Windows
	wcout << setw(w) << setfill(L'\304') << L'\304' << L'\n';
#else
	for(int i=0; i<w; i++)	cout << "─";	cout << EOL;
	//cout << setw(w) << setfill(HPH) << HPH << EOL;
#endif
}

#if defined _WIGREG || defined _BIOCC

// Align position to the up or down resoluation level
// f.e. by resoluation==5 pos 102 -> 100+relative, pos 104 -> 105++relative
//	@pos: chromosome's position
//	@res: resoluation
//	@relative: 0 or 1
//	1 used for 1-relative position (the first base is 1, WIG)
//	0 used for 0-relative position (BED)
//	return: aligned position
chrlen AlignPos(chrlen pos, BYTE res, BYTE relative)
{
	short rest = pos % res;
	return rest - relative ?
		pos + relative - rest + (rest > res<<1 ? res : 0) :
		pos;
}

#endif

//size_t getAvailSystemMemory()
//{
//#ifdef OS_Windows
//	MEMORYSTATUSEX status;
//	status.dwLength = sizeof(status);
//	GlobalMemoryStatusEx(&status);
//	return (size_t)status.ullAvailPhys;
//#else
//	long pages = sysconf(_SC_PHYS_PAGES);
//	long page_size = sysconf(_SC_PAGE_SIZE);
//	return pages * page_size;
//#endif
//}

//#ifdef OS_Windows
//string Wchar_tToString(const wchar_t* wchar)
//{
//	string str = strEmpty;
//	while(*wchar)	str += char(*wchar++);
//	return str;
//}
//
//wchar_t* StringToWchar_t(const string &str, wchar_t* wchar)
//{
//	int i = (int)str.size();
//	wchar[i] = 0;
//	for(--i; i >= 0; i--)
//		wchar[i] = (wchar_t)str[i];
//	return wchar;
//}
//#endif

/************************ end of common Functions ************************/

const char* Gr::title[] = {"FG","BG"};

/************************ class Options ************************/
#define ENUM_REPLACE '?'	// symbol in description that is replaced by enum value

// prints value in parentheses
#define PRINT_IN_PRTHS(v)	cout<<" ["<<(v)<<']'

const char*	OptTitle = "Option ";
const char*	optTitle = " option ";
const char* Ambiguous = "Ambiguous";
const char*	Default = " Default: ";
const char*	Missing = "missing ";


const char* Options::Booleans [] = {"OFF", "ON"};

const char* Options::_TypeNames [] = {
	NULL, "<name>", "<char>", "<int>", "<float>", "<long>", NULL, NULL,  "<int;int>", NULL, NULL
};

const char Options::Option::EnumDelims [] = {'|', ',', ';'};

const string sValue = "value";
//const string sWrongValue = "wrong " + sValue;

// Checks is string represents digital value
bool isValidFloat(const char *str)
{
	char c = *str;
	BYTE dotCnt = 0, eCnt = 0;

	// str==NULL is checked before
	if(!isdigit(c))						// check first char
		if(c == DOT)	dotCnt++;
		else if(c != HPH && c != PLUS)	return false;
	for(str++; *str; str++)				// check next chars
		if((c=*str) == DOT)	
			if(dotCnt)	return false;	// more than one dot
			else		dotCnt++;
		else if(tolower(c) == 'e')
			if(eCnt)	return false;	// more than one 'e'
			else		eCnt++;
		else if(!isdigit(c))	return false;

	return true;
}

// Recursively prints string with replaced ENUM_REPLACE symbol by enum/combi value.
//	@buff: external buffer to copy and output temporary string
//	@vals: enum/combi values or NULL for other types
//	@cnt: external counter of enum/combi values
void PrintTransformDescr(char* buff, const char** vals, short* cnt)
{
	if(vals) {			// enum/combi?
		const char* subStr = strchr(buff, ENUM_REPLACE);
		if(subStr) {	// something to replace by enum value
			size_t strLen = subStr - buff;
			buff[strLen] = 0;
			cout << buff << vals[(*cnt)++];	// output substring and enum value
			buff[strLen] = ENUM_REPLACE;
			PrintTransformDescr(buff+strLen+1, vals, cnt);
		}
		else	cout << buff;
	}
	else	cout << buff;
}

// Recursively prints string with EOL inside as a set of left-margin strings
// Used to output aligned option descriptions
// First string is printed from current stdout position.
// Last substring doesn't include EOL.
//	@buff: external buffer to copy and output temporary string
//	@str: input string with possible EOLs
//	@subStr: substring of input string to the first EOL, or NULL if input string is not ended by EOL
//	@vals: enum/combi values or NULL for other types
//	@cnt: external counter of enum/combi values
void PrintSubLine(char* buff, const char* str, const char* subStr, const char** vals, short* cnt)
{
	if(subStr) {	// is substring ended by EOL exist?
		// form substring
		size_t strLen = subStr - str;
		strncpy(buff, str, strLen);
		buff[strLen] = 0;
		PrintTransformDescr(buff, vals, cnt);	// output enum values
		cout << EOL;
		for(BYTE t=0; t<OPT_DESCF_TSHIFT; t++)	cout << TAB;
		str = subStr + 1;		// skip EOL
		subStr = strchr(str, EOL);
		PrintSubLine(buff, str, subStr, vals, cnt);
	}
	else {			// output rest of initial string without EOL
		strcpy(buff, str);
		PrintTransformDescr(buff, vals, cnt);
	}
}



// Sets option value.
//	@isWord: true if option is a word, false if option is a char
//	@opt: option without HYPHENs
//	@val: value of option
//	@isNextOpt: true if next parameter is option
//	Return: 0 if success, -1 if not found, 1 if option or value is wrong
int Options::Option::SetVal(bool isWord, char* opt, char* val, bool isNextOpt)
{
	if( isWord ) { if(Str==NULL || strcmp(Str, opt))	return -1; }
	else if(Char != *opt)	return -1;

	if(Sign.Is(Signs::Trimmed))	return PrintAmbigOpt(isWord, opt, "Duplicated", NULL);
	// check value existence
	bool noVal = val==NULL || (val[0]==HPH && !isdigit(val[1]));	// true if no value, including negative
	if(ValRequired())	{ if(noVal)	return PrintWrongOpt(isWord, opt, NULL, sValue + " required"); }
	else if(!noVal&&isNextOpt)	return PrintWrongOpt(isWord, opt, NULL, sValue + " prohibited");

	Sign.MarkAs(Signs::Trimmed);
	switch(ValType) {
		case tCHAR:	if(NVal != vUNDEF)
				return strlen(val) > 1 ? 
					PrintWrongOpt(isWord, opt, val) :
					SetTriedFloat(*val, isWord);	// value is treated as int,
		case tNAME: SVal = val;			return 0;	// otherwise it is treated as string
		case tENUM:	return SetEnum(val) ? 0 : PrintWrongOpt(isWord, opt, val);
		case tCOMB:	return SetComb(val) ? 0 : PrintWrongOpt(isWord, opt, val);
		case tPAIR:	return !SetPair(val);
		case tHELP:	return PrintUsage(true);
		case tVERS:	return PrintVersion();
		default:
			if(isValidFloat(val))
				return SetTriedFloat(atof(val), isWord);	// numerical value
	}
	return PrintWrongOpt(isWord, opt, val);
}

// Check option for obligatory.
//	return: -1 if option is obligatory but not stated, otherwise 1
int Options::Option::CheckOblig()
{
	if( Sign.Is(Signs::Oblig) && ValRequired()
	&& ( (ValType == tNAME && SVal == NULL)
		|| (ValType != tNAME && NVal == vUNDEF) ) ) {
			cerr << Missing << "required option " 
					<< OptToStr(Char==HPH, Char==HPH ? Str : &Char) << EOL;
			return -1;
	}
	return 1;
}

// Returns string represented pair of value's separated by delimiter.
const string Options::Option::PairValsToStr(const pairVal* vals) 
{
	return static_cast<ostringstream & >( ostringstream() << dec 
		<< vals->first << EnumDelims[2] << vals->second ).str();
}

// Checks limits and set numerical value
//	@val: numerical value
//	@isWord: true if option enters by long name
//	return: 1 if limits are exceeded, otherwise 0
int Options::Option::SetTriedFloat(float val, bool isWord)
{
	int outOfLimit = -1;
	if( val < MinNVal )			outOfLimit = 0;
	else if( val > MaxNVal )	outOfLimit = 1;
	if( outOfLimit >= 0 ) {
		const char* sign[] = {"less", "more"};
		cerr << OptTitle << OptToStr(isWord, isWord ? Str : &Char);
		cerr << SepSCl << sValue << setprecision(4) << BLANK << val
				<< " is " << sign[outOfLimit] << " than permissible " 
				<< (!outOfLimit ? MinNVal: MaxNVal ) << EOL;
		return 1;
	}
	NVal = val;
	return 0;
}

// Checks and sets enum option value.
//	@val: input value as string
//	return: true if success
bool Options::Option::SetEnum(char* val)
{
	if( ValRequired() ) {				// non-boolean template
		int ind = GetEnumInd(val);
		if( ind < 0 )	return false;
		NVal = ind + MinNVal;	// add possible minimum value as a shift
	}
	else
		NVal = !NVal;			// invers boolean value
	return true;
}

// Checks and sets enum option value.
//	@val: input value as C string
//	return: true if success
bool Options::Option::SetComb(char* vals)
{
	char* pdelim;	// a pointer to the first occurrence of delimiter COMMA in vals
	int	ind;
	bool ret = true;

	NVal = 0;	// reset default value
	// run through given values
	for(char* val = vals; true; val = pdelim + 1) {
		pdelim = strchr(val, EnumDelims[1]);
		if(pdelim)	*pdelim = '\0';		// temporary cut 'val' to C string with single value
		ind = GetEnumInd(val);
		if(ind < 0)		ret = false;
		else// set bitwise val, in which running number of each bit (from right end)
			// corresponds to number of finded value in enum
			NVal = int(NVal) ^ (1<<ind);
		if(pdelim)	*pdelim = EnumDelims[1];	// restore 'vals' string
		else	break;							// no delimiter: last or single value
	}
	return ret;
}

// Checks and sets pair option value.
//	@val: input pair value as C string
//	return: true if success
bool Options::Option::SetPair(char* vals)
{
	char* pdelim = strchr(vals, EnumDelims[2]);	// a pointer to the delimiter ';' in vals

	if(!pdelim)	return !PrintWrongOpt(true, Str, vals);
	short val = atoi(vals);
	if(SetTriedFloat(val, true))	return false;
	((pairVal*)SVal)->first = val;
	val = atoi(pdelim+1);
	if(SetTriedFloat(val, true))	return false;
	((pairVal*)SVal)->second = val;
	return true;
}

inline void Options::Option::PrintCharDefValue() {
	PRINT_IN_PRTHS(char(NVal));
}

// Prints option in full or short way.
//	@descr: if true, prints in full way: signature, description (marks as Required if needed), default value,
//	otherwise signature only
void Options::Option::Print(bool descr)
{
	USHORT	len = 0;		// first len is used as counter of printed chars
	bool	fixValType = ValType==tENUM || ValType==tCOMB;
	char*	buffer;

	if(descr)	{ cout << BLANK;	len++;	}	// full way: double blank first
	// ** signature
	cout << BLANK << HPH << Char;	len += 3;
	if(Str)	{
		if(Char != HPH) {	cout << "|--";	len += 3; }
		cout << Str;
		len += strlen(Str);
	}
	// ** option values
	if(fixValType)		len += PrintEnumVals();		// print enum values
	else if((buffer = const_cast<char*>(_TypeNames[ValType])) != NULL) {	// print value type
		cout << BLANK << buffer;
		len += 1 + strlen(buffer);
	}
	// ** description
	if(!descr)	return;
	short cnt = OPT_DESCF_TSHIFT - len / 8;	// 3*8: right boundary of descriptions
	// align description 
	if(!cnt || (cnt<0 && !(len % 8)))	cnt = 1;
	for(BYTE i=0; i<cnt; i++)	cout << TAB;

	// print description
	cnt = 0;	// use as external enum counter
	len = strlen(Descr);	// from now len is used as length of description string
	buffer = new char[len+1];
	PrintSubLine(
		buffer,
		Descr,
		strchr(Descr, EOL), 
		fixValType ? (const char**)SVal : NULL, 
		&cnt
	);
	delete [] buffer;
	if(AddDescr) {
		if(Descr[len-1] != EOL)	cout << BLANK;
		cout << AddDescr;
	}
	if(Sign.Is(Signs::Oblig))	cout << " Required";
	else if(ValType >= tHELP)	cout << " and exit";

	// print default value
	if(ValRequired() && NVal != vUNDEF)
		switch(ValType) {
			case tENUM:
			case tCOMB:
				if(NVal >= MinNVal)	// do not print default if it is never set by user
					PRINT_IN_PRTHS(	((char**)SVal)
						[int(NVal)-int(MinNVal)] );	// offset by min enum value
				break;
			case tPAIR:	PRINT_IN_PRTHS(PairValsToStr((pairVal*)SVal)); break;
			case tCHAR:	PRINT_IN_PRTHS(char(NVal)); break;
			default:	PRINT_IN_PRTHS(NVal);
		}
	else if(SVal != NULL)	PRINT_IN_PRTHS(ValType == tENUM ? "NONE" : SVal);
	cout << EOL;
}

// Prints enum or combi values
//	return: number of printed symbols
BYTE Options::Option::PrintEnumVals()
{
	if( !ValRequired() )	return 0;

	char** vals = (char**)SVal;			// array of val images
	BYTE len  = BYTE(strlen(vals[0]));	// number of printed chars
	BYTE vCnt = BYTE(MaxNVal);			// number of values
	if(MinNVal)							// is range of values limited from below?
		vCnt -= BYTE(MinNVal) - 1;
	cout << " <" << vals[0];			// first val image from array of val images
	for(BYTE i=1; i<vCnt; i++) {
		cout << EnumDelims[ValType-tENUM] << vals[i];
		len += strlen(vals[i]);
	}
	cout << '>';
	return len + vCnt + 2;	// here vCnt denotes the number of printed delimiters
}

// Performs a case-insensitive search of given string value among enum values.
//	@val: input value as string
//	return: index of finded value in enum, or -1 if the value is not present in enum
int Options::Option::GetEnumInd (const char* val)
{
	for(int i=0; i<MaxNVal; i++)
		if( !_stricmp(val, ((const char**)SVal)[i]) )
			return i;
	return -1;
}

// Prints Usage params
void Options::Usage::Print(Option* opts) const
{
	if(Opt != vUNDEF)	// output option value
		opts[Opt].Print(false);
	else if(Par) {		// output parameter
		if(IsParOblig)	cout << BLANK << Par;
		else			PRINT_IN_PRTHS(Par);
		if(ParDescr)	// output parameter description
			cout << "\n  " << Par << " - " << ParDescr;
	}
	cout << endl;
}


// Check obligatory options and output message about first absent obligatory option.
//	return: -1 if some of obligatory options does not exists, otherwise 1
int Options::CheckObligs()
{
	for(int i=0; i<_OptCount; i++)
		if( _Options[i].CheckOblig() < 0 )	return -1;
	return 1;
}

// Set value of option 
//	@isWord: true if option is a word, false if option is a char
//	@opt: option with HYPHENs
//	@val: value of option
//	@isNextOpt: true if next parameter is option
//	@argIndex: the current index in argc; is increased by 1 in case of required value
//	return: 0 if success, 1 otherwise
int Options::SetOption (bool isWord, char* opt, char* val, bool isNextOpt, int *argIndex)
{
	int res;

	opt += int(isWord)+1;
	if(isWord) { if(strlen(opt) == 1)
		return PrintAmbigOpt(isWord, opt, Ambiguous, "excess '-'?"); }
	else if(strlen(opt) > 1)
		return PrintAmbigOpt(isWord, opt, Ambiguous, "forgot '-'?");

	for(int i=0; i<_OptCount; i++) {
		res = _Options[i].SetVal(isWord, opt, val, isNextOpt);
		if(!res) {
			if( _Options[i].ValRequired() )	(*argIndex)++;
			return 0;
		}
		if(res>0)		return 1;
	}
	cerr << "wrong option: " << OptToStr(isWord, opt) << EOL;
	return 1;
}

// Ouptuts option with error message to cerr
//	@isWord: true if option is long
//	@opt: option
//	@val: value or NULL
//	@msg: error message about value
//	@return: always 1
int Options::PrintWrongOpt(bool isWord, const char* opt, const char* val, const string msg)
{
	cerr << OptTitle << OptToStr(isWord, opt)
		 << SepSCl << (msg == strEmpty ? "wrong " + sValue : msg);
	if( val ) cerr << BLANK << val;
	cerr << EOL;
	return 1;
}

// Ouptuts ambiguous option with error message to cerr
//	@isWord: true if option is long
//	@opt: option
//	@headMsg: message at the beginning
//	@tailMsg: message at the end or NULL
int Options::PrintAmbigOpt(bool isWord, const char* opt, const char* headMsg, const char* tailMsg)
{
	cerr << headMsg << optTitle << OptToStr(isWord, opt);
	if( tailMsg )	cerr << SepSCl << tailMsg;
	cerr << EOL;
	return 1;
}

// Prints version
//	return: always 1
int	Options::PrintVersion()
{
	cout<< Product::Version
#ifndef _NO_ZLIB
		<< "\tzlib "ZLIB_VERSION
#endif
		<< endl;
	return 1;
}

// Prints 'usage' information
//	@title: if true prints title before information
//	return: 1 if title is settinf to true, 0 otherwise
int Options::PrintUsage (bool title)
{
	BYTE i, k;
	if( title )		cout << Product::Descr << endl << endl;
	
	// output 'Usage' section
	cout << "Usage:";
	for(k=0; k<_UsageCount; k++) {
		cout << TAB << Product::Title;	PRINT_IN_PRTHS("options");
		// output required options
		for(i=0; i<_OptCount; i++)
			if( _Options[i].Sign.Is(Signs::Oblig) )
				_Options[i].Print(false);
		// output parameters
		_Usages[k].Print(_Options);
	}
	cout << endl;

	// output options section
	cout << "Options:\n";
	for(k=0; k<_GroupCount; k++) {
		if(_OptGroups[k])	cout << _OptGroups[k] << ":\n";
		for(i=0; i<_OptCount; i++)
			if( _Options[i].OptGroup == k )
				_Options[i].Print(true);
	}
	return int(title);
}

// Returns command line.
//	@argc: count of main() parameters
//	@argv: array of main() parameters
string const Options::CommandLine(int argc, char* argv[])
{
	ostringstream oss;
	int i;
	for (i = 0; i < argc-1; i++)
		oss << argv[i] << BLANK;
	oss << argv[i];
	return oss.str();
}



// Parses and checks main() parameters and their values. Output message if some of them is wrong.
//	@argc: count of main() parameters
//	@argv: array of main() parameters
//	@obligPar: name of required application parameter or NULL if not required
//	return: index of first parameter (not option) in argv[], argc if it is absent, negative if  tokenize complets wrong
int Options::Tokenize(int argc, char* argv[], const char* obligPar)
{
	int i, res = 1;
	for (i = 1; i < argc; i++) {	// argv[0] is the path to the program
		if( argv[i][0] != HPH )	{
			if( i < argc-1 				// not a last option or parameter
			&& argv[i+1][0] == HPH ) {	// next word is an option
				cerr << argv[i] << ": neither option nor parameter"  << EOL;
				res = -1;
			}		
			break;
		}
		if( SetOption(
			argv[i][1] == HPH, argv[i], 
			i+1 < argc ? argv[i+1] : NULL,
			i+2 < argc ? argv[i+2][0]==HPH : false,
			&i) )
		{	res = -1; break; }
	}
	if( res > 0 )	res = CheckObligs();	// check required options
	if( res > 0 && obligPar && i == argc ) {
		cerr << Missing << obligPar << EOL;
		res = -1;
	}
	return i * res;
}

//void Options::GetOpt(int i)
//{
//	_Options[i].Print(true);
//}
/************************ end of class Options ************************/

/************************ class Err ************************/

//const char* Err::TREAT_BED_EXT = "after extension";	// clarifying message in the bed stretch operation

const char* Err::_msgs [] = {
/* NONE */		"WARNING",//WARNING,
/* MISSED */	"missing",
/* F_NONE */	"no such file",
/* FD_NONE */	"no such file or directory",
/* F_MEM */		"memory exceeded",
/* F_OPEN */	"could not open",
/* F_CLOSE */	"could not close",
/* F_READ */	"could not read",
/* F_BIGLINE */	"buffer is less than length of line",
/* FZ_MEM */	"not enough internal gzip buffer",
/* FZ_OPEN */	"wrong reading mode READ_ANY for zipped file",
/* FZ_BUILD */	"this build does not support zipped files",
/* F_WRITE */	"could not write",
#ifndef _FQSTATN
/* TF_FIELD */	"number of fields is less than expected",
///* TF_SPEC */	"wrong line format",
/* TF_EMPTY */	"no",
/* BP_BADEND */		"'start' position is equal or more than 'end'",
/* BP_NEGPOS */		"negative position",
#ifdef _BEDR_EXT
/* BR_RNAME */	"wrong read name format:",
#endif
///* FA_LONGLEN */"length of chromosome is more than ULONG_MAX",
#endif
#if defined _ISCHIP || defined _FQSTATN
/* FQ_HEADER */	"no '@' marker; missed header line",
/* FQ_HEADER2 */"no '+' marker; missed second header line",
#elif defined _DENPRO || defined _BIOCC
/* ARR_OUTRANGE */	"out of range",
/* SUM_EXCEED */	"exceeded digital limit while S calculated. Anormous density. Calculate P only",
#endif
/* EMPTY */		""
};

// Initializes _outText by cstring contained message kind of "<sender>: <text> <specifyText>".
void Err::set_message(const char* sender, const char *text, const char *specifyText)
{
	size_t senderLen = sender!=NULL ? strlen(sender) : 0;
	size_t textLen = strlen(text);
	size_t outLen = senderLen + textLen + 1 + strlen(SepSCl);
	if(specifyText)	outLen += strlen(specifyText) + 1;
	_outText = new char[outLen];
	memset(_outText, '\0', outLen);
	if(sender) {
		if(senderLen)
			strcpy(_outText, sender);
		strcat(_outText, SepCl);
	}
	strcat(_outText, text);
	if(specifyText) {
		strcat(_outText, sBLANK);
		strcat(_outText, specifyText);
	}
}

// Returns string containing file name and issue number.
//	@issName: name of issue
//	@issNumb: number of issue
//	@fName: file name
//const string Err::IssueNumbToStr(const string& issName, ULONG issNumb, const string& fName)
//{
//	string res = fName;
//	if(fName != strEmpty)	res += SepSCl;
//	return res + issName + BLANK + NSTR(issNumb);
//}

// Gets message "no @fileName.@fileExt[.gz] files in this directory"
const string Err::MsgNoFiles (const string & fileName, const string fileExt)
{
	return string("no " + fileName + fileExt + "[" + ZipFileExt + "] files in this directory");
}


Err::Err(const Err & src)
{
	_code = src._code;
	int size = strlen(src._outText) + 1;
	_outText = new char[size];
	strcpy(_outText, src._outText);
	//_specifyText = src._specifyText;
}

// Throws exception or outputs Err message.
//	@throwExc: if true then throws exception, otherwise outputs Err message
//	@eol: if true then carriage should be return while output Err message
void Err::Throw(bool throwExc, bool eol) {
	if(throwExc)	throw *this;
	else {		
		dout << _outText;
		if(eol)		dout << EOL;
		fflush(stdout);
	}
}

// Outputs warning with prefix "WARNING" and additional text, if it is setting.
void Err::Warning(string const& addText) {
	dout << _msgs[ErrWARNING];
	if(_outText[0] != ':') dout << SepCl;	// check if sender is not recorded
	dout << _outText;
	if( !addText.empty() )	
		dout << addText;
	dout << EOL;
	fflush(stdout);
}

/************************ end of class Err ************************/

/************************ class FileSystem ************************/

// Returns true if file system's object exists
//	@name: object's name
//	@st_mode: object's system mode
bool FS::IsExist(const char* name, int st_mode)
{
	struct_stat64 st;
	return ( !_stat64(name, &st) && st.st_mode & st_mode );
}

// Checks if file system's object doesn't exist
//	@name: object's name
//	@st_mode: object's system mode
//	@throwExcept: if true throws excwption,
//	otherwise outputs Err message as warning without EOL
//	@ecode: error's code
//	return: true if file or directory doesn't exist
bool FS::CheckExist	(const char* name,  int st_mode, bool throwExcept, Err::eCode ecode)
{
	if( IsExist(name, st_mode) )	
		return false;
	Err(ecode, name).Throw(throwExcept);
	return true;
}

// Searches through a file name for the any extention (slash|back-slash insensible).
//	@fname: file name
//	return: the index of the DOT mathed extention; otherwise npos
size_t FS::GetExtPos(const string &fname) {
	size_t pos = fname.find_last_of(DOT);
	if( pos != string::npos ) {
		if( pos==0 )					return string::npos;	// ./name
		if( pos==1 && fname[0]==DOT )	return string::npos;	// ../name
	}
	return pos;
}

// Gets size of file or -1 if file doesn't exist
LLONG FS::Size (const char* fname)
{
	struct_stat64 st;
	return _stat64(fname, &st) == -1 ? -1 : st.st_size;
}

// Gets real size of zipped file  or -1 if file cannot open; limited by UINT
LLONG FS::UncomressSize	(const char* fname)
{
	FILE *file = fopen(fname, "rb");	// "read+binary"
	if( file == NULL )		return -1;
	BYTE sz[4] = {0,0,0,0};
	_fseeki64(file, -4, SEEK_END);
	fread(sz, 1, 4, file);
	fclose(file);
	return (sz[3] << 3*8) + (sz[2] << 2*8) + (sz[1] << 8) + sz[0];
}

// Returns true if file has a specified  extension.
//	@fname: file name
//	@ext: extension includes dot symbol and can be composite
bool FS::HasExt	(const string& fname, const string& ext)
{ 
	size_t pos = fname.find(ext);
	return pos == string::npos ? false : fname.size() - pos == ext.size();
}

// Returns string containing real file extension (without zip extention).
//	@fname: pointer to the file name
//	return: string containing real file extension or empty string if no real extention
string const FS::GetExt(const char* fname) {
	const char * pdot = strrchr(fname, DOT);
	if( !pdot )		return strEmpty;
	if( strcmp(pdot, ZipFileExt.c_str()) )
		return string(pdot+1);				// no zip extention
	const char * pprevdot = pdot - 1;
	for(; pprevdot >= fname; pprevdot--)	// find previous DOT
		if( *pprevdot == DOT )	
			break;
	return pprevdot+1 == fname ? strEmpty : string(pprevdot+1, pdot-pprevdot-1);
}

// Returns file name without extentiom (slash|back-slash insensible)
string const FS::FileNameWithoutExt (const string& fname)
{
	size_t pos = GetExtPos(fname);
	return pos == string::npos ? fname : fname.substr(0, pos);
}

// Returns short file name by long one
//	@fname: long file name
string const FS::ShortFileName (const string& fname)
{
#ifdef OS_Windows
	if(fname.find(REAL_SLASH) != string::npos) {
		string tmp(fname);
		replace(tmp.begin(), tmp.end(), REAL_SLASH, SLASH);
		return tmp.substr(tmp.find_last_of(SLASH) + 1);
	}
#endif
	return fname.substr(fname.find_last_of(SLASH) + 1);
}

// Returns directory name by long file name
//	@fname: long file name
//	@addSlash: true if slash sould be added at the end
string const FS::DirName (const string& fname, bool addSlash)
{
#ifdef OS_Windows
	if(fname.find(REAL_SLASH) != string::npos) {
		string tmp(fname);
		replace(tmp.begin(), tmp.end(), REAL_SLASH, SLASH);
		return tmp.substr(0, tmp.find_last_of(SLASH) + int(addSlash));
	}
#endif
	return fname.substr(0, fname.find_last_of(SLASH) + int(addSlash));
}

// Returns the name of last subdirectory by long file name
//	@fname: long file name
string const FS::LastSubDirName (const string& fname)
{
	const string& dir = FS::DirName(fname, false);
	size_t pos = dir.find_last_of(SLASH);
	return pos == string::npos ? dir :
		dir.substr(dir.substr(0, pos).length() + 1);
}

// Returns the name ended by slash without checking the name
string const FS::MakePath(const string& name)
{
#ifdef OS_Windows
	if(name.find(REAL_SLASH) != string::npos) {
		string tmp(name + SLASH);
		replace(tmp.begin(), tmp.end(), REAL_SLASH, SLASH);
		return tmp;
	}
#endif
	return name + SLASH;
}

#ifndef _WIGREG
// Fills external vector of strings by file's names found in given directory
// Implementation depends of OS.
//	@files: external vector of strings that should be filled by file's names
//	@dirName: name of directory
//	@fileExt: file's extention as a choosing filter
//	@all: true if all files with given extention should be placed into external vector,
//	otherwise only one (any)
//	return: true if files with given extention are found
bool FS::GetFiles	(vector<string>& files, const string& dirName,
	const string& fileExt, bool all)
{
	BYTE count = 0;
#ifdef OS_Windows
	string fileTempl = FS::MakePath(dirName) + '*' + fileExt;
	WIN32_FIND_DATA ffd;

	HANDLE hFind = FindFirstFile( fileTempl.c_str(), &ffd );
	if( hFind == INVALID_HANDLE_VALUE )
		return false;
	if( all ) {
		// count files to reserve files capacity
		do	count++;
		while (FindNextFile(hFind, &ffd));
		files.reserve(count);
		// fill files
		hFind = FindFirstFile( fileTempl.c_str(), &ffd );
		do	files.push_back( string(ffd.cFileName) );
		while (FindNextFile(hFind, &ffd));
	}
	else
		files.push_back( string(ffd.cFileName) );
	FindClose(hFind);
	return true;
#else
	struct dirent *entry;
	DIR *dir = opendir(dirName.c_str());	// doesn't need to check because of programmes options inspection
	if( all ) {
		// count all files to reserve files capacity
		while( readdir(dir) )	count++;
		closedir(dir);
		if( !count )	return false;
		files.reserve(count);
		dir = opendir(dirName.c_str());
	}
	// fill files
	string name;
	while( entry = readdir(dir) ) {
		name = string(entry->d_name);
		if( HasExt(name, fileExt) ) {
			files.push_back(name);
			if( !all )	break;
		}
	}
	closedir (dir);
	return files.size() > 0;
#endif	// OS_Windows
}
#endif	// _WIGREG
/************************ end of class FileSystem ************************/

/************************ class FT ************************/

const FT::fType FT::Types[] = {
	{ "", strEmpty, strEmpty,			TabFilePar( 0, 0, '\0', NULL) },
	{ "bed", "read",	"reads",		TabFilePar( 6, 6, HASH, Chrom::Abbr) },
	{ "bed", "feature", "features",		TabFilePar( 3, 6, HASH, Chrom::Abbr) },
	{ "wig", "interval", "intervals",	TabFilePar( 2, 2, HASH, NULL) },
	{ "sam", strEmpty, strEmpty,		TabFilePar( 0, 0, HASH, NULL) },
	{ "fq", strEmpty, strEmpty,			TabFilePar( 0, 0, '\0', NULL) },
	{ "fa", strEmpty, strEmpty,			TabFilePar( 0, 0, '\0', NULL) }
};
const BYTE FT::Count = sizeof(FT::Types)/sizeof(FT::fType);

// Returns file format
//	@fName: file name (with case insensitive extension)
FT::eTypes FT::GetType(const char* fName)
{
	const string ext = FS::GetExt(fName);
	const char* c_ext = ext.c_str();
	for(int i = 2; i<Count; i++)	// start from ordinary bed
		if(!_stricmp(c_ext, Types[i].Extens))	return eTypes(i);
	return UNDEF;
}

// Validates file extension
//	@fName: file name (with case insensitive extension and [.gz])
//	@t: file type
//	@printfName: true if file name should be ptinted
//	@throwExc: true if throw exception, otherwise throw warning
//	return: true if file extension correspondes to file type
bool FT::CheckType(const char* fName, eTypes t, bool printfName, bool throwExc) { 
	if( GetType(fName) != (t == ABED ? BED : t) ) {
		Err("wrong extension", printfName ? fName : NULL).Throw(throwExc);
		return false;
	}
	return true;
}

#ifdef _ISCHIP
// Gets file extension, beginning at DOT and adding .gz if needed
//	@t: file type
//	@isZip: true if add ".gz"
const string FT::RealExt(eTypes t, bool isZip)
{
	string ext = string(".") + Types[t].Extens;
	if(isZip)	ext += ZipFileExt;
	return ext;
}
#endif

/************************ end of class FT ************************/

// Prints elapsed time
//	@elapsed: elapsed time in seconds
//	@title: string printed before time outpu
//	@watch: true if time should be printed as a stopwatch (with decimal places and without empty minutes)
//	@parentheses: if true then output time in parentheses
//	@isEOL: if true then ended output by EOL
void PrintTime(float elapsed, const char *title, bool watch, bool parentheses, bool isEOL)
{
	if(title)		dout << title;
	int hours = long(elapsed)/60;
	int mins = hours%60;
	hours /= 60;
	if(parentheses)	dout << '(';
	dout << setfill('0') << right;		// right couse it may be chanched by previuos output
	if(hours)	dout << setw(2) << hours << COLON;			// hours
	if(mins || !watch)	dout << setw(2) << mins << COLON;	// mins
	// secs
	if(watch)	dout << setw(5) << fixed << setprecision(2) << (elapsed - mins*60);
	else		dout << setw(2) << long(elapsed)%60;
	if(parentheses)	dout << ')';
	if(isEOL)	dout << EOL, fflush(stdout);
}

/************************  class TimerBasic ************************/
bool	TimerBasic::Enabled = false;

// Stops enabled timer and return elapsed wall time in seconds
long TimerBasic::GetElapsed() const
{
	time_t stopTime;
	time( &stopTime );
	return (long)difftime(stopTime, _startTime);
}

// Prints elapsed time interval
//	@elapsed: elapsed time in seconds 
//	@title: string printed before time output
//	@parentheses: if true then output time in parentheses
//	@isEOL: if true then ended output by EOL
void TimerBasic::Print(long elapsed, const char *title, bool parentheses, bool isEOL)
{
	PrintTime(elapsed, title, false, parentheses, isEOL);
}

/************************  end ofclass TimerBasic ************************/

/************************  class Timer ************************/
clock_t	Timer::_StartCPUClock;
/************************  end of class Timer ************************/
#ifdef _TEST
/************************  class Stopwatch ************************/

// Stops Stopwatch
//	@title: if not empty, and if instance was launched, output sum wall time with title
//	'const' to apply to constant objects
void Stopwatch::Stop(const string title) const
{
	if(!_isStarted)		return;
	_sumTime += GetElapsed();
	if(title!=strEmpty)	PrintTime(_sumTime, (title + sBLANK).c_str(), false, false, true);
}

/************************  end of class Stopwatch ************************/
#endif	// _TEST
/************************  class StopwatchCPU ************************/

// Stops StopwatchCPU
//	@title: string printed before time output
//	@print: if true time should be printed
//	@isEOL: if true then ended output by EOL
void StopwatchCPU::Stop(const char* title, bool print, bool isEOL)
{
	_sumclock += clock() - _clock;
	if(print)	PrintTime(float(_sumclock)/CLOCKS_PER_SEC, title, true, false, isEOL);
}

/************************  end of class StopwatchCPU ************************/


#ifdef _MULTITHREAD
/************************  class Mutex ************************/
pthread_mutex_t	Mutex::_mutexes[Mutex::Count];

void Mutex::Init() {
	for(BYTE i=0; i<Count; i++)
#ifdef OS_Windows
		InitializeCriticalSection(&_mutexes[i]);
#else
		pthread_mutex_init(&_mutexes[i], NULL);
#endif
}

void Mutex::Finalize() {
	for(BYTE i=0; i<Count; i++)
#ifdef OS_Windows
		DeleteCriticalSection(&_mutexes[i]);
#else
		pthread_mutex_destroy(&_mutexes[i]);
#endif
}

void Mutex::Lock(const eType type) {
#ifdef OS_Windows
	EnterCriticalSection(&_mutexes[type]);
#else
	pthread_mutex_lock(&_mutexes[type]);
#endif
}

void Mutex::Unlock(const eType type) {
#ifdef OS_Windows
	LeaveCriticalSection(&_mutexes[type]);
#else
	pthread_mutex_unlock(&_mutexes[type]);
#endif
}

/************************  end of class Mutex ************************/

/************************  class Thread ************************/

//Thread::Thread(thrRetValType(
//	#ifdef OS_Windows
//	__stdcall 
//	#endif
//	*proc)(void*), void *arglist)
//{
//#ifdef OS_Windows
//	//_thread = (HANDLE) _beginthread(proc, 0, arglist);	// gets unstable call WaitForSingleObject
//	_thread = (HANDLE) _beginthreadex(NULL, 0, proc, arglist, 0, NULL);
//#else
//	pthread_create(&_thread, NULL, proc, arglist);
//	//int code = pthread_create(&_thread, NULL, proc, arglist);
//	//if (code) { 
//	//	char buf[256]; 
//	//	strerror_r(code, buf, sizeof buf);
//	//	cout << "Thread constructor: " << buf << endl;
//	//}
//#endif
//}

/********************  end of class Thread *********************/
#endif	// _MULTITHREAD

/************************ class Chrom ************************/

const char*		Chrom::Abbr = "chr";
#ifndef _FQSTATN
const string	Chrom::UndefName = "UNDEF";
const string	Chrom::Short = "chrom";
const string	Chrom::Title = "chromosome";
const BYTE		Chrom::MaxAbbrNameLength = BYTE(strlen(Chrom::Abbr)) + MaxMarkLength;
const BYTE		Chrom::MaxShortNameLength = BYTE(Chrom::Short.length()) + MaxMarkLength;
const BYTE		Chrom::MaxNamedPosLength = 
	BYTE(strlen(Chrom::Abbr)) + MaxMarkLength + CHRLEN_CAPAC + 1;

chrid Chrom::_cID = UnID;	// user-defined chrom ID

// Sets chromosome's ID stated by user with validation.
//	@cID: chromosome's ID
//	return: true if the check has passed, otherwise do not set and print message to cerr
bool Chrom::SetStatedID(chrid cID)
{
	if( cID>64 && cID!=M && cID!=X && cID!=Y ) {
		cerr << Mark(cID) << ": wrong " << Title << "'s name\n";
		return false;
	}
	_cID = cID;
	return true;
}

// Returns a pointer to the first occurrence of C sunstring. Recurcive.
//	@str: C string to find in
//	@templ: C string to be located
//	@templLen: length of templ (extern because of avoiding recursive recalculate)
const char* SubStr(const char* str, const char* templ, int templLen)
{
	str = strchr(str, *templ);
	if( str )
		for(short i=1; i<templLen; i++)
			if( *++str != templ[i] )
				return SubStr(str, templ, templLen);
	return str;
}

chrid Chrom::ID(const char* cName, size_t prefixLen)
{
	if( !cName )		return UnID;
	cName += prefixLen;									// skip prefix
	//if(*cName == M || strchr(cName, USCORE))	return UnID;
	if( strchr(cName, USCORE))	return UnID;			// exclude chroms with '_'
	if(isdigit(*cName))			return atoi(cName);		// numeric chromosome
	//if(*cName <= '9')			return atoi(cName);		// numeric chromosome
	return islower(*cName) ? toupper(*cName) : *cName;	// letter's chromosome
}

// Locate chrom mark in string.
//	@str: string checked for chrom number
//	return: pointer to the chrom number in str, or a null pointer if Chrom::Abbr is not part of str.
const char* Chrom::FindMark(const char* str) {
	const char* substr = strstr(str, Abbr);
	return substr ? (substr+strlen(Abbr)) : NULL;
}

// Returns the length of prefix (substring before short chromosome's name)
//	@cLongName: long chromosome's name
//	return: length of substring before short chromosome's name, or -1 if short name is not finded
short Chrom::PrefixLength(const char* cLongName)
{
	// search from the beginning of the cName because
	// it simpler to find the first digit for multidigit number

	// start search from the occurrence of 'chr'
	const char* str = SubStr(cLongName, Abbr, strlen(Abbr));
	if( str )
		for(; *str; str++)
			if( isdigit(*str) || isupper(*str) )
				return str - cLongName;
	return -1;
}

//char* Chrom::LongToShortName(char* name) 
//{
//	short shift = PrefixLength(name);
//	return shift > 0 ? name + shift : NULL;
//}
#endif	// _FQSTATN

/************************ end of class Chrom ************************/

#if !defined _WIGREG && !defined _FQSTATN

/************************ class Read ************************/
readlen	Read::Len;				// length of Read
#if defined _ISCHIP || defined _BEDR_EXT

//readlen	Read::Len;				// length of Read
const char	Read::Strands[] = { '+', '-' };
const char	Read::NmPos1Delimiter = ':';
const char	Read::NmPos2Delimiter = '-';
//const char* Read::NmSuffMate1 = "/1";
//const char* Read::NmSuffMate2 = "/2";

#endif

#ifdef _ISCHIP

char	Read::SeqQuality;		// the quality values for the sequence (ASCII)
bool	Read::ExtReadName;		// true if Read name includes a position
short	Read::LimitN = vUNDEF;	// maximal permitted number of 'N' in Read or vUNDEF if all
const char Read::ToUp	= 'a' - 'A';
const char Read::Complements[] = {'T',0,'G',0,0,0,'C',0,0,0,0,0,0,'N',0,0,0,0,0,'A'};

void Read::Init(readlen rLen, bool extName, char seqQual, short limN)
{
	Len = rLen;
	ExtReadName = extName,
	SeqQuality = seqQual;
	if(limN < rLen)		LimitN = limN;
}

// Copies complemented Read.
void Read::CopyComplement(char* dst, const char* src)
{
	//for(dst += Read::Len-1; *dst != EOL; dst--, src++)
	//	*dst = Complements[*src - 'A' - (*src >= 'a' ? ToUp : 0)];
	for(char i=Read::Len-1; i>=0; i--, src++)
		dst[i] = Complements[*src - 'A' - (*src >= 'a' ? ToUp : 0)];
	//dst[Read::Len] = EOL;
}

// Checks Read for number of 'N'
//	@read: checked Read
//	return: -1 if Read is NULL, 0 if N limit is exceeded, 1 if success
//	return:	1: NULL Read
//			0: success
//			-2: N limit is exceeded
int Read::CheckNLimit(const char* read)
{
	if( !read )		return 1;
	if( LimitN != vUNDEF ) {
		readlen cntN = 0, i = Len;
		for(read += i-1; i; i--, read--)
			if( *read == cN	&& ++cntN > LimitN )
				return -2;
	}
	return 0;
}

// Prints Read values - parameters.
void Read::Print()
{
	cout << "Read" << SepDCl << "length = " << int(Len);
	if(ExtReadName)		cout << SepSCl << "name includes a position";
	cout << SepSCl << "'N'-limit" << SepCl;
	if( LimitN == vUNDEF )	cout << Options::BoolToStr(false);
	else					cout << LimitN;
	cout << EOL;
}

#endif

/************************ end of struct Read ************************/

#endif
