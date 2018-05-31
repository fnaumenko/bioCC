/*
	bioCC is a fast advanced correlation calculator for basics bioinformatics file formats.
	It computes Signal and Pearson correlation coefficients for densities, coverages and features.
	Program allows to know correlation coefficients for the whole genome, for each chromosome
	separately and for predefined regions inside chromosomes:
	again for the all regions and for each region separately as well.
	
	bioCC is designed to treat a bunch of files at once.

	Copyright (C) 2017 Fedor Naumenko (fedor.naumenko@gmail.com)

	This program is free software. It is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY;
	without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
	See the	GNU General Public License for more details.
 */

#include "Calc.h"
#include <fstream>

using namespace std;

const string Product::Title = "bioCC";
const string Product::Version = "1.0";
const string Product::Descr = "Correlation's calculator for bed/wig files";

const string OutFile = string(Product::Title) +  "_out.txt";
const string HelpOutFile = "duplicate standard output to " + OutFile + " file";
const string InFiles = "input files";

enum eOptGroup	{ oINPUT, oTREAT, oTREAT_R, oOUTPUT, oOTHER };	// oOTHER should be the last 

const char* Options::_OptGroups [] = {
	"Input", "Processing", "Region processing", "Output", "Other"
};
const BYTE Options::_GroupCount = sizeof(Options::_OptGroups)/sizeof(char*);

// --cc option: correlation coefficient notations
const char* CCs [] = { "P", "S" };			// corresponds to CCkey::eCC
// --total option: total coefficients notations
const char* prCCs [] = { "IND", "TOT" };	// corresponds to eTotal; totalOFF is hidden
// --sort option: sorting type notations
const char* sorts [] = { "RGN", "CC" };		// corresponds to eRS; rsOFF is hidden
// --info option: types of info notations
const char* infos [] = { "LAC", "NM", "CNT", "STAT" };	// corresponds to eInfo; iNONE is hidden

const char* ForAligns = "For the alignments only";
const char* IgnoreBed = "Ignored for the ordinary beds";

//	{ char,	str,	Signs,	type,	group,	defVal,	minVal,	maxVal,	strVal,	descr, addDescr }
// field 7: vUNDEF if value is prohibited
// field 6: vUNDEF if no default value should be printed
Options::Option Options::_Options [] = {
	{ 'a', "align",		0,	tENUM,	oINPUT,	FALSE,	vUNDEF, 2, NULL,
	"input bed files are alignments", NULL },
	{ 'g', "gen",		1,	tNAME,	oINPUT, vUNDEF, 0, 0, NULL,
	"chromosome sizes file, reference genome library,\nor single nucleotide sequence.", NULL },
	{ HPH, "gap-len",	0,	tINT,	oINPUT,	1000, 50, 100000, NULL,
	"minimal length of undefined nucleotide region in genome\nwhich is declared as a gap.\nIgnored for the chromosome sizes file and for the ordinary beds", NULL },
	{ 'd', "dupl",		0,	tENUM,	oINPUT, TRUE,	0, 2, (char*)Options::Booleans,
	"accept duplicate reads.", ForAligns },
	//{ HPH, "diff-sz",	0,	tENUM,	oINPUT, FALSE,	0, 2, (char*)Options::Booleans,
	//"allow to ignore reads with different size.", ForAligns },
	{ 'l', "list",		0,	tNAME,	oINPUT, vUNDEF, 0, 0, NULL,
	"list of multiple input files.\nFirst (primary) file in list is comparing with others (secondary)", NULL },
	{ 'c', Chrom::Abbr,	0,	tNAME,	oTREAT, vUNDEF, 0, 0, NULL,	"treat specified chromosome only", NULL },
	{ 'r', "cc",		0,	tCOMB,	oTREAT,	CCkey::ccP, CCkey::ccP, CCkey::ccS, (char*)CCs,
	"correlation coefficient, in any combination: ? - Pearson, ? - signal", NULL },
	{ 's', "space",		0,	tINT,	oTREAT,	100, 2, 1e4, NULL,
	"resolution: span in bp by which reads will be counted\nto define a density.", ForAligns },
	{ 'p',	"pr-cc",	0,	tCOMB,	oOUTPUT, Results::cIND, Results::cIND, Results::cTTL, (char*)prCCs,
	"print coefficient, in any combination:\n? - for each chromosome individually, ? - total", NULL },
	{ 'f', "fbed",		0,	tNAME,	oTREAT_R, vUNDEF,	0, 0, NULL,
	"'template' ordinary bed file which features define compared regions.\n", IgnoreBed},
	{ 'e', "ext-len",	0,	tINT,	oTREAT_R,0, 0, 1e4, NULL,
	"length by which the features in primary file (for ordinary beds) or in\n'template' (for alignments and wigs) will be extended in both directions\nbefore treatment", NULL },
	{ HPH, "ext-step",	0,	tINT,	oTREAT_R,0, 0, 500, NULL,
	"step of extending features in primary bed file;\nif 0 then no step calculation. For the ordinary beds only", NULL },
	{ 'b', "bin-width",	0,	tFLOAT,	oTREAT_R,0, 0, 1.0F, NULL,
	"width of the histogram bin", NULL },
	{ HPH, "sort",	0,	tENUM,	oOUTPUT, rsOFF,	rsR, rsC, (char*)sorts,
	"print region coefficients, sorted by:\n? - regions, ? - coefficients", NULL },
	{ HPH, "norm",	0,	tENUM,	oTREAT_R,TRUE,	0, 2, (char*)Options::Booleans,
	"normalize regions before calculation.", IgnoreBed },
	//{ HPH, "cross",	0, tBOOL, FALSE, 0, 0, NULL, "cross-correlation" },
	//{ HPH, "circ",	0, tBOOL, FALSE, 0, 0, NULL, "circular cross-correlation" },
	//{ 's', "step",	0, tINT, 2e5, 1, (float)INT_MAX, NULL, "minimal shift between outputted coefficient during cross-correlation." },
	{ 'i', "info",	0,	tENUM, oOUTPUT,	Obj::iNM, Obj::iLAC, Obj::iSTAT, (char*)infos,
	"print information about file:\n? - laconic, ? - name only, ? - number of items, ? - statistics", NULL },
	{ 'w', "warn",	0,	tENUM, oOUTPUT,	FALSE,	vUNDEF, 2, NULL,
	"print each file's item ambiguity, if they exist.", NULL },
	{ 'o', "out",	0,	tENUM,	oOUTPUT,FALSE,	vUNDEF, 2, NULL, HelpOutFile.c_str(), NULL },
	{ 't', "time",	0,	tENUM,	oOTHER,	FALSE,	vUNDEF, 2, NULL, "print run time", NULL },
	{ 'v', Version,	0,	tVERS,	oOTHER,	vUNDEF, vUNDEF, 0, NULL, "print program's version", NULL },
	{ 'h', "help",	0,	tHELP,	oOTHER,	vUNDEF, vUNDEF, 0, NULL, "print usage information", NULL }
};
const BYTE	Options::_OptCount = sizeof(Options::_Options)/sizeof(Options::Option);

const Options::Usage Options::_Usages[] = {	// content of 'Usage' variants in help
	{ vUNDEF, "file1 file2 ...", true, NULL },
	{ oFILE_LIST, NULL, true, NULL }
};
const BYTE Options::_UsageCount = sizeof(Options::_Usages)/sizeof(Options::Usage);

ofstream outfile;		// file ostream duplicated cout; inizialised by file in code
dostream dout(cout, outfile);	// stream's duplicator

/*****************************************/
int main(int argc, char* argv[])
{
	if (argc < 2)	return Options::PrintUsage(false);			// output tip
	int fileInd = Options::Tokenize(argc, argv);
	if( fileInd < 0 )	return 1;								// wrong option
	if(!Chrom::SetStatedID(Options::GetSVal(oCHROM))) return 1;	// wrong chrom name

	int ret = 0;					// main() return code
	ChromSizes* cSizes = NULL;
	const BedF* templ = NULL;
	const FileList* fList = NULL;
	if( Options::GetBVal(oOUTFILE) )	outfile.open( OutFile.c_str() );
	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		char **inFiles;
		short cntInFiles;
		const char* gName = FS::CheckedFileDirName(oGFILE);
		DefRegions gRgns(gName, &cSizes, Options::GetIVal(oGAPLEN));

		// check file-list is set & exist
		const char* fListName = Options::GetSVal(oFILE_LIST);
		if( fListName ) {
			fList = new FileList( FS::CheckedFileName(fListName) );
			inFiles = fList->Files();
			cntInFiles = fList->Count();
			if( cntInFiles == 0 )
				Err(Err::MISSED, NULL, InFiles + " (no uncommented line in " + 
					string(fListName) + ")").Throw();
		}
		else {
			inFiles = argv + fileInd;
			cntInFiles = argc - fileInd;
		}
		if( cntInFiles < 2 )			// check input files
			Err(Err::MISSED, NULL, cntInFiles ? "secondary " + InFiles : InFiles).Throw();

		CorrPair cPair(inFiles[0], cSizes, gRgns, Options::GetSVal(oFBED), cntInFiles > 2);
		for(short i=1; i<cntInFiles; i++)
			cPair.CalcCC(inFiles[i]);
	}
	catch(const Err &e)			{ ret = 1;	dout << e.what() << EOL; }
	catch(const exception &e)	{ ret = 1;	dout << e.what() << EOL; }
	catch(...)					{ ret = 1;	dout << "Unregistered error\n"; }
	if(cSizes)	delete cSizes;
	if(templ)	delete templ;
	if(fList)	delete fList;
	timer.Stop(false);
	if( outfile.is_open() )		outfile.close();	// in case of holding execution by user
//#ifdef OS_Windows
//	system("pause");
//#endif
	return ret;
}


/************************ class FileList ************************/
#ifdef OS_Windows
// Returns true if 'name' is file's pattern
bool	IsFilePattern	(const char* name)
{
	return strchr(name, '*') != NULL || strchr(name, '?') != NULL;
}

string GetPath	(const LPCTSTR name)
{
	const char* pch = strrchr(name, '/');
	return pch ? string(name, pch-name+1) : strEmpty;
}

// Fills vector by filenames according name template.
// Works only if _UNICODE is not defined.
//	@files: vector of file names
//	@templ: name template, which can include '*' and '?' marks
void FillFilesByTemplate(vector<string>& files, const LPCTSTR templ)
{
	string path = GetPath(templ);
	WIN32_FIND_DATA ffd;
 
	// check directory and estimate listFileNames capacity
	HANDLE hFind = FindFirstFile( templ, &ffd );
	if( hFind == INVALID_HANDLE_VALUE )		
		Err("bad file or content", templ).Throw();
	if( files.capacity() == 0 ) {		// count files to reserve capacity
		short count = 1;
		for(; FindNextFile(hFind, &ffd); count++);
		files.reserve(count);
		hFind = FindFirstFile( templ, &ffd );
	}
	// fill the list
	do	files.push_back(path + string(ffd.cFileName));	//  works only if _UNICODE isn't defined
	while (FindNextFile(hFind, &ffd));
	FindClose(hFind);
}
#endif	// OS_Windows

FileList::FileList(char* files[], short cntFiles) : _files(NULL), _memRelease(true)
{
#ifdef OS_Windows
	short i;
	bool hasTemplate;
	for(i=0; i<cntFiles; i++)
		if( hasTemplate = IsFilePattern(files[i]) )
			break;
	if( hasTemplate ) {
		// First we fill vector of file names, second initialize _files by this vector.
		// It needs because files[] may contain file names and file template (pattern) as well,
		// so in case of Windows we don't know the total amount of files beforehand.
		vector<string> tmpFiles;
		if( cntFiles > 1 )
			tmpFiles.reserve(cntFiles);
		// else if it's a single template name, capacity will be reserved in FillFilesByTemplate(),
		// or if it's a single common name, capacity will not be reserved at all

		for(i=0; i<cntFiles; i++)
			if( IsFilePattern(files[i]) )
				FillFilesByTemplate(tmpFiles, files[i]);
			else		// real list of files
				tmpFiles.push_back(files[i]);

		_files = new char*[_count=tmpFiles.size()];
		for(i=0; i<_count; i++) {
			_files[i] = (char*)malloc(tmpFiles[i].length()+1);
			strcpy(_files[i], tmpFiles[i].c_str());
		}
	}
	else 
#endif	// OS_Windows
	{
		_files = files;
		_count = cntFiles;
		_memRelease = false;
	}
}

FileList::FileList(const char* fName) : _files(NULL), _memRelease(true)
{
	TabFile file(fName);
	ULONG cntLines;
	char *dstStr;
	const char *srcStr;
	// no needs to check since aborting invalid file is set
	const char *currLine = file.GetFirstLine(&cntLines);
	vector<char*> tmpFiles;	// temporary vector because cntLines is not proof, but estimated capacity
	
	_count = 0;
	tmpFiles.reserve(cntLines);
	while(currLine!=NULL) {
		srcStr = file.StrField(0);
		dstStr = (char*)malloc(strlen(srcStr)+1);
		strcpy(dstStr, srcStr);
		tmpFiles.push_back(dstStr);
		_count++;
		currLine=file.GetLine();
	}
	_files = new char*[_count];
	for(short i=0; i<_count; i++)
		_files[i] = tmpFiles[i];
}

/************************ end of class FileList ************************/

