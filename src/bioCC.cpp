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
const string Product::Descr = "Correlation's calculator";

const string OutFile = string(Product::Title) +  "_out.txt";
const string HelpOutFile = "duplicate standard output to " + OutFile + " file";
const string InFiles = "input files";

enum eOptGroup	{ oINPUT, oTREAT, oTREAT_R, oOUTPUT, oOTHER };	// oOTHER should be the last 
const BYTE	Options::_GroupCount = oOTHER + 1;	// count of option groups in help

const char* Options::_OptGroups [] = {
	"Input", "Processing", "Regions processing", "Output", "Other"
};

// --cc option: correlation coefficient notations
const char* CCs [] = { "P", "S" };			// corresponds to CCkey::eCC
// --total option: total coefficients notations
const char* totals [] = { "ADD", "ONLY" };	// corresponds to eTotal; totalOFF is hidden
// --sort option: sorting type notations
const char* sorts [] = { "RGN", "CC" };		// corresponds to eRS; rsOFF is hidden
// --info option: types of info notations
const char* infos [] = { "NOTE", "STAT" };	// corresponds to eInfo; iOFF is hidden


//	{ char,	str,	Signs,	type,	group,	defVal,	minVal,	maxVal,	strVal,	descr }
// field 7: vUNDEF if value is prohibited
// field 6: vUNDEF if no default value should be printed
Options::Option Options::_Options [] = {
	{ 'a', "align",		0,	tENUM,	oINPUT,	FALSE,	vUNDEF, 2, NULL,
	"input bed files are alignments" },
	{ 'g', "gen",		1,	tNAME,	oINPUT, vUNDEF, 0, 0, NULL,
	"reference genome size file, genome library,\nor single nucleotide sequence. Required" },
	{ HPH, "gap-len",	0,	tINT,	oINPUT,	1000, 50, 100000, NULL,
	"minimal length of undefined nucleotide region in genome\nwhich is declared as a gap.\nIgnored for genome size file and for ordinary beds" },
	{ 'd', "dupl",		0,	tENUM,	oINPUT, TRUE,	0, 2, (char*)Options::Booleans,
	"accept duplicate reads. For alignments only" },
	{ HPH, "diff-sz",	0,	tENUM,	oINPUT, FALSE,	0, 2, (char*)Options::Booleans,
	"allow to ignore reads with different size. For alignments only" },
	{ 'l', "list",		0,	tNAME,	oINPUT, vUNDEF, 0, 0, NULL,
	"list of multiple input files.\nFirst (primary) file in list is comparing with others (secondaries)" },
	{ 'c', Chrom::Abbr,	0,	tCHAR,	oTREAT, vUNDEF, 0, 0, NULL,	"treat stated chromosome only" },
	{ 'r', "cc",		0,	tCOMB,	oTREAT,	CCkey::ccP, CCkey::ccP, CCkey::ccS, (char*)CCs,
	"correlation coefficient, in any combination: ? - Pearson, ? - signal" },
	{ 's', "space",		0,	tINT,	oTREAT,	50, 1, 1e4, NULL,
	"resolution: span in bp in which reads will be counted\nto define a density. For alignments only" },
	{ 'T',	Total,		0,	tENUM,	oTREAT,	totalOFF, totalADD, totalONLY, (char*)totals,
	"output total coefficients: ? - in addition, ? - solely" },
	{ 'f', "fbed",		0,	tNAME,	oTREAT_R,vUNDEF,	0, 0, NULL,
	"'template' ordinary bed file which features define compared regions.\nIgnored for ordinary beds" },
	{ 'e', "ext-len",	0,	tINT,	oTREAT_R,0, 0, 1e4, NULL,
	"length by which the features in primary file (for ordinary beds) or\nin 'template' (for alignments and wigs) extend in both directions" },
	{ HPH, "ext-step",	0,	tINT,	oTREAT_R,0, 0, 500, NULL,
	"step of extending features in primary bed file;\nif 0 then no step calculation. For ordinary beds only" },
	{ 'b', "bin-width",	0,	tFLOAT,	oTREAT_R,0, 0, 1.0F, NULL,
	"width of bins of histogram" },
	{ HPH, "sort",	0,	tENUM,	oTREAT_R,rsOFF,	rsR, rsC, (char*)sorts,
	"output region coefficients, sorted by:\n? - regions, ? - coefficients" },
	{ HPH, "norm",	0,	tENUM,	oTREAT_R,TRUE,	0, 2, (char*)Options::Booleans,
	"normalize regions before calculation. Ignored for ordinary beds" },
	//{ HPH, "cross",	0, tBOOL, FALSE, 0, 0, NULL, "cross-correlation" },
	//{ HPH, "circ",	0, tBOOL, FALSE, 0, 0, NULL, "circular cross-correlation" },
	//{ 's', "step",	0, tINT, 2e5, 1, (float)INT_MAX, NULL, "minimal shift between outputted coefficient during cross-correlation." },
	{ 'i', "info",	0,	tENUM, oOUTPUT,	Bed::iOFF,	Bed::iNOTE, Bed::iSTAT, (char*)infos,
	"output summary information about feature ambiguities, if they exist:\n? - notice, ? - statistics. Ignored for wigs" },
	{ 'w', "warn",	0,	tENUM, oOUTPUT,	FALSE,	vUNDEF, 2, NULL,
	"output each feature ambiguity, if they exist. Ignored for wigs" },
	{ HPH, "lac",	0,	tENUM,	oOUTPUT,FALSE,	vUNDEF, 2, NULL, "laconic output" },
	{ 'o', "out",	0,	tENUM,	oOUTPUT,FALSE,	vUNDEF, 2, NULL, HelpOutFile.c_str() },
	{ 't', "time",	0,	tENUM,	oOTHER,	FALSE,	vUNDEF, 2, NULL, "output run time" },
	{ 'v', Version,	0,	tVERS,	oOTHER,	vUNDEF, vUNDEF, 0, NULL, "print program's version and quit" },
	{ 'h', "help",	0,	tHELP,	oOTHER,	vUNDEF, vUNDEF, 0, NULL, "print usage information and quit" }
};
const BYTE	Options::_OptCount = oHELP + 1;
const BYTE	Options::_UsageCount = 2;		// count of 'Usage' variants in help
const Options::Usage Options::_Usages[] = {	// content of 'Usage' variants in help
	{	vUNDEF,	" file1 fileN..."	},
	{	oFILE_LIST, StrEmpty	}
};

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
	const ChromSizes* cSizes = NULL;
	const BedF* templ = NULL;
	const FileList* fList = NULL;
	if( Options::GetBVal(oOUTFILE) )	outfile.open( OutFile.c_str() );
	Timer::Enabled = Options::GetBVal(oTIME);
	Timer timer;
	try {
		char **inFiles;
		short cntInFiles;
		const char* gName = FS::CheckedFileDirName(oGFILE);
		GenomeRegions gRgns(gName, cSizes, Options::GetIVal(oGAPLEN));

		// check file-list is set & exist
		const char* fListName = Options::GetSVal(oFILE_LIST);
		if( fListName ) {
			fList = new FileList( FS::CheckedFileName(fListName) );
			inFiles = fList->Files();
			cntInFiles = fList->Count();
			if( cntInFiles == 0 )
				Err(Err::P_MISSED, StrEmpty, InFiles + " (no uncommented line in " + 
					string(fListName) + ")").Throw();
		}
		else {
			inFiles = argv + fileInd;
			cntInFiles = argc - fileInd;
		}
		if( cntInFiles < 2 )			// check input files
			Err(Err::P_MISSED, StrEmpty, cntInFiles ? "secondaries "+InFiles : InFiles).Throw();

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
