/*
bioCC is a fast advanced correlation calculator for basics bioinformatics file formats.
It computes signal- and Pearson correlation coefficients for densities, coverages and features. 
Program allows to know correlation coefficients for the whole genome, for each chromosome separately and for predefined regions inside chromosomes: again for the all regions and for each separately as well.
bioCC is designed to treat a bunch of files at once.
Copyright (C) 2016 Fedor Naumenko
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

// correlation coefficient notations
const char* CCs [] = { "P", "S" };	// add first empty str to print valid default
// output total coefficients notations
const char* totals [] = { "ADD", "ONLY" };
// sorting coefficients notations
const char* sorts [] = { "RGN", "CC" };

//{ char, str, optOblig, valRequired, type, group, defVal, minVal, maxVal, strVal, descr }
Options::Option Options::_Options [] = {
	{ 'a', "align",	0, false,tENUM, oINPUT,	FALSE,	0, 2, NULL,
	"input bed files are alignments" },
	{ 'd', "dupl",	0, true, tENUM,	oINPUT, TRUE,	0, 2, (char*)Options::Booleans,
	"accept duplicate reads. For alignments only" },
	{ 'g', "gen",	1, true, tNAME, oINPUT, vUNDEF, 0, 0, NULL,
	"genome size file, or reference genome (file name or directory). Required" },
	{ HPH, "gap-len",0,true, tINT,	oINPUT, 1000, 100, 100000, NULL,
	"minimal length of undefined nucleotides region in genome\nwhich is declared as a gap.\nIgnored for genome size file" },
	{ 'l', "list",	0, true, tNAME,	oINPUT, vUNDEF, 0, 0, NULL,
	"list of multiple input files.\nFirst (primary) file in list is comparing with others (secondaries)" },
	{ 'c',Chrom::Abbr,0,true,tCHAR, oTREAT, vUNDEF, 0, 0, NULL,
	"treat stated chromosome only" },
	{ 'r', "cc",	0, true, tCOMB,	oTREAT, ccP, ccP, ccS, (char*)CCs,
	"correlation coefficient, in any combination: ? - Pearson, ? - signal" },
	{ 's',"space",	0, true, tINT,	oTREAT, 50, 1, 1e4, NULL,
	"resolution: span in bps in which reads will be counted\nto define a density. For alignments only" },
	{ 'T', Total,	0, true, tENUM, oTREAT,totalOFF,totalADD,totalONLY, (char*)totals,
	"output total coefficients: ? - in addition, ? - solely" },
	{ 'f', "fbed",	0, true, tNAME, oTREAT_R, vUNDEF,	0, 0, NULL,
	"'template' bed file which features defines compared regions.\nIgnored for ordinary beds" },
	{ 'e',"exp-len",0, true, tINT,	oTREAT_R, 0, 0, 1e4, NULL,
	"length of expanding features in first ordinary bed file\nor in 'template' bed file" },
	{ HPH,"exp-step",0, true,tINT,	oTREAT_R, 0, 0, 500, NULL,
	"step of expanding features in first ordinary bed file;\nif 0 then no step calculation" },
	{ 'b', "bin-width",0,true,tFLOAT, oTREAT_R, 0, 0, 1.0F, NULL,
	"width of bins of histogram" },
	{ HPH, "sort",	0, true, tENUM,	oTREAT_R,rsOFF,rsR,rsC, (char*)sorts,
	"output coefficients for each region, sorted by:\n? - regions, ? - coefficients" },
	{ HPH, "norm",	0, true, tENUM,	oTREAT_R,TRUE, 0, 2, (char*)Options::Booleans,
	"normalize regions before calculation" },
	//{ HPH, "cross",	0, false,tBOOL, FALSE, 0, 0, NULL, "cross-correlation" },
	//{ HPH, "circ",	0, false,tBOOL, FALSE, 0, 0, NULL, "circular cross-correlation" },
	//{ 's', "step",	0, true, tINT, 2e5, 1, (float)INT_MAX, NULL, "minimal shift between outputted coefficient during cross-correlation." },
	{ HPH, "alarm",	0, false,tENUM, oOUTPUT,FALSE,	0, 2, NULL,
	"output features ambiguities, if they exist. Ignored for wigs" },
	{ HPH, "lac",	0, false,tENUM, oOUTPUT,FALSE,	0, 2, NULL, "laconic output" },
	{ HPH, "stat",	0, false,tENUM, oOUTPUT,FALSE,	0, 2, NULL,
	"output features ambiguities statistics, if they exist. Ignored for wigs" },
	{ 'o', "out",	0, false,tENUM, oOUTPUT,FALSE,	0, 2, NULL, HelpOutFile.c_str() },
	{ 't', "time",	0, false,tENUM, oOTHER,	FALSE,	0, 2, NULL, "output run time" },
	{ 'v', Version,	0, false,tVERS,	oOTHER, vUNDEF, 0, 0, NULL, "print program's version and quit" },
	{ 'h', "help",	0, false,tHELP,	oOTHER, vUNDEF, 0, 0, NULL, "print usage information and quit" }
};
const BYTE	Options::_OptCount = oHELP + 1;
const BYTE	Options::_UsageCount = 2;		// count of 'Usage' variants in help
const Options::Usage Options::_Usages[] = {	// content of 'Usage' variants in help
	{	vUNDEF,	" file0 fileN..."	},
	{	oFILE_LIST, StrEmpty	}
};

ofstream outfile;		// file ostream duplicated cout; inizialised by file in code
dostream dout(cout, outfile);	// stream's duplicator

/*****************************************/

int main(int argc, char* argv[])
{
	if (argc < 2)
	{ Options::PrintUsage(false);	return 0; }		// output tip

	short fileInd = Options::Tokenize(argc, argv);
	if( fileInd < 0 )	return 1;	// wrong option
	int ret = 0;					// main() return code

	FileList* flList = NULL;
	Timer timer(Options::GetBVal(oTIME));
	if( Options::GetBVal(oOUTFILE) )	outfile.open( OutFile.c_str() );
	timer.Start();
	Timer::StartCPU();
	try {
		char **inFiles;
		short cntInFiles;
		// check file-list of bed-files is set & exist
		const char *lsFileName = Options::GetSVal(oFILE_LIST);
		if( lsFileName ) {
			flList = new FileList( FS::CheckedFileName(lsFileName) );
			inFiles = flList->Files();
			cntInFiles = flList->Count();
		}
		else {
			inFiles = argv + fileInd;
			cntInFiles = argc - fileInd;
		}
		if( cntInFiles < 2 )			// check input files
			Err(Err::P_MISSED, StrEmpty, cntInFiles ? "secondaries "+InFiles : InFiles).Throw();
		chrid cID = Chrom::ID(Options::GetSVal(oCHROM));
		GenomeRegions gRgns(FS::CheckedFileDirName(Options::GetSVal(oGFILE)),
			&cID, Options::GetIVal(oGAPLEN));
		FileList secfList(inFiles+1, cntInFiles-1);	// list of second bed-files names

		CorrPair cPair(cID, inFiles[0], gRgns, Options::GetSVal(oFBED), secfList.Count()>1);
		for(short i=0; i<secfList.Count(); i++)
			cPair.CalcCC(secfList[i]);
	}
	catch(const Err &e)			{ ret = 1;	dout << e.what() << EOL; }
	catch(const exception &e)	{ ret = 1;	dout << e.what() << EOL; }
	catch(...)					{ ret = 1;	dout << "Unregistered error\n"; }
	if( flList )	delete flList;
	timer.Stop(false);
	Timer::StopCPU();
	if( outfile.is_open() )		outfile.close();	// in case of holding execution by user
//#ifdef OS_Windows
//	system("pause");
//#endif
	return ret;
}
