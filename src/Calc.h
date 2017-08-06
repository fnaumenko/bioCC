#pragma once

#include "Data.h"
#include "bioCC.h"

enum eTotal {	// defines types of outputting total CC
	totalOFF = vUNDEF,	// not set: is not pointed in command line
	totalADD = 1,		// add total CC
	totalONLY= 2,		// only total CC
};

enum eRS {	// defines types of outputting template regions (features) CC
	rsOFF = vUNDEF,		// not set: is not pointed in command line
	rsR = 1,			// sorted by regions
	rsC	= 2				// sorted by coefficients
};

typedef pair<Regions::Iter, Regions::Iter> RegionsRange;

class ShellGenomeRegions : public Chroms<RegionsRange>
/*
 * Wrapper class for  Chroms<RegionsRange> object to implement the interface for the common * functionality of Bedf and GenomeRegions.
 * Allows pass through these objects via iterator in a common way.
 */
{
public:
	// Returns an iterator pointing to the first Region of pointed chrom
	//	@it: chromosome's iterator
	const Regions::Iter ChromBegin(Chroms<RegionsRange>::cIter it) const {
		return it->second.first; }

	// Returns an iterator referring to the past-the-end Region of pointed chrom
	//	@it: chromosome's iterator
	const Regions::Iter ChromEnd(Chroms<RegionsRange>::cIter it)	const { 
		return it->second.second; }

	// Constructor by GenomeRegions
	ShellGenomeRegions(const GenomeRegions& gRgns) {
		for(GenomeRegions::cIter it=gRgns.cBegin(); it!=gRgns.cEnd(); it++)
			AddVal(CID(it), RegionsRange(it->second.Begin(), it->second.End()));
	}

	// Constructor by BedF
	ShellGenomeRegions(const BedF& bedF) {
		for(BedF::cIter it=bedF.cBegin(); it!=bedF.cEnd(); it++)
			AddVal(CID(it), RegionsRange(bedF.FeaturesBegin(it), bedF.FeaturesEnd(it)));
	}

	// Gets number of regions of pointed chrom
	//	@it: chromosome's iterator
	inline chrlen RegionsCount(Chroms<RegionsRange>::cIter it) const { 
		return it->second.second - it->second.first; }

	// Gets summary length of regions of pointed chrom
	//	@it: chromosome's iterator
	chrlen RegionsLength(Chroms<RegionsRange>::cIter it) const { 
		// need singleton?
		chrlen len = 0;
		for(Regions::Iter iter=it->second.first; iter<it->second.second; iter++)
			len += iter->Length();
		return len;
	}
};

class Results : public Chroms<CC>
/*
 * 'Results' accumulates calculated coeficients for each chromosome and total,
 * and provides theirs output.
 */
{
private:
	CC _total;
	bool _giveLocal,	// true if local CC should be calculated
		 _giveTotal;	// true if total CC should be calculated

public:
	// creates inctance with defined local/total inquires
	inline Results(eTotal etotal, bool notSingleChrom) :
		_giveLocal(etotal != totalONLY),
		_giveTotal(etotal != totalOFF && notSingleChrom) {}

	// Gets inquire about local CC: true if should be calculated
	inline bool GiveLocal() const { return _giveLocal; }
	// Gets inquire about total CC: true if should be calculated
	inline bool GiveTotal() const { return _giveTotal; }

	// Sets total results
	inline void	SetTotal(CC total) { _total = total; }

	// Prints results for each chromosome and total or total only (if it was defined).
	//	@printTitles: if true, print titles, otherwise does not only in case of single result
	void Print(bool printTitles);
};

struct ChromMap
/*
 * 'ChromMap' is a chromosome wig map with Treated sign.
 * chrom map is a simle array, contents coverage/debsity value for each minimal span
 */
{
	bool	Treated;	// true if chromosome is treating
	arrchrlen Map;		// chrom map

	inline ChromMap() : Treated(true) {}
	
	inline void Reserve(chrlen size) { Map.Reserve(size); }
};

class ChromsMap : public Obj, public Chroms<ChromMap>
/*
 * 'ChromsMap' is a container of ChromMap's arrays (one array for one chromosome).
 * Provides methods to calculate R.
 */
{
protected:
	class Pocket
	/* 
	 * keeps data pointers and temporary variables needed for reading constructor only
	 */
	{
	protected:
		GenomeRegions& __gRgns;	// initial genome regions
		arrchrlen*  _map;		// current chroms map (to avoid call indexator each time)
		chrlen		_cSize;		// size of current adding chromosome

	public:
		inline Pocket(GenomeRegions& gRgns) : __gRgns(gRgns) {}

		// sets and reserves current chrom's map and size
		inline void Reserve(chrid cID, arrchrlen* map) {
			(_map = map)->Reserve((_cSize = __gRgns[cID].LastEnd())/_Space);
		}
	};

	static BYTE	_Space;

	double	_binWidth;	// width of bins of histogram; if 0, no histogram
	eRS		_printFRes;	// sign to print results for each feature from 'template' and how to sort it

private:
	// Gets a pair of total genome means between this and second ChromsMap
	//	@gRgn: genome regions
	//	@map: second ChromsMap
	//	return: a pair of means for each set of chroms
	pairDbl GetGenomeMean(const GenomeRegions& gRgn, const ChromsMap& map) const;

public:
	inline ChromsMap() : 
		_binWidth (Options::GetFVal(oBINWIDTH)),
		_printFRes(eRS(Options::GetIVal(oFRES))) {}

	//inline arrchrlen & operator[] (chrid cID) { return At(cID).Map; }
	inline const arrchrlen & operator[] (chrid cID) const { return At(cID).Map; }

	// Calculates r for each region and fills results
	//	@cc: type of correlation coefficient
	//	@wig: ChromsMap object to correlate with
	//	@shGRgns: shell of treated chrom's regions
	//	@results: object to fill results
	void CalcRegionsR(eCC ecc, const ChromsMap& wig,
		const ShellGenomeRegions& shGRgns, Results& results);

	// Calculates r and fills results
	//	@cc: type of correlation coefficient
	//	@wig: ChromsMap object to correlate with
	//	@gRgns: real chrom's regions
	//	@bedF: template with defined regions or NULL
	//	@results: object to fill results
	void CalcR(eCC ecc, const ChromsMap& wig, const GenomeRegions& gRgns, const BedF* bedF, Results & results);

#ifdef DEBUG
	// Prints regions.
	//	@rgnCnt: number of regions to print or all if not specified
	void Print(chrlen rgnCnt=0) const;
#endif
};

class WigMap : public ChromsMap
/*
 * 'WigMap' encapsulates wig variableStep format.
 * It's a wig-specialized shell of ChromsMap.
 */
{
private:
	class WigPocket : public ChromsMap::Pocket
	{
	public:
		WigPocket(GenomeRegions& gRgns) : Pocket(gRgns) {}

		// Adds region for current chromosome.
		//	@start:	region's start position
		//	@size:	region's length
		//	@val:	region's value
		//	return: 1 if region is added; otherwise 0
		BYTE AddRegion(chrlen start, chrlen size, chrlen val);
	};

	// Adds chromosome and set it as current
	//	@cID: adding chromosome's ID
	//	@pocket: initializing temporary variables
	inline void AddChrom (chrid cID, Pocket& pocket) {
		pocket.Reserve(cID, &(AddEmptyClass(cID).Map));
	}

public:
	// Creates wigMap object.
	//	@fName: file name
	//	@incID: readed chromosome's ID or Chrom::UnID if all
	//	@primary: if true object is primary
	//	@printFileName: if true print file name in exception
	//	@gRgns: genome's regions
	WigMap(const char* fileName, chrlen incID, bool primary, bool printFileName,
		GenomeRegions& gRgns);
};

class DensMap : public ChromsMap
/*
 * 'DensMap' encapsulates density map.
 * It's a density-specialized shell of ChromsMap.
 */
{
private:
	class DensPocket : public ChromsMap::Pocket
	{
	private:
		const BedR& _bedR;
		const readlen _halfRLen;
		chrlen	_wCurrLen,		// length of uncsanning rest of current window
				_wi;			// current window's index
		//ULONG	_rtotalCnt;		// total number of reads in distribution
		BedR::cItemsIter	_currIt,
							_endIt;

		// Fills number or reads from window started with @start position
		void ScanWindow(chrlen start);

	public:
		DensPocket(const BedR& bedR, GenomeRegions& gRgns) : 
			_bedR(bedR), _halfRLen(bedR.ReadLen()/2), 
				//_rtotalCnt(0),
				Pocket(gRgns) {}

		// Fills number or reads from region @rgn
		//void ScanRegion(const Region&);
		void ScanChrom();

		// Adds chrom by BedR iter
		//	@it: BedR iter
		//	@map: pointer to the chroms map
		void AddChrom (BedR::cIter it, arrchrlen* map);

		//inline ULONG TotalReadsCnt() const { return _rtotalCnt; }
	};

	// Adds chromosome and set it as current
	//	@it: chromosome's bed iterator
	//	@pocket: initializing temporary variables
	inline void AddChrom (BedR::cIter it, DensPocket& pocket) {
		pocket.AddChrom(it, &(AddEmptyClass(CID(it)).Map));
	}

public:
	//DensMap(const BedR& bedR, bool primary, bool printFileName, GenomeRegions& gRgns);
	DensMap(const BedR& bedR, GenomeRegions& gRgns);
};

struct ChromRanges : public ChromItemsInd	
// Structure represented a set of chromosome's ranges.
// To keep chromosome's number and first/last ranges indexes BedF::ChromItemsInd is used;
// FirstInd/LastInd are keeping first/last ranges indexes,
// and lengths of all features from bed1, bed2 are keeping additionally for efficiency.
{
	chrlen FeatrsLen1;		// length of all features of chromosome in bed1
	chrlen FeatrsLen2;		// length of all features of chromosome in bed2

	inline ChromRanges() : ChromItemsInd(), FeatrsLen1(0), FeatrsLen2(0) {}

	inline ChromRanges(chrlen firstInd, chrlen lastInd, chrlen len1, chrlen len2) : 
		ChromItemsInd(firstInd, lastInd), FeatrsLen1(len1), FeatrsLen2(len2) {}
};

class JointedBeds : Chroms<ChromRanges>
/*
 * 'JointedBeds' represents two bed-files as a chromosomes collection
 * and theirs joint features (ranges).
 * This is fast but complicated implementation of calculating algorithm.
 */
{
private:
	class PairR
	/*
	 * 'PairR' pepresetns a pair of R classes, 
	 * which accumulates variances & covariance and calculated r.
	 * First in pare is signal R, second - Pearson R.
	 */
	{
	private:
		// Class 'R' accumulates variances & covariance and calculated r
		class R {
		private: 
			double	_var1, _var2,	// variances 
					_cov;			// covariance: SUM( X-Xmean)^2 * (Y-Ymean)^2 )
	// Mean value for every range may accept only one of two values: 0-mean or 1-mean,
	// and they are static within chromosome.
	// So for efficiency we can keep theirs in 3 arrays:
	// two square's arrays: arrays of squared values from each bed;
	// there are only 2 combinations, but they are duplicated for the simplicity acces,
	// and one crossing array: array of all combinations
	// of multiplications different values from bed1 and bed2
			double	_sqMeans1[4], _sqMeans2[4],	// square's arrays for bed, bed2
					_crossMeans[4];				// crossing array
		public:
			inline R () : _var1(0), _var2(0), _cov(0) {}
			// Keeps mean values for bed1, bed2. 
			//	@clear: if true, clear instance for treatment of new chromosome
			void Init(double mean1, double mean2, bool clear);
			// Accumulates next length of range
			void Increment(chrlen len, char val);
			// Returns coefficient of correlation
			inline double Get() const {	return CCkey::CalcR(_var1, _var2, _cov); }
		};

		R	_s;	// signal R
		R	_p;	// Pearson R
		eCC	_ecc;

	public:
		inline PairR (eCC ecc) : _ecc(ecc) {}
		// Keeps mean values for bed1, bed2. 
		//	@clear: if true, clear instance for treatment of new chromosome
		void Init(double mean1, double mean2, bool clear) {
			if( CCkey::IsS(_ecc) )	_s.Init(0, 0, clear);
			if( CCkey::IsP(_ecc) )	_p.Init(mean1, mean2, clear);
		}
		// Accumulates next length of range
		void Increment(chrlen len, char val) {
			if( CCkey::IsS(_ecc) )	_s.Increment(len, val);
			if( CCkey::IsP(_ecc) )	_p.Increment(len, val);
		}
		// Returns coefficient of correlation
		CC Get() const {
			CC res;
			if( CCkey::IsS(_ecc) )	res.SetS(_s.Get());
			if( CCkey::IsP(_ecc) )	res.SetP(_p.Get());
			return res;
		}
	};

	struct Range 
	// Structure represented a range.
	// Range is a joint feature with start position and joint (combined) value.
	// All features from bed1, bed2 are splitted by contiguous ranges with predefined value:
	// VAL1 (only the first BedF has a feature here), or
	// VAL2 (only the second BedF has a feature here), or
	// VAL1 & VAL2 (both of the Beds have a feature here), or
	// 0 (no features for both of the Beds)
	{
		chrlen	Start;	// start position of range in standard chromosomal coordinates
		char	Val;	// value of range

		inline Range() : Start(0), Val(0) {}
		inline Range(chrlen pos, char val) : Start(pos), Val(val) {}
	};

	vector<Range>	_ranges;

public:
	// Fills ChromRanges & Range by given two beds.
	// Beds chromosomes should be checked as Treated.
	// Both bed1 & bed2 must be valid: no duplicated, crossed, adjacent, coverage features;
	// in other case R may be wrong
	JointedBeds(BedF & bed1, BedF & bed2);
	
	// Calculates r and fills results
	//	@cc: type of correlation coefficient
	//	Rgns: treated chroms regions
	//	@results: object to fill results
	void CalcR(eCC ecc, GenomeRegions & gRgns, Results & results);
#ifdef DEBUG
	void	Print();
#endif
};

#ifdef DEBUG
class BedMap : public Chroms<arrbmapval>
/*
 * Class 'BedMap' represents bed-file as a list of byte's arrays (one array for one chromosome),
 * where each byte in array represents one nucleotide.
 * Bytes corresponding to the nucleotides included in features have value 1, not included - 0.
 * This is rather slow and simple implementation of calculating algorithm.
 */
{
public:
	//BedMap	(const BedF& bed, const ChromSizes& cSizes);
	BedMap	(const BedF& bed, const GenomeRegions& gRgns);
	
	// Calculates r and fills results
	//	@cc: type of correlation coefficient
	//	@bMap: BedMap object to correlate with
	//	@results: object to fill results
	void CalcR(eCC ecc, BedMap & bMap, Results & results);

	//void  Write(string fileName);
};
#endif	// DEBUG

class CorrPair
/*
 * 'CorrPair' represents pair of objects to compare,
 * and methods for recognizing types and calculation CC
 */
{
private:
	struct FileType	// keeps pointers to the common methods
	{
		typedef Obj*	(CorrPair::*CreateObj)	(const char*, bool);
		typedef void	(*DeleteObj)				(Obj*);	// static function
		typedef bool	(CorrPair::*FillGenRegns)(GenomeRegions&);
		//typedef void	(CorrPair::*CalculateCC)	(GenomeRegions&, Results&);

		const string Ext;				// file name extention
		CreateObj		Create;			// constructor of type
		DeleteObj		Delete;			// destructor of type
		FillGenRegns	FillGenRgns;	// filler GenomeRegions by common chroms
		//CalculateCC		CalcCC;			// calculator
	};
	static const int _FileTypesCnt = 3;		// count of treatment file's types
	static FileType _FileTypes[];

	Obj*	_firstObj;
	Obj*	_secondObj;
	GenomeRegions&	_gRgns;		// initial genome regions to correlate
	//GenomeRegions*	_currgRgns;	// actual genome regions to correlate
	BedF*	_templ;
	chrid	_cID;
	eCC		_ecc;
	bool	_printAlarm;
	bool	_printStat;
	bool	_printTitle;
	bool	_printName;		// print file's names
	BYTE	_typeInd;
	Timer	_timer;

	// Returns true if BedR are treating
	inline bool IsWig()		const { return !_typeInd; }
	// Returns true if BedF are treating
	inline bool IsBedF()	const { return _typeInd == 1; }
	// Returns true if BedE are treating
	inline bool IsAlign()	const { return _typeInd == 2; }

	// Creates features bed object.
	//	@fName: file name
	//	@primary: if true object is primary
	Obj* CreateBedF	(const char* fName, bool primary) {
		BedF* bedF = new BedF(fName, _cID, _printName, primary, _printAlarm, _printStat);
		if( bedF->SameFeaturesLength() ) {
			if( !bedF->EOLPrinted() )	dout << EOL;
			Err("looks like an alignment!", fName).Warning();
		}
		if( primary && !Options::GetIVal(oEXPSTEP) )
			// expand primary bedF f.e. to get true correlation
			// between narrow TFBS (primary) and wide recovered peaks (secondary)
			bedF->Expand(Options::GetIVal(oEXPLEN), _printStat);
		return bedF;
	}
	
	// Creates alignment object
	//	@fName: file name
	//	@primary: if true object is primary
	Obj* CreateBedR	(const char* fName, bool primary) {
		const BedR bedR(NULL, fName, _cID, _printName, Timer::Enabled,
			primary, false, _printStat, Options::GetBVal(oDUPL));
		return new DensMap(bedR, _gRgns);
	}

	// Creates wigMap object.
	//	@fName: file name
	//	@primary: if true object is primary
	Obj* CreateWig	(const char* fName, bool primary) {
		if( _printName ) {
			dout << fName << MSGSEP_BLANK;
			if( primary )	dout << EOL;
		}
		return new WigMap(fName, _cID, primary, false, _gRgns);
	}

	// Fill GenomeRegions by BedF common chroms
	//	return: true if there are comon chroms
	bool FillComnonChromsBedF(GenomeRegions& gRgns);

	// Fill GenomeRegions by coverages or density common chroms
	//	return: true if there are comon chroms
	bool FillComnonChromsMap(GenomeRegions& gRgns);

	// Calculates r for genome features.
	//	@firstBed: first BedF to correlate
	//	@gRgns: GenomeRegions with common chroms in it
	//	@results: results holder
	inline void CalcCCBedF(BedF& firstBed, GenomeRegions& gRgns, Results& results) {
		//((BedF*)_secondObj)->Print();
		JointedBeds(firstBed, *((BedF*)_secondObj)).CalcR(_ecc, gRgns, results);
		// control check by BedMap
		//BedMap bMap(*((BedF*)_secondObj), gRgns);
		//BedMap(firstBed, gRgns).CalcR(_ecc, bMap, results);
	}

	// Calculates r for genome coverages or density.
	//	@gRgns: GenomeRegions with common chroms in it
	//	@results: results holder
	inline void CalcCCMap(GenomeRegions& gRgns, Results& results) {
		((ChromsMap*)_firstObj)->CalcR(_ecc, *((ChromsMap*)_secondObj), gRgns, _templ, results);
	}

	// Checks file extisting and extention validity
	//	@fName: file's name
	//	@abortInvalid: if true throw extention if checking is false
	//	return: index in _FileTypes or _FileTypesCnt if checking is false
	BYTE CheckFileExt(const char * fName, bool abortInvalid);

public:
	// Creates an instance with checking primary object.
	//	@cID: chromosome's ID
	//	@fName: primary file's name
	//	@gRgns: genome regions
	//	@templName: template bed file, or NULL if undefined
	//	@multiFiles: true if more then one secondary files are placed
	CorrPair(chrlen cID,
		const char* fName, GenomeRegions& gRgns,
		const char* templName,
		bool multiFiles
	);

	~CorrPair();
	
	// Adds secondary object, calculates and prints CCkey.
	void CalcCC(const char * fName);
};

#ifdef DEBUG
class TestCC
{
	class Sample
	{
		static const int ArrLen = 100;
		arrchrlen _arr;
	
	public:
		Sample(const string& fname);
	
		inline void GetR(const Sample& sample)
		//{ _arr.GetR(1, eCC::ccSP, sample._arr).Print(); }
		{ _arr.GetR(1, ccP, sample._arr).Print(); }

		void Print()
		{ for(int i=0; i<ArrLen; i++)	cout << _arr[i] << EOL;	}
	};

public:
	TestCC();
};
#endif	// DEBUG
