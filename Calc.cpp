#include "Calc.h"

const char* AND = " and";
const string ForCorrelation = " for correlation";
const string sFormat = " format";

/************************ ChromsMap ************************/

#define	MAP(it)	(it)->second.Map
BYTE ChromsMap::_Space = 0;

struct FeatureR : pair<chrlen, CC>
/*
 * 'FeatureR' keeps feature's ID and R for this feature
 */
{
	inline FeatureR(chrlen i, const CC & res) {	first = i; second = res; }

	inline bool operator < (const FeatureR & rccr) const { 
		return second < rccr.second; 
	}

	// returns single value
	inline double GetSingleVal() const	{ return second.GetSingleVal(); }

	// set single negative value to absolute value
	inline void SetSingleAbsVal()		{ second.SetSingleAbsVal(); } 

	inline void Print() const { dout << first << TAB; second.Print(); }
};

class FeatureRs : vector<FeatureR>
{
private:
	// Replace negative values by positive
	//	return: true if even one value had been replaced
	bool SetAbsVals()
	{
		bool holdNegative = false;
		for(vector<FeatureR>::iterator it=begin(); it!=end(); it++)
			if( it->GetSingleVal() < 0 ) {
				it->SetSingleAbsVal();
				holdNegative = true;
			}
		return holdNegative;
	}

public:
	inline FeatureRs(chrlen cnt)	{ reserve(cnt); }

	inline void AddVal(chrlen ind, CC val) {	push_back(FeatureR(ind+1, val)); }

	void Print(int printFRes) {
		if( printFRes == rsOFF )	return;
		if( printFRes == rsC )		// soretd by feature; are sorted initially
			sort(begin(), end());	// by increase
		dout << "#ftr\tr\n";
		for(vector<FeatureR>::iterator it=begin(); it!=end(); it++)
			it->Print();
	}

	// Creates and prints histogram
	void PrintHist(double binWidth) {
		if( !binWidth )		return;
		SetAbsVals();
		sort(begin(), end());	// by increase
		// define factor
		// factor is a divisor of binWidth: 0.1--0.9=>10, 0.01--0.09=>100 etc
		short F = 10;
		for( ; binWidth * F < 1; F *=10);
		vector<FeatureR>::iterator it=begin();
		// then float instead of double because of wrong consolidation by round double
		double minBin = float(int(F*it->GetSingleVal()))/F;
		//float minBin = F*it->GetSingleVal();
		//		minBin = float(int(minBin))/F;
		double maxBin = F*(end()-1)->GetSingleVal();
		int	maxdecBin = int(maxBin);
		if( maxBin - maxdecBin )	maxdecBin++;	// round up
		if( maxdecBin % 2 )			maxdecBin++;	// get even bin
		maxBin = double(maxdecBin)/F;
		Array<int> hist(int((maxBin - minBin) / binWidth) + 1);		// histogram
		// consolidation: fill histogram
		for(; it!=end(); it++)
			hist[ int((maxBin-it->GetSingleVal()) / binWidth) ]++;
		// print histogram
		//dout << "bin up\tcount\t" << hist.Length() << EOL;
		dout << "bin up\tcount\n";
		for(BYTE k=0; k<hist.Length(); k++)
			dout << (maxBin-k*binWidth) << TAB << hist[k] << EOL;
		//dout << "finish\n";
	}
};

// Gets a pair of total genome means between this and second ChromsMap
//	@gRgn: genome regions
//	@map: second ChromsMap
//	return: a pair of means for each set of chroms
pairDbl ChromsMap::GetGenomeMean(GenomeRegions& gRgn, const ChromsMap& map) const
{
	// calculate total mean through the set of arrays:
	// total_mean = SUM(arr(i).mean * arr(i).relative_count) / SUM(arr(i).relative_length
	//	where:
	//	arr(i).mean - mean of each array i
	//	arr(i).relative_length = arr(i).length / min.arr.length

	chrlen	minSz = gRgn.MinSize(),	// minimal chrom size
			relSz,				// chrom relative size: chrom_size/min_chrom_size
			sumRelSz = 0;		// sum of chrom relative sizes
	double	sumRelMean1 = 0,	// sum of relative length for first object
			sumRelMean2 = 0;	// sum of relative length for second object
	pairDbl	means;

	for(GenomeRegions::cIter it=gRgn.cBegin(); it!=gRgn.cEnd(); it++) {
		relSz = gRgn.Size(it)/minSz;
		means = At(CID(it)).Map.SetMeans(map[CID(it)]);
		sumRelMean1 += means.first	* relSz;
		sumRelMean2 += means.second * relSz;
		sumRelSz += relSz;
	}
	return pairDbl(sumRelMean1 / sumRelSz, sumRelMean2 / sumRelSz);
}

// Calculates r for each region and fills results
//	@cc: type of correlation coefficient
//	@wig: ChromsMap object to correlate with
//	@shGRgns: shell of treated chrom's regions
//	@results: object to fill results
void ChromsMap::CalcRegionsR(CCkey::eCC ecc, const ChromsMap& wig,
	const ShellGenomeRegions& shGRgns, Results& results)
{
	chrid	cID;
	bool	norm = Options::GetBVal(oFNORM);	// true if normalize regions before calc 
	chrlen	rCnt, rLen,
			maxVal,	
			currMaxVal,
			i, currStart,
			regStart, regEnd,
			regLen;
	ChromsMap::cIter cit1, cit2;		// iterators pointing to the chrom's ChromMap of this and wig
	Regions::Iter rit;				// iterator pointing to the chrom's Regions
	ShellGenomeRegions::cIter cit;	// iterator pointing to the chrom's RegionsRange
	arrchrlen	arr1,		// aggregate of regions1
				arr2,		// aggregate of regions2
				arrMax1,	// max values of regions1
				arrMax2;	// max values of regions2

	for(cit1=cBegin(); cit1!=cEnd(); cit1++) {
		if( !TREATED(cit1) )	continue;
		cID = CID(cit1);
		cit = shGRgns.GetIter(cID);
		if( cit == shGRgns.cEnd() ) {	// no cID found: possible only for bedF
			Err("no " + Chrom::TitleName(cID), Template).
				Warning(": skip this " + Chrom::Title + ForCorrelation);
			continue;
		}
		cit2 = wig.GetIter(cID);
		
		rCnt = shGRgns.RegionsCount(cit);
		rLen = shGRgns.RegionsLength(cit); 
		rLen = shGRgns.RegionsLength(cit)/_Space + 1;	// add 1 since arounds by division
		arr1.Reserve(rLen);
		arr2.Reserve(rLen);
		arrMax1.Reserve(rCnt);
		arrMax2.Reserve(rCnt);
		FeatureRs fResults(_binWidth ? rCnt : 0);	// create histogram
		currStart = maxVal = 0;
		if( norm )
			// normilize data within regions
			// 1. get max values
			for(rit=shGRgns.ChromBegin(cit), i=0; rit!=shGRgns.ChromEnd(cit); rit++, i++) {
				regStart = rit->Start/_Space;
				regEnd = rit->End/_Space;
				if( maxVal < (currMaxVal = arrMax1[i] = MAP(cit1).GetMaxVal(regStart, regEnd)) )
					maxVal = currMaxVal;
				if( maxVal < (currMaxVal = arrMax2[i] = MAP(cit2).GetMaxVal(regStart, regEnd)) )
					maxVal = currMaxVal;
		}
		// concatenate subarrays
		bool fillfResults = _binWidth || _printFRes >= 0;
		for(rit=shGRgns.ChromBegin(cit), i=0; rit!=shGRgns.ChromEnd(cit); rit++, i++) {
			regStart = rit->Start/_Space;
			regLen = rit->Length()/_Space;
			arr1.Concat(MAP(cit1), currStart, regStart, regLen);
			arr2.Concat(MAP(cit2), currStart, regStart, regLen);
			if( norm ) {
				// 2. multiply regions by ratio to normilize
				arr1.MultiplyVal(float(maxVal)/arrMax1[i], currStart, currStart+regLen);
				arr2.MultiplyVal(float(maxVal)/arrMax2[i], currStart, currStart+regLen);
				//arrchrlen::SynchMultiplyVal(arr1, arr2, 
				//	float(maxVal)/arrMax1[i], float(maxVal)/arrMax2[i],
				//	currStart, currStart+regLen);
			}
			if( fillfResults )		// fill results for each region
				fResults.AddVal(i, arr1.GetR(cID, ecc, arr2, currStart, currStart+regLen));
			currStart += regLen;
		}
		results.AddVal(cID, arr1.GetR(cID, ecc, arr2));
		fResults.Print(_printFRes);
		fResults.PrintHist(_binWidth);	// print histogram
	}
}

// Calculates r and fills results
//	@ecc: type of correlation coefficient
//	@wMap: ChromsMap object to correlate with
//	@gRgns: real chrom's regions
//	@bedF: template with defined regions or NULL
//	@results: object to fill results
void ChromsMap::CalcR(CCkey::eCC ecc, const ChromsMap& wig,
	GenomeRegions& gRgns, const BedF* bedF, Results& results)
{
		if( bedF ) {						// calculate R by regions from bedF
			const ShellGenomeRegions shGRgns(*bedF);
			CalcRegionsR(ecc, wig, shGRgns, results);
		}
		else if( !gRgns.SingleRegions() ) {	// calculate R by regions from gRgns
			const ShellGenomeRegions shGRgns(gRgns);
			CalcRegionsR(ecc, wig, shGRgns, results);
		}
		else {
			Iter it;
			if( results.GiveLocal() )		// calculate R for each whole chrom
				for(it=Begin(); it!=End(); it++)
					if( TREATED(it) )
						results.AddVal(CID(it), MAP(it).GetR(CID(it), ecc, wig[CID(it)]));
		
			if( results.GiveTotal() ) {		// calculate total R
				CCaggr ccaggr;
				if( CCkey::IsP(ecc) )
					ccaggr.SetMeans(GetGenomeMean(gRgns, wig));
				for(it=Begin(); it!=End(); it++)	// accumulate ccaggr
					if( TREATED(it) ) {
						MAP(it).AccumVars(CID(it), ecc, ccaggr, wig[CID(it)]);
						if( ccaggr.IsUndefS() )		// check for esceeding
							break; 
					}
				results.SetTotal(ccaggr.GetR());		// calculate total R
			}
		}
}

#ifdef DEBUG
// Prints regions.
//	@rgnCnt: number of regions to print or all if not specified
void ChromsMap::Print(chrlen rgnCnt) const
{
	chrlen i, k=0, pos=0;
	chrlen val=0, newVal;

	if( !rgnCnt )	rgnCnt = CHRLEN_UNDEF;
	cout << "ChromsMap:\n";
	for(cIter it = cBegin(); it != cEnd(); it++) {
		cout << Chrom::AbbrName(CID(it)) << EOL;
		const arrchrlen& map = MAP(it);
		for(i=0; i<map.Length(); i++) {
			newVal = map[i];
			if( val != newVal ) {
				if( k >= rgnCnt )	break;
				cout << pos << TAB << val << EOL;
				val = newVal;
				pos = i;
				k++;
			}
		}
	}
}
#endif

/************************ end of class ChromsMap ************************/

/************************ class WigMap ************************/

const char* kyeTrack	= "track type=";
const char* kyeWiggle	= "wiggle_0";
const char* keyVarStep	= "variableStep";
const char* keyFixStep	= "fixedStep";
const char* keyChrom	= "chrom=";
const char* keySpan		= "span=";
const char* keySpace	= "space=";
const char* progSpec	= "regulated";
const char* sPR			= "PeakRanger";

const BYTE	lenKeyVarStep = strlen(keyVarStep);
const BYTE	lenKeyChrom = strlen(keyChrom);
const BYTE	lenKeySpan = strlen(keySpan);

// Adds region for current chromosome and increases record counters.
//	@start:	region's start position
//	@size:	region's length
//	@val:	region's value
//	@ambig: ambiguities to collect statistics
//	return: 1 if region is added; otherwise 0
void WigMap::WigPocket::AddRegion(chrlen start, chrlen size, chrlen val, Ambig& ambig)
{
	if( !val )	return;
	_recCnt.first++;
	chrlen end = start + size;
	if( end > _cSize ) {
		// !! in case of regions _cSize is not a chrom size, but end of last region. Надо бы сделать 2 раздельных
		ambig.TreatCase(Ambig::EXCEED);
		return;							// only omit exceeding record
	}
	if( size < ChromsMap::_Space )	{
		start = AlignPos(start, _Space, 1);
		end = AlignPos(end, _Space, 1);
	}
	// reduce start and end because wig chroms positions are 1-relative
	_map->Fill(--start/_Space, --end/_Space, val);
	_recCnt.second++;
}


// Returns a pointer to the substring defined by key.
//	@str: null-terminated string to search
//	@key: null-terminated string to search for
//	return: a pointer to the substring after key, or NULL if key does not appear in str
const char* KeyStr(const char* str, const char* key)
{
	const char* strKey = strstr(str, key);
	return strKey ? (strKey + strlen(key)) : NULL;
}

// Checks definition or declaration line for key
//	return: point to substring followed after the key
inline const char* CheckSpec(const char* line, const char* key, const TabFile& file)
{
	const char* strKey = KeyStr(line, key);
	if( !strKey )
		file.ThrowLineExcept(string("absent or wrong '") + key + "' key");
	return strKey;
}

// Returns span value
//	@str: null-terminated string to search span
//	return: span value, or 1 if span does not appear in str
chrlen GetSpan(const char* str, const TabFile& file)
{
	return *str == keySpan[0] ?		// span is optional key
		atoi(CheckSpec(str, keySpan, file)) :
		1;
}

// Initializes instance from tab file
//	@ambig: ambiguities
//	@pgRgns: pointer to GenomeRegions inctance
//	return: numbers of all and initialied items for given chrom
dchrlen WigMap::InitChild	(Ambig& ambig, void* pgRgns)
{
	ULONG	pos;				// position of current line. ULONG, since it used in GetFirstLine()
	chrlen	newPos = 0,			// positions of new line
			startPos = 0,		// current region position
			currSpan,			// current span (count of data with the same value
			startSpan = 0;		// span from last declaration line
								// for check in AddRegion() only
	chrlen	val = 0, newVal;	// current, new readed values. Should be float ??
	chrid	cID = Chrom::UnID,	// current chrom ID
			newcID;				// chrom ID from declaration line
	bool	skipLine = false;	// if true skip data line
	const char* line;			// current readed line
	TabFile&	file = ambig.File();
	WigPocket	pocket(*(GenomeRegions*)pgRgns);
		
	if( !(line = file.GetFirstLine(&pos)) )	return make_pair(0, 0);
	Reserve(Chrom::StatedAll() ? Chrom::Count : 1);
	pos = 0;
	line = CheckSpec(line, kyeTrack, file);		// check track type key
	currSpan = strchr(line, BLANK) - line;		// temp using: the length of wiggle type in definition
	if( strncmp(line, kyeWiggle, currSpan) )	// not a wiggle_0.  use _stricmp ?
		file.ThrowExcept("type '" + string(line, currSpan) + "' does not supported");
	if( KeyStr(line, sPR) && !KeyStr(line, progSpec) )	// check regulated PR
		file.ThrowExcept("unregulated " + string(sPR) + " wiggle");
	if( file.IsAborting() ) {				// define static space
		line = KeyStr(line, keySpace);
		if(line)	_Space = atoi(line);
	}
	while( line = file.GetLine() )
		if( isdigit(line[0]) ) {	// ** data line
			if( skipLine )	continue;
			newPos = file.IntField(0);
			newVal = chrlen(file.IntField(1));
			if( val == newVal && newPos-pos == _Space)
				currSpan += _Space;		// unregulated (fixed span): accumulate current span
			else {
				pocket.AddRegion(startPos, currSpan, val, ambig);
				currSpan = startSpan;
				startPos = newPos;
				val = newVal;
			}
			pos = newPos;
		}
		else {							// ** declaration line
			if( !KeyStr(line, keyVarStep) )		// not started from 'variableStep'?
				if( KeyStr(line, keyFixStep) )	// fixedStep
					file.ThrowExcept(string(keyFixStep) + " is not acceptable");	// wrong type
				else	// wrong or not stated; double KeyStr() call, but in uniformly order
					CheckSpec(line, keyVarStep, file);
			// add last region if val > 0
			pocket.AddRegion(startPos, currSpan, val, ambig);
			startSpan = val = 0;
			// define chromosome
			line = CheckSpec(line, keyChrom, file);
			newcID = Chrom::IDbyAbbrName(line);

			if( newcID == Chrom::UnID )	continue;	// skip additional chroms
			// set chromosome
			line += strlen(Chrom::Abbr) + Chrom::NameLength(newcID) + 1;	// stay at "span="
			if( cID != newcID ) {	// new chromosome
				cID = newcID;
				if( skipLine = (Chrom::StatedAll()		// read all chromosomes
				|| Chrom::StatedID() == cID) ) {		// read given chromosomes
					if( !_Space )	// static space wasn't defined in definition line
						_Space = startSpan = currSpan = GetSpan(line, file);
					AddChrom(cID, pocket);	// _Space should be defined at this time
					ambig.SetTreatedChrom(cID);
					pos = 0;
				}
				else if(ChromsCount())	break;	// given is readed already
				skipLine = !skipLine;
			}
			if( !(skipLine || startSpan) )		// define current span
				startSpan = currSpan = GetSpan(line, file);
			startPos = pos;
		}	
	pocket.AddRegion(startPos, currSpan, val, ambig);	// add region at last line

	return pocket.RecCount();
}

/************************ end of class WigMap ************************/

/************************ class DensMap ************************/

DensMap::DensMap(const BedR& bedR, GenomeRegions& gRgns)
{
	_Space = Options::GetIVal(oSPACE);
	DensPocket pocket(bedR, gRgns);
	Reserve(bedR.ChromsCount());
	for(BedR::cIter it=bedR.cBegin(); it!=bedR.cEnd(); it++) {
		AddChrom(it, pocket);
		pocket.ScanChrom();
	}
	_isBad = bedR.IsBad();
	_EOLneeded = bedR.EOLNeeded();
}

// Adds chrom by BedR iter
//	@it: BedR iter
//	@map: pointer to the chroms map
void DensMap::DensPocket:: AddChrom (BedR::cIter it, arrchrlen* map)
{
	Reserve(CID(it), map);
	_currIt = _bedR.ReadsBegin(it);
	_endIt	= _bedR.ReadsEnd(it);
	_wi = 0;
	_wCurrLen = _Space;
}

// Fills number or reads from window started with @start position
void DensMap::DensPocket::ScanWindow(chrlen start)
{
	chrlen	rCentre;
	ULONG rCnt = 0;

	for(; _currIt!=_endIt; _currIt++) {		// loop through window for Reads
		//rCentre = *_currIt + _halfRLen;
		rCentre = _currIt->Pos + _halfRLen;
		if( rCentre >= start ) {				// pass left-of-window Reads
			if( rCentre >= start + _wCurrLen )
				break;					// exit on first right-of-window Read
			rCnt++;
			//_rtotalCnt++;
		}
	}
	if( rCnt > UINT_MAX )
		Err("Density value " + NSTR(rCnt) + " exceeded limit or " + NSTR(UINT_MAX)
		+ ". Try to reduce space").Throw();
	(*_map)[_wi] = chrlen(rCnt);
}

// Fills number or reads from region @rgn
//void DensMap::DensPocket::ScanRegion(const Region& rgn)
void DensMap::DensPocket::ScanChrom()
{
	//chrlen	start = rgn.Start,
	//		end = rgn.End;
	chrlen start = 0, end = _cSize;

	if(start + _wCurrLen <= end) {		// current window belong to region entirely
		do {
			ScanWindow(start);
			_wi++;
			//if( ++_wi == static_cast<chrlen>(_map->Length()) )//*_Space) )
			//	Err("window's index "+NSTR(_wi)+" is out of range "+NSTR(_map->Length()),
			//		"ReadDistrib::ScanRegion").Warning();
			start += _wCurrLen;
			_wCurrLen = _Space;			// set user's win length
		}
		while(start + _wCurrLen < end);
		if(start + _wCurrLen >= end)	// window's and region's ends are equal or output
			return;
	}
	// current window exceeds region
	_wCurrLen -= end - start;			// decrease current window by the rest of region
	ScanWindow(start);					// scan part of current window
}

/************************ end of class DensMap ************************/

/************************ class Results ************************/

// Prints results for each chromosome and total or total only (if it was defined).
//	@printTitles: if true, print titles, otherwise doesn't only in case of single result
void Results::Print(bool printTitles)
{
	cIter it = cBegin();
	chrid cCnt = ChromsCount();
	if( !cCnt && !_total.NotEmpty())		
		dout << " no " << Chrom::TitleName() << ForCorrelation << EOL;
	else if( cCnt == 1  ) {
		if( printTitles )	
			dout << Chrom::AbbrName(CID(it)) << TAB;
		it->second.Print();
	}
	else {
		Sort();
		for(; it!=cEnd(); it++) {
			dout << Chrom::AbbrName(CID(it)) << TAB;
			it->second.Print();
		}
		if( _total.NotEmpty() ) {
			dout << Total << SepSCl;
			_total.Print();
		}
	}
}

/************************ end of class Results ************************/

/************************ class JointedBeds ************************/

#define VAL1	0x1	// value represented first BedF's feature
#define VAL2	0x2	// value represented second BedF's feature

// Fills ChromRanges & Range by given two beds.
// Beds chromosomes should be checked as Treated.
// Both bed1 & bed2 must be valid: no duplicated, crossed, adjacent, coverage features;
// in other case R may be wrong
JointedBeds::JointedBeds(BedF& bed1, BedF& bed2)
{
	char val;			// current joint range value
	chrlen	pos, pos2,	// current positions and in bed2. Initially are equal
			fi1, fi2,				// features indexes in bed1, bed2
			fCnt1, fCnt2,			// count of features in bed1, bed2
			firstInd=0, lastInd=0;	// current first, last feature indexes in JointedBeds
	Region	f1, f2,					// dedicated feature used for detecting
			fEnd = Region(CHRLEN_UNDEF, CHRLEN_UNDEF-1);// last chromosome's joint feature
	//chrid	cID;
	Bed::cIter cit1, cit2;

	//bed2.Print();
	Reserve(min(bed1.ChromsCount(), bed2.ChromsCount()));
	_ranges.reserve( 2*(bed1.FeaturesCount() + bed2.FeaturesCount()));	// ranges

	for(cit1 = bed1.Begin(); cit1 != bed1.End(); cit1++) {
		if( !cit1->second.Treated )	continue;
		fi1 = fi2 = val = 0;
		//cID = bed1.ChromID(cit1);
		cit2 = bed2.GetIter(CID(cit1));
		fCnt1 = cit1->second.ItemsCount();
		fCnt2 = cit2->second.ItemsCount();
		f1 = bed1.Feature(cit1);
		f2 = bed2.Feature(cit2);
		// loop through current chromosome's features 
		while( fi1 < fCnt1 || fi2 < fCnt2 ) {
			pos = val & VAL1 ? (f1.End + 1) : f1.Start;
			pos2 = val & VAL2 ? (f2.End + 1) : f2.Start;
			if( pos < pos2 ) {
				val ^= VAL1;		// flip val for bed1
				if( !(val & VAL1) )	// true when bed1 feature is closed (every second range)
					f1 = ++fi1 < fCnt1 ? bed1.Feature(cit1, fi1) : fEnd;
			}
			else if( pos > pos2 ) {
				pos = pos2;
				val ^= VAL2;		// flip val for bed2
				if( !(val & VAL2) )	// true when bed2 feature is closed (every second range)
					f2 = ++fi2 < fCnt2 ? bed2.Feature(cit2, fi2) : fEnd;
			}
			else {
				val ^= VAL1 ^ VAL2;	// flip val for both beds
				if( !(val & VAL1) )	// true when bed1 feature is closed 
					f1 = ++fi1 < fCnt1 ? bed1.Feature(cit1, fi1) : fEnd;
				if( !(val & VAL2) )	// true when bed2 feature is closed 
					f2 = ++fi2 < fCnt2 ? bed2.Feature(cit2, fi2) : fEnd;
			}
			_ranges.push_back( Range(pos, val) );	// add new joint feature
			lastInd++;
		}
		AddVal(CID(cit1), ChromRanges(
			firstInd,
			lastInd,
			bed1.FeaturesLength(cit1),
			bed2.FeaturesLength(cit2)
			));
		firstInd = lastInd;
	}
}

// Calculates r and fills results
//	@cc: type of correlation coefficient
//	@cSizes: chrom sizes
//	@results: object to fill results
void JointedBeds::CalcR(CCkey::eCC ecc, const ChromSizes& cSizes, Results& results)
{
	genlen	gSize;			// genome's size
	chrlen	cSize,			// chromosome's size
			start, stop,	// range's boundary positions
			ri,				// index of range
			len;			// length of range
	char	val;			// value of current range
	double	featrsLen1, featrsLen2;
	PairR	localCC(ecc),
			totalCC(ecc);
	ChromRanges cRanges;	// current ChromRanges

	if( results.GiveTotal() )
		gSize = cSizes.GenSize();
	for(cIter it=cBegin(); it!=cEnd(); it++) {
		cRanges = it->second;
		cSize = cSizes[CID(it)];
		featrsLen1 = double(cRanges.FeatrsLen1);
		featrsLen2 = double(cRanges.FeatrsLen2);

		if( results.GiveLocal() )
			localCC.Init(featrsLen1/cSize, featrsLen2/cSize, true);
		if( results.GiveTotal() )
			totalCC.Init(featrsLen1/gSize, featrsLen2/gSize, false);
		
		start = val = 0;	// first range
		for(ri=cRanges.FirstInd; ri<=cRanges.LastInd; ri++) {
			stop = _ranges[ri].Start;
			len = stop - start;
			if( results.GiveLocal() )	localCC.Increment (len, val);
			if( results.GiveTotal() )	totalCC.Increment (len, val);
			// next range
			val = _ranges[ri].Val;
			start = stop;
		}
		len = cSize - stop;		// last range
		if( results.GiveLocal() ) {
			localCC.Increment (len, 0);
			results.AddVal(CID(it), localCC.Get());
		}
		if( results.GiveTotal() )
			totalCC.Increment(len, 0);
	}
	if( results.GiveTotal() )
		results.SetTotal( totalCC.Get() );
}

#ifdef DEBUG
void	JointedBeds::Print()
{
	chrlen	ri;				// index of range
	cout << "JointedBeds:\n";
	for(cIter it=cBegin(); it!=cEnd(); it++) {
		cout << Chrom::AbbrName(CID(it)) << ":";
		for(ri=it->second.FirstInd; ri<=it->second.LastInd; ri++)
			cout << '\t' << _ranges[ri].Start << '\t' << int(_ranges[ri].Val) << EOL;
	}
}
#endif

// Keeps mean values for bed1, bed2. 
//	@clear: if true, clear instance for treatment of new chromosome
void JointedBeds::PairR::R::Init(double mean1, double mean2, bool clear)
{
	double	d1 = 1 - mean1,
			d2 = 1 - mean2;

	_sqMeans1[0] = _sqMeans1[2] = mean1 * mean1;
	_sqMeans1[1] = _sqMeans1[3] = d1 * d1;
	_sqMeans2[0] = _sqMeans2[1] = mean2 * mean2;
	_sqMeans2[2] = _sqMeans2[3] = d2 * d2;
	_crossMeans[0] = mean1 * mean2;
	_crossMeans[1] = -mean2 * d1;
	_crossMeans[2] = -mean1 * d2;
	_crossMeans[3] = d1 * d2;
	if( clear )
		_cov = _var1 = _var2 = 0;
}

// Accumulates next length of range
void JointedBeds::PairR::R::Increment(chrlen len, char val)
{
	_cov	+= len * _crossMeans[val];
	_var1	+= len * _sqMeans1[val];
	_var2	+= len * _sqMeans2[val];
}

/************************ end of class JointedBeds ************************/

#ifdef DEBUG
/************************ class BedMap ************************/

#define YES	1
#define NO	0

//BedMap::BedMap(const BedF& bed, const ChromSizes& cSizes)
BedMap::BedMap(const BedF& bed, GenomeRegions& gRgns)
{
	chrid cID;
	chrlen i, fCnt;

	Reserve(bed.ChromsCount());
	for(Bed::cIter it=bed.cBegin(); it!=bed.cEnd(); it++)
		if( TREATED(it) ) {
			cID = CID(it);
			fCnt = bed.FeaturesCount(it);
			//arrbmapval map(cSizes[cID]);
			arrbmapval map(gRgns.Size(cID));
			for(i = 0; i < fCnt; i++) {
				const Region & rgn = bed.Feature(it, i);
				map.Fill(rgn.Start, rgn.End, YES);
			}
			AddClass(cID, map);
		}
}

// Calculates r and fills results
//	@cc: type of correlation coefficient
//	@bMap: BedMap object to correlate with
//	@results: object to fill results
void BedMap::CalcR	(CCkey::eCC ecc, BedMap & bMap, Results & results)
{
	for(Iter it=Begin(); it!=End(); it++)
		results.AddVal(CID(it), it->second.GetR(CID(it), ecc, bMap.At(CID(it))));
}

//void BedMap::Write(string fileName)
//{
//	TxtFile oFile(fileName, TxtFile::WRITE, 1);
//	long i,
//		size = _maps.Length(),
//		cntLines =  size / LINE_LEN;
//	char cntChroms =  (char)_maps.Length();
//	for(char k=0; k<cntChroms; k++) {
//		for(i = 0; i<cntLines; i++)
//			oFile.AddRecord(_maps[k].Map + i * LINE_LEN, LINE_LEN);
//		if( size % LINE_LEN )		// last line
//			oFile.AddRecord(_maps[k].Map + i * LINE_LEN, size%LINE_LEN);
//	}
//	oFile.Write();
//}

/************************ end of class BedMap ************************/
#endif	// DEBUG

/************************ class CorrPair ************************/

static inline void DelWig (Obj* obj) { if(obj) delete (WigMap*)obj;	}
static inline void DelBedF(Obj* obj) { if(obj) delete (BedF*)obj;	}
static inline void DelBedR(Obj* obj) { if(obj) delete (DensMap*)obj;}

CorrPair::FileType CorrPair::_FileTypes[_FileTypesCnt] = {
	{ &CorrPair::CreateWig,  DelWig, &CorrPair::FillComnonChromsMap	},
	{ &CorrPair::CreateBedF, DelBedF,&CorrPair::FillComnonChromsBedF}, 
	{ &CorrPair::CreateBedR, DelBedR,&CorrPair::FillComnonChromsMap,}
};

// Creates an instance with checking primary object.
//	@cID: chromosome's ID
//	@primefName: primary file's name
//	@cSizes: chrom sizes to check excdeeding chrom length
//	@gRgns: genome regions
//	@templName: template bed file, or NULL if undefined
//	@multiFiles: true if more then one secondary files are placed
CorrPair::CorrPair(
	const char* primefName, const ChromSizes* cSizes, GenomeRegions& gRgns,
	const char* templName, bool multiFiles
) :
	_firstObj(NULL), _secondObj(NULL), _templ(NULL),
	_cSizes		(cSizes),
	_gRgns		(gRgns),
	_ecc		(CCkey::eCC(Options::GetIVal(oCC))),
	_info		(Bed::eInfo(Options::GetIVal(oINFO))),
	_typeInd	(CheckFileExt(primefName, true)),
	_printWarn	(Options::GetBVal(oWARN)),
	_absolPrintName	(multiFiles)
{
	if(templName)
		if(IsBedF()) {
			if(IsNotLac())		Err("ignored", string(Template) + sBLANK + templName).Throw(false);
		}
		else {
			_templ = new BedF(Template, FS::CheckedFileName(templName),
				_cSizes, _info, _absolPrintName, true, _printWarn);
			if(!_templ->Extend(Options::GetIVal(oEXTLEN), _cSizes, _info))
				_templ->CheckFeaturesLength(Options::GetIVal(oSPACE), "space", Template);
		}
	if( IsNotLac() ) {
		if( CCkey::IsP(_ecc) )		dout << "Pearson";
		if( CCkey::IsBoth(_ecc) )	dout << AND << BLANK;
		if( CCkey::IsS(_ecc) )		dout << "Signal";
		dout << " CC between ";
	}
	_firstObj = (this->*_FileTypes[_typeInd].Create)(primefName, true);
	if( _absolPrintName || IsNotLac() ) {
		dout << AND;
		if( multiFiles )	dout << "...";
		dout << EOL;
	}
}

CorrPair::~CorrPair() {
	if(_templ)		delete _templ;
	//if(_currgRgns)	delete _currgRgns;
	_FileTypes[_typeInd].Delete(_firstObj);
	_FileTypes[_typeInd].Delete(_secondObj);
}

// Adds secondary object, calculates and prints CCkey.
void CorrPair::CalcCC(const char* fName)
{
	_FileTypes[_typeInd].Delete(_secondObj);
	_secondObj = NULL;

	// check type
	BYTE typeInd = CheckFileExt(fName, false);
	if( typeInd == _FileTypesCnt )	return;
	if( typeInd != _typeInd) {
		Err("different" + sFormat, fName).Throw(false);
		return;
	}
	// create object
	_secondObj = (this->*_FileTypes[_typeInd].Create)(fName, false);
	if( _secondObj->IsBad() )	return;
	// set common (treated) chroms and regions
	GenomeRegions gRgns(_gRgns);	// actually treated chroms regions
	if( !(this->*_FileTypes[_typeInd].FillGenRgns)(gRgns) )	return;

	// calculate r
	Results results(Options::GetIVal(oPRCC));
	if( IsBedF() ) {
		int expStep	= Options::GetIVal(oEXTSTEP);
		if(expStep) {	// calculation r by step increasing expanding length
			int expLen	= Options::GetIVal(oEXTLEN);
			for(int i=0; i<=expLen; i += expStep) {
				BedF corrBedF(*((BedF*)_firstObj));
				corrBedF.Extend(i, _cSizes, _info);
				CalcCCBedF(corrBedF, results);
				//dout << i << TAB;
				results.Print(IsNotLac());
				results.Clear();
			}
			return;
		}
		else			// primary bedF is already extended
			CalcCCBedF(*((BedF*)_firstObj), results);
	}
	else
		CalcCCMap(gRgns, results);
	results.Print(IsNotLac());
}

// Creates features bed object.
//	@fName: file name
//	@primary: if true object is primary
Obj* CorrPair::CreateBedF	(const char* fName, bool primary)
{
	BedF* bedF = new BedF(NULL,fName, _cSizes, _info, _absolPrintName, primary, _printWarn);
	// if primary is bad, it throws exception and quit
	if( !bedF->IsBad() && bedF->SameFeaturesLength() ) {
		if( bedF->EOLNeeded() )	dout << EOL;
		Err("looks like an alignment!", fName).Warning();
	}
	if(primary)
		// expand primary bedF f.e. to get true correlation
		// between narrow TFBS (primary) and wide recovered peaks (secondary)
		bedF->Extend(Options::GetIVal(oEXTLEN), _cSizes, _info);
	return bedF;
}
	
// Creates alignment object
//	@fName: file name
//	@primary: if true object is primary
Obj* CorrPair::CreateBedR	(const char* fName, bool primary)
{
	const BedR bedR(NULL, fName, _cSizes, _info, _absolPrintName, primary, _printWarn,
		Options::GetBVal(oDUPL));
	return new DensMap(bedR, _gRgns);
}

// Creates wigMap object.
//	@fName: file name
//	@primary: if true object is primary
inline Obj* CorrPair::CreateWig	(const char* fName, bool primary)
{
	return new WigMap(fName, _cSizes, _info, _absolPrintName, primary, _gRgns);
}

// Checks first & second bedF for common chroms.
//	@gRgns: GenomeRegions: not used
//	return: true if there are common chroms
bool CorrPair::FillComnonChromsBedF(GenomeRegions& gRgns)
{
	BedF& obj1 = *((BedF*)_firstObj);

	if( !SetCommonChroms(obj1, *((BedF*)_secondObj), _printWarn) )
	{	Err(sNoCommonChroms).Throw(false); return false; }
	//const Regions rgns;
	//for(BedF::cIter it=obj1.cBegin(); it!=obj1.cEnd(); it++)
	//	if( TREATED(it) )
	//		gRgns.AddChrom(CID(it), rgns);
	return true;
}

// Checks first & second ChromsMap for common chroms
// and fill GenomeRegions by coverages or density common chroms
//	@gRgns: GenomeRegions to fill
//	return: true if there are common chroms
bool CorrPair::FillComnonChromsMap(GenomeRegions& gRgns)
{
	ChromsMap& obj1 = *((ChromsMap*)_firstObj);
	
	// set common (treated) chroms and regions
	if( !SetCommonChroms(obj1, *((ChromsMap*)_secondObj), _printWarn) )
	{	Err(sNoCommonChroms).Throw(false); return false; }
	for(ChromsMap::cIter it=obj1.cBegin(); it!=obj1.cEnd(); it++)
		if( TREATED(it) )
			gRgns.AddChrom(CID(it), _gRgns[CID(it)]);
	return true;
}

// Checks file extisting and extention validity
//	@fName: file's name
//	@abortInvalid: if true throw extention if checking is false
//	return: index in _FileTypes or _FileTypesCnt if checking is false
BYTE CorrPair::CheckFileExt(const char * fName, bool abortInvalid)
{
	if( FS::CheckFileExist(fName, abortInvalid) )	return _FileTypesCnt;
	FT::eTypes format = FT::GetType(fName);
	if(format==FT::WIG)	return 0;
	if(format==FT::BED)	return Options::GetBVal(oALIGN) ? 2 : 1;
	Err("unpredictable" + sFormat, fName).Throw(abortInvalid);
	return _FileTypesCnt;
}

/************************ end of class CorrPair ************************/

#ifdef DEBUG
TestCC::Sample::Sample(const string& fname)
{
	_arr.Reserve(ArrLen);
	TabFile file(fname);

	for(int i=0; file.GetLine() || i<ArrLen; i++ )
		_arr[i] = file.IntField(0);
}

TestCC::TestCC()
{
	const string path = "..\\..\\TestR\\";
	Sample sample1(path + "testR1.txt");
	Sample sample2(path + "testR2.txt");
	sample1.GetR(sample2);
}
#endif	// DEBUG
