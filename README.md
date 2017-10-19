# bioCC
Fast advanced **C**orrelation **C**alculator for basic **bio**informatic file formats.<br>
It computes Pearson’s and signal’s correlation coefficients for densities, coverage and features.<br>
Program allows to obtain correlation coefficients for the whole genome, for each chromosome separately and for predefined regions within the chromosomes.

**bioCC** is designed to treat a bunch of files at once.

The program runs on the command line under Linux and Windows.

## Installation
### Executable file

**Linux**<br>
Go to the desire directory and type commands:<br>
```wget -O bioCC.gz https://github.com/fnaumenko/bioCC/releases/download/1.0/bioCC-Linux-x64.gz```<br>
```gzip -d bioCC.gz```<br>
```chmod +x bioCC```

**Windows**<br>
Download archive from [here](https://github.com/fnaumenko/bioCC/releases/download/1.0/bioCC-Windows-x64.zip) 
and unzip by any archiver, for instance [WinRar](https://www.win-rar.com/download.html?&L=0).

### Compiling in Linux
Required libraries:<br>
g++<br>
zlib (optionally)

Go to the desired directory and type commands:<br>
```wget -O bioCC.zip https://github.com/fnaumenko/bioCC/archive/1.0.zip```<br>
```unzip bioCC.zip```<br>
```cd bioCC-1.0```<br>
```make```

If **zlib** is not installed on your system, a message will be displayed from the linker.<br>
In that case you can compile the program without the ability to work with .gz files. 
To do this, open *makefile* in any text editor, uncomment last macro in the second line, comment third line, save *makefile*, and try again ```make```.<br>
To be sure about **zlib** on your system, type ```whereis zlib```.

## Usage
```
bioCC [options] -g|--gen <name> file0 file1 ...
bioCC [options] -g|--gen <name> -l|--list <file>
```

### Help
```
Input:
  -a|--align            input bed files are alignments
  -g|--gen <name>       chromosome sizes file, reference genome library, or single nucleotide sequence. Required
  --gap-len <int>       minimal length of undefined nucleotide region in genome which is declared as a gap.
                        Ignored for the chromosome sizes file and for the ordinary beds [1000]
  -d|--dupl <OFF|ON>    accept duplicate reads. For the alignments only [ON]
  -l|--list <name>      list of multiple input files.
                        First (primary) file in list is comparing with others (secondary)
Processing:
  -c|--chr <chars>      treat specified chromosome only
  -r|--cc <P,S>         correlation coefficient, in any combination: P - Pearson, S - signal [P]
  -s|--space <int>      resolution: span in bp by which reads will be counted to define a density.
                        For the alignments only [100]
  -p|--pr-cc <IND,TOT>  print coefficient, in any combination:
                        IND - for each chromosome individually, TOT - total [IND]
Region processing:
  -f|--fbed <name>      'template' ordinary bed file which features define compared regions.
                        Ignored for the ordinary beds
  -e|--ext-len <int>    length by which the features in primary file (for ordinary beds) or in 'template'
                        (for alignments and wigs) will be extended in both directions before treatment [0]
  --ext-step <int>      step of extending features in primary bed file; if 0 then no step calculation.
                        For the ordinary beds only [0]
  -b|--bin-width <float>        width of the histogram bin [0]
  --sort <RGN|CC>       print region coefficients, sorted by: RGN - regions, CC - coefficients
  --norm <OFF|ON>       normalize regions before calculation. Ignored for the ordinary beds [ON]
Output:
  -i|--info <LAC|NM|CNT|STAT>   print information about file:
                        LAC - laconic, NM - name only, CNT - number of items, STAT - statistics [NM]
  -w|--warn             print each file's item ambiguity, if they exist.
  -o|--out              duplicate standard output to bioCC_out.txt file
Other:
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit
```

## Details

### Input data types
Alignment read densities should be presented by [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format.<br>
Sequencing coverage should be presented by [WIG](https://genome.ucsc.edu/goldenpath/help/wiggle.html) format.<br>
Features are also presented by [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format.

**BED**<br>
Formally bed file as a result of peak calling process, and bed file represents alignment have the common required fields, but they differ in the interpretation. 
Therefore, for greater clarity, we will refer to the files of the first type as *ordinary bed*, and the second – *alignment bed*.<br>
For more information about read’s density see option ```-s|--space```.

**WIG**<br>
*bedGraph* type is not supported in this version, and the only *wiggle* in *variableStep* format is acceptable.<br>
Note, that one of the fastest *wiggle* generator [PeakRanger](http://ranger.sourceforge.net/manual1.18.html) 
(it also supports the separate strand generation) produces data with some peculiarity, 
namely each its interval has initial span = 1. 
This does not allow to compare coverage.<br>
To level this feature, use the [**wigReg**](https://github.com/fnaumenko/wigReg) software.<br>
**wigReg** also allows to compress *wig* files generated by [MACS](http://liulab.dfci.harvard.edu/MACS/00README.html) 2-10 times.

All input files should have chromosomes items (features, data lines) been clustering together.<br>
For faster processing, items belonging to the same chromosome also should be sorted in position ascending order.<br>
*Wig* files usually meet these requirements. 
Opposite, *bed* files often have messy initialization. 
The simplest way to prepare *bed* files is to sort them, for example by using **sortBed** utility from [bedtools](http://bedtools.readthedocs.io/en/latest/) package 
(though for **bioCC** the order of the chromosomes themselves does not matter).<br>
It concerns 'template' *bed* file as well (see ```-f|--fbed``` option).

**bioCC** recognizes the file formats automatically by their extension. To distinguish between *ordinary beds* and *alignments*, a special option is provided.

Compressed files in gzip format (.gz) are acceptable.

### Input data order
Comparable files can be represented both through program parameters, and by means of a file containing their names. 
In both cases the first file in a list – *primary* – is compared to the others – *secondary* – by turns.
The *primary* file specifies a set of compared chromosomes (unless it is limited by the option ```-c|--chr``` or ```-g|--gen```). 
If one or more chromosomes are represented in the *primary* file, only they will be compared.<br>
Also, in the case of *wiggle*, the *primary* defines a resolution (see option ```-s|--space``` for more information).

Be careful by using standard naming conventions like *abc?.bed*, *a\*.bed*. 
You must be sure that first of the lexically recognized files is really primary, and that all other files really need to be compared.<br>
For more details about list of file names as a file see option ```-l|--list```.

### Options description
Enumerable option values are case insensitive.

```-a|--align```<br>
Indicate that input bed files are *alignments*, so density correlation would be performed.<br>
Since *alignment* may be considered as an 'instance' of *ordinary bed*, in theory both cases should give the same result. 
It is true while the *alignments* as well as the *ordinary beds* have no *ambiguous* reads/features (see ```–i|--info``` option for definition). 
But in practice this condition is almost never achieved, and these ambiguities are resolved by different way. 
Consequently the results can be dramatically different.<br>
In addition, *ordinary beds* are compared using a separate, ultra-fast algorithm.<br>
If alignments are treated without this option, bioCC will print a warning message.<br>
If *ordinary bed* files have features of different sizes and are treated with this option, 
**bioCC** will print a cancel message and complete.
.

```-g|--gen <name>```<br>
Chromosome sizes file, reference genome library, or single nucleotide sequence.<br>
Genome library is a directory contained nucleotide sequences for each chromosome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.<br>
If ```name``` is a *.fa[.gz]* file, **bioCC** accepts the corresponding chromosome as the only treated.<br>
Otherwise first the program searches for *.fa* files in the directory ```name```. 
If there are no such files in this directory, **bioCC** searches for *.fa.gz* files.<br>
If chromosome is specified by option ```–c|--chr```, the program searches for the corresponding *.fa[.gz]* file.<br>
The chromosome sizes file is recognized by *.sizes* extension.<br>
The difference between chromosome sizes file and genome library/file is that in the latter all the undefined regions in the reference genome (gaps), will be excluded from calculation.<br>
Undefined regions are regions with only ambiguous reference characters ‘N’ in them.<br>
It can make sense in some special cases. Synchronous gaps in the alignments slightly increase the degree of orderliness. 
For instance, Pearson correlation coefficient for two sequences, recalled by different aligners (MOSAIK and SMALT in this case) from one model control sequence (‘input’), 
is equal to 0.86 when ignoring the gaps, and 0.70 when they are taken into account (with the default settings of resolution and gap length).<br>
The minimal length of accounting gaps is managed by ```--gap-len``` option.<br>
For example, chromosome 1 from mm9 library contains 14 regions, separated by gaps with length more then 400 bps, and 10 regions, separated by gaps with length more then 1000.<br>
Indicating genome library has the same effect as ```-f|--fbed``` option, where 'template' is a set of defined regions.<br>
The program searches for gaps once. On subsequent calls with the same length of gaps, it uses the search results stored in the specified directory.<br>
The program also generates once a chromosome sizes file in the same directory.<br>
One can obtain a genome library in  UCSC: ftp://hgdownload.soe.ucsc.edu/goldenPath/ or in Ensemble: ftp://ftp.ensembl.org/pub/release-73/fasta storage. 
In the second case please copy genomic sequences with the same masked type only, f.e. unmasked (‘dna'), since program does not recognise mask’s types.<br>
This option is required.

```--gap-len <int>```<br>
Minimal length of undefined nucleotides region which is taken as a gap. 
For more details see ```-g|--gen``` option.<br>
Ignored for chromosome sizes file, specified by ```-g|--gen```.<br>
Default: 1000

```-d|--dupl <OFF|ON>```<br>
Accept or deny duplicated reads to participate in density calculation.<br>
This option only makes sense for the *alignments*.<br>
Default: ```ON```

```-l|--list <file>```<br>
External list of compared files: a plain text file, in which each file name is located on a separate line.<br>
Lines beginning with ‘#’ are considered as comments and will be skipped.<br>
Empty lines are allowed.<br>
This option abolishes input files as parameters.

```-c|--chr <chars>```<br>
Treat specified chromosome only. Samples of option’s value: 1, 20, X.<br>
Reduces run time on 1.5-20 times depending on how far this chromosome is placed in an input *alignment* or *wig*.<br>
For *ordinary beds* it has no time-improvement effect: any result appears quickly.<br>

```-r|--cc <P,S>```<br>
Correlation coefficient, in any combination: ```P``` – Pearson, ```S``` – signal.<br>
For more details see [Pearson and signal correlation](#pearson-and-signal-correlation).<br>
Default: Pearson

```-s|--space <int>```<br>
Resolution: the length of windows (span) in bp by which reads will be counted to define a density.<br>
Read is considered belonging to span if it`s centre is placed within the span. 
Thus each read is counted once. 
Then, when using the reference genome from the input sequences, the reference undefined regions (gaps) are preliminarily removed from the compared sequences: 
the regions separated by a gap are merged.<br>
As a result, the program compares the actual read density distributions.<br>
This option is topical for the *alignments* only.<br>
*Wiggles* are already presented according to their resolution.
**bioCC** takes the resolution of *primary wiggle* as the basic for treatment. 
If *secondary wiggle* resolution is differ, it will be reduced to the base. 
For example, if the primary resolution is 10 and the secondary resolution is 1, then every 10 partitions in the secondary data will be merged so finally the secondary resolution will also be 10. 
Contrariwise, if secondary resolution is 100, then each partition of the secondary data will be divided into 10 with the same density.<br>
The best way is to use *wiggles* with the same resolutions.<br>
Default: 100.

```-p|--pr-cc <IND,TOT>```<br>
Print coefficients, in any combination: ```IND``` - for each chromosome individually, ```TOT``` - total.<br>
If only one chromosome is specified, the total coefficient is not output as an identical.<br>
Default: ```IND```.

```-f|--fbed <file>```<br>
'Template' *ordinary* bed file with features that defines compared regions.<br>
In practice the most interesting case is comparing detected peaks across the distributions, for what this option is mainly constructed.<br>
If options ```-b|--bin-width``` and ```--f``` are not defined, only total coefficients for the all regions are printed.<br>
If this option is established, predefined reference regions are ignored (assuming that the features found by the peak callers are always setting in the defined regions).<br>
This option is ignored for the *ordinary bed* files.

```-e|--ext-len <int>```<br>
Value by which all features in a 'template' *bed* file or in a *primary ordinary bed* file should be stretched in both directions before comparison.<br>
If set, all the features from 'template' or *primary bed* file will be stretched before processing: *start* positions are decreased for this value, *end* positions are increased.<br>
If stretched features become intersected, they are joined.<br>
This option is constructed for the special cases, for instance, 
ChIP-seq treatment of transcription factor binding sites (TFBS). 
Such binding sites have a well-defined uniform length (8-20 bp), while the length of corresponded enriched regions can reach 500-2000 bp. 
A direct comparison of the 'real' *bed* file and *bed* file after peak caller can give a low result even in case of complete correspondence of 'real' and recovered peaks.<br>
Stretching initial TFBS lengths allow correlate these data correctly.<br>
Default: 0

```--ext-step <int>```<br>
If set, activates the mode of consecutive calculation of the coefficients for stretching features in *primary ordinary* bed file with the stated step. 
The maximum value of the extension is limited by the value of ```--e|--ext-len``` option.<br>
This option is topical for *ordinary bed* files only.<br>
Default: 0 (no step calculation)

```-b|--bin-width <float>```<br>
If set, forces to consolidate coefficients into bins, and print histogram values.<br>
Histogram shows count of coefficients (frequency) within some value ranges (bins). 
It is printed as a list of pairs *\<bin upper bound\>\<count in bin\>*. 
Negative coefficients have been turning to absolute during consolidation.<br>
This option defines the width of bin as a part of 1, so it takes a value in the range 0-1.<br>
Empty bins at the edges are not printed.<br>
For example, if –b value is 0.1, and all coefficients are placing in the range 0.65 to 0.85, only three bins [0.6-] 0.7, [0.7-] 0.8, [0.8-] 0.9 would be printed.<br>
This option is topical only with option ```-f|--fbed``` or for reference genome library.<br>
Default: 0 (no consolidation).

```--sort <RGN|CC>```<br>
If set, force to output coefficients calculated for each region as a list of pairs *\<number-of-region\>\<coefficient\>*. 
```RGN``` value prescribes list to be sorted by region’s number, ```CC``` – by coefficient.<br>
If both of the coefficients are declared, list is sorted by Pearson.<br>
First region number is 1.<br>
This option is topical only with option ```-f|--fbed``` or for reference genome library.<br>

```--rn <OFF|ON>```<br>
Normalize regions before calculation. <br>
Normalization means levelling distribution’s patterns by their maximal values across regions.<br>
In spite the fact that of each pair of regions will be always compared correctly, the common result for all the regions may been unfairly falling down. It occurs when patterns levels are appreciably differ from region to region.<br>
This option is assigned to eliminate such effect.<br>
This option is topical only with option ```-f|--fbed``` or for reference genome library.<br>
Default: ```ON```

```-i|--info <LAC|NM|CNT|STAT>```<br>
Output information about number of items (features/reads/intervals).<br>
```LAC```: laconic output. This value minimizes the output as possible to remain clear. It is constructed mainly for using in batch file. 
But it does not disable warning output, stated by option ```-w|--warn```.<br>
```NM```:  brief output. Prints file names without any additional information.<br>
```CNT```:  prints file names and number of all and accepted items, if they are different.<br>
```STAT```: prints item ambiguities statistics, if they exist.<br>
When we are speaking about *bed* files, there are some issues which can be normal or may be not – that depends of file’s destination. 
For instance, duplicated and crossed features are normal for the *alignments*, but rather unusual for the *ordinary beds*. 
On the other hand, features with different length are typical for the *ordinary beds*, but they are a rare case for the *alignments*. 
All these situations we call a conditional term *ambiguities*.<br>
It concerns to 'template' *bed* as well.<br>
**bioCC** treats ambiguities the most natural way. 
For *ordinary bed* it merges crossed or adjacent features, and ignores submerging and duplicated features. 
For *alignment* it accepts duplicated reads by default (in fact, all the ambiguities are permitted for the alignment, except different length of reads).<br>
Thus, not all records present in the file can be accepted.<br>
In some circumstances you need to be aware of these issues. 
There are two methods used to identify them: detailed and summary.<br>
The ```STAT``` value provides the summary method. 
It forces to display number of all recognized certain type ambiguities, and appropriate treatment.<br>
In particular, this is a simple way to know the number of duplicated reads.<br>
The detailed method is managed by option ```-w|--warn```.<br>
Note, that in contrast to *bed* with their features/reads, the number of intervals does not match the number of lines in the *wiggle* file. 
All contiguous bases with the same data value are considered to be a single interval, although it may correspond to several file lines.<br>
Default: ```NM```

```-w|--warn```<br>
Output ambiguity for each feature, if they exist.<br>
If it is specified, information about type of ambiguity, number of line where it occurs, and resulting treatment would be printed each time it happens.<br>
Duplicated reads in *alignment* are not printed when using this method as they are considered normal for the alignment, but they are reported in the summary method, see ```-i|--info``` option.

```-o|--out```<br>
Duplicate standard output to **bioCC_out.txt** file.<br>
It is analogue of **tee** Linux command and is rather useful by calling **bioCC** under Windows.


## Pearson and signal correlation
For the data under consideration, the Pearson comparison is correct, since the data are linear in nature.<br>
While we consider coverage/density as the data distributed along regular dimension (scaled chromosome’s length), signal's method is appropriate as well.<br>
The only difference between them is the subtraction of the mean value in covariance on Pearson’s coefficient.

What consequences does it entail, and which coefficient is better to choose?

A. Coverage/density distributions.<br>
To be more clear there are 3 illustrations of pair of functions which are quite near to the fragments of real coverage.<br>
Both of coefficients demonstrate value’s normalization independence (![fig 1,3](https://github.com/fnaumenko/bioCC/tree/master/pict/Signal-Pearson.png)).<br>
signal method is a bit more sensible. But if we are interesting of measure of the degree of distribution patterns (shapes) only, this method becomes inappropriate when some distribution have non-zero background (*DC offset* in terms of signal function), as it shows in ![fig 2](https://github.com/fnaumenko/bioCC/tree/master/pict/Signal-Pearson.png)<br>
Briefly, Pearson’s method is mainly recommended for the distributions comparison.<br>
But if background’s level is considered part of the measure of similarity, in this case signal method would be preferable (it is true if resolution is sufficiently low to average background’s irregularity).

B.  Bed features.
Bed features can be accounted as discrete function accepted one of two values: 0 or 1. In that case signal method becomes inappropriate due to obvious reason: it counts intersections only. 
For example, by comparison two mutually complementary beds (while function1 have zero value every time when function2 have non-zero and vice versa), signal coefficient would be 0. 
Although the correct answer is -1.<br>
Thus, in the case of features, only the Pearson method is correct.

##
If you face to bugs, incorrect English, or have commentary/suggestions, please do not hesitate to write me on fedor.naumenko@gmail.com
