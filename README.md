### What is CAMMiQ?

CAMMiQ is a software tool for microbial identification and quantification. Specifically, it builds a compact index to manage a database of microbial genomes (in particular, bacterial/viral genomes in NCBI [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) database) and answers queries of the following form: given a set of high throughput sequencing (HTS) reads obtained from a mixture of genomes in the database, identify the corresponding genomes (in the database) and compute their relative abundance.

**_To cite CAMMiQ, please mention_ https://www.biorxiv.org/content/10.1101/2020.06.12.149245v2** 

### How to install CAMMiQ?
**Dependencies**: Though several scripts (to compile the code and download the database) are written using the Bash shell and python, our core programs to build the index and query the set of reads are written in C++11, and need to be compiled using a relatively recent version (The current build of CAMMiQ was tested with version 5.2.0) of gcc that will support C++11. Multithreading is handled using pthread and [OpenMP](https://en.wikipedia.org/wiki/OpenMP).

In addition, you'll need the following components to compile the sources:
* https://github.com/jlabeit/parallel-divsufsort which constructs suffix arrays in a parallelized and lightweight fashion.
* https://github.com/martinus/robin-hood-hashing which is a faster and more memory efficient alternative of STL unordered_map.
* [IBM ILOG CPLEX Optimization Studio](https://www.ibm.com/products/ilog-cplex-optimization-studio) - CAMMiQ requires its c++ interface for solving a (mixed) integer linear program (ILP) to figure out the most likely composition of genomes in the mixture based on the distribution of querying reads. 

To install CAMMiQ, just clone our repository and run the script install_CAMMiQ.sh we have prepared - it will automatically download the related repos above and compile their codes. The only exception is, you'll have to first download and install the latest version (e.g., 12.9.0) of IBM ILOG CPLEX Optimization Studio by yourself.  
```
git clone https://github.com/algo-cancer/CAMMiQ
./install_CAMMiQ.sh --cplex-dir <CPLEX_DIR>
```
where ```<CPLEX_DIR>``` should be replaced with the directory of your CPLEX_Studio.

### How do I use CAMMiQ?
To begin using CAMMiQ, you will first need to index the input genomes. 

#### What is CAMMiQ index composed of?
CAMMiQ index is composed of three parts - all of them are necessary to query the input reads: 
* (i) A set of representative (i.e., maximally sparsified) *shortest unique substrings* on each input genome, organized as a prefix tree (trie) and hashed based on their common prefixes, encoded in binary file ```index_u.bin1```.
* (ii) A set of representative (i.e., maximally sparsified) *shortest doubly unique substrings* on each input genome, organized as a prefix tree (trie) and hashed based on their common prefixes, encoded in binary file ```index_d.bin2```.
* (iii) Other meta-information of the input genomes, including ```genome_lengths.out```, a text file containing the genome lengths; ```unique_lmer_count_u.out```, a text file containing the number of unique *L*-mers on each genome; and ```unique_lmer_count_d.out```, a text file containing the number of doubly unique *L*-mers on each genome.

#### How do I construct the index?
You'll need to run ```./cammiq --build [options] [parameters]``` from command line, where ```[options]``` include 
  * ```--unique``` CAMMiQ will build the set of *unique substrings* from each input genome, consisting of ```index_u.bin1```, ```genome_lengths.out``` and ```unique_lmer_count_u.out```. Note that ```--unique``` is the **default** option if nothing is specified in the command line instructions.
  * ```--doubly_unique``` CAMMiQ will build the set of *doubly unique substrings* from each input genome, consisting of ```index_d.bin2```, ```genome_lengths.out``` and ```unique_lmer_count_d.out```.
  * ```--both``` CAMMiQ will build its complete indices consisting of all above files.

On the other hand, ```[parameters]``` include the following list of (possibly mandatory) parameters.
* ```-f <MAP_FILE>``` **Mandatory**. ```<MAP_FILE>``` gives a list of reference genomes, e.g., all/selected complete genomes in RefSeq for the bacterial, archaeal, and viral domains (downloaded with ```CAMMiQ-download```), which constitute CAMMiQ's database, possibly alongwith NCBI's taxonomic information. The input lines in ```<MAP_FILE>``` should contain at least 4 tab-delimited fields; from left to right, they are: 
  * File names
  * Genome IDs (encoded in the index files)
  * NCBI taxonomic IDs
  * Organism names
  
  Here is an example format of ```<MAP_FILE>```:
  ``` 
  GCF_000010525.1_ASM1052v1_genomic.fna	1	7	Azorhizobium caulinodans ORS 571
  GCF_000007365.1_ASM736v1_genomic.fna	2	9	Buchnera aphidicola str. Sg (Schizaphis graminum)
  GCF_000218545.1_ASM21854v1_genomic.fna	3	11	Cellulomonas gilvus ATCC 13127
  GCF_000020965.1_ASM2096v1_genomic.fna	4	14	Dictyoglomus thermophilum H-6-12
  GCF_000012885.1_ASM1288v1_genomic.fna	5	19	Pelobacter carbinolicus DSM 2380
  ......
  ```
  
  (As a shortcut, ```-f``` can alternatively take a list of fasta files to build an index on these files. However, to query the index you will again need to organize the information of these fasta files in a ```<MAP_FILE>``` and use it as the input in the ```--query``` mode.
* ```-d <FASTA_DIR>``` **Mandatory**. ```<FASTA_DIR>``` should contain the list of (fasta) file names given in ```<MAP_FILE>```.
* ```-k <int>``` **Optional**. The minimum length of a unique or doubly-unique substring to be considered in CAMMiQ index. Default value is ```k = 26```.
* ```-L <int>``` **Optional (but strongly recommended)**. Potential read length in a query supported by CAMMiQ index. Default value is ```L = 100```, which fits best for reads with length ```100```; if, for instance, the reads in your query have length ```75```, then you are expected to build an index by specifying ```L = 75```. 
* ```-Lmax <int>``` **Optional (but recommended)**. The maximum length of a unique or doubly-unique substring to be considered in CAMMiQ index. Default value is ```Lmax = 50``` or ```Lmax = 0.5 * L```.
* ```-h <int>|<int1 int2>``` **Optional**. Length of the common prefixes of the unique or doubly-unique substrings to be hashed. Default value is ```h = 26```.
  * Note a: The value of ```h``` is required to be *less than or equal to* ```k```. 
  * Note b: ```-h``` option can take in one or two integer values. With one input value ```h0```, CAMMiQ will set both hash lengths (for the collection of unique substrings and for the collection of doubly-unique substrings) ```h0```; with two input values ```h1``` and ```h2```, CAMMiQ will set the hash length for unique substrings ```h1``` and the hash length for doubly-unique substrings ```h2```.
* ```-t <int>``` **Optional**. Number of threads used during CAMMiQ's index construction. Note that CAMMiQ uses OpenMP during its index construction, which, by default, is 'auto-threaded' (i.e. attempts to use all available CPUs on a computer).

#### How do I query the collection of (metagenomic) reads?
You'll need to run ```./cammiq --query [options]``` from command line, where ```[options]``` specifies the following list of (possibly mandatory) parameters.
* ```-f <MAP_FILE>``` **Mandatory**.
* ```-q (-Q) <QUERY_FILE>``` **Mandatory**.
* ```-I <INDEX_DIR>``` **Mandatory**.
* ```-o <OUTPUT_FILE>``` **Optional (but strongly recommended)**.
Here is an example of ```<OUTPUT_FILE>```
```
TAXID	ABUNDANCE	NAME
1795	0.052824	Mycolicibacterium neoaurum VKM Ac-1815D
2105	0.045993	Mycoplasma leachii PG50
547143	0.051095	Hydrogenobaculum sp. 3684
547145	0.048998	Hydrogenobaculum sp. SHO
547146	0.042636	Hydrogenobaculum sp. SN
1476577	0.044080	Candidatus Saccharibacteria oral taxon TM7x
......
```

* ```-e <int>``` **Optional (but strongly recommended)**.
* ```-h <int>|<int1 int2>``` **Optional (but need special attention)**. 
* ```-t <int>``` **Optional**.


