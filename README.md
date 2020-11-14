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
where ```<CPLEX_DIR>``` should be replaced with the directory of your **CPLEX_Studio**.

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
* ```-f <MAP_FILE>``` **Mandatory**. ```<MAP_FILE>``` gives a list of reference genomes in fasta format, e.g., all/selected complete genomes in RefSeq for the bacterial, archaeal, and viral domains (downloaded with ```CAMMiQ-download```), which constitute CAMMiQ's database, possibly alongwith NCBI's taxonomic information. The input lines in ```<MAP_FILE>``` should contain at least 4 tab-delimited fields; from left to right, they are: 
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
  
  (As a shortcut, ```-f``` can alternatively take a list of fasta files to build an index on these files. However, to query the index you will again need to organize the information of these fasta files in a ```<MAP_FILE>``` and use it as the input in the ```--query``` mode.)

  *Note.* It is in fact **recommended** that genome IDs grow from 1 to the number of genomes to be indexed. This feature was designed to minimize the memory usage when constructing CAMMiQ's index. See "[What is the expected computational cost of CAMMiQ?](###-What-is-the-expected-computational-cost-of-CAMMiQ?)" for more details.

  *A final note about ```<MAP_FILE>```.* In general, genome IDs in two different lines are *expected* to be different. However, sometimes genome IDs in two distinct lines *can be the same*; in other words, two fasta files can have *the same* genome ID, which means that CAMMiQ will treat them as contigs from the same genome. This functionality turns out to be useful **in some special cases**, e.g., in the ```--read_cnts``` style query, you may want to count the number of reads originated from each genus, yet the corresponding information of each fasta file is given at species/strain level - in this case you may want to use the same genome ID for the genomes belonging to the same species/strain to build the index. You need to be extremely careful about the query scenario when you have two identical genome IDs for different genomes (fasta files) just to avoid any **misuse** of CAMMiQ. 

* ```-D <FASTA_DIR>``` **Mandatory**. ```<FASTA_DIR>``` should contain the list of (fasta) file names given in ```<MAP_FILE>```.
* ```-k <int>``` **Optional**. The minimum length of a unique or doubly-unique substring to be considered in CAMMiQ index. Default value is ```k = 26```.
* ```-L <int>``` **Optional (but strongly recommended)**. Potential read length in a query supported by CAMMiQ index. Default value is ```L = 100```, which fits best for reads with length ```100```; if, for instance, the reads in your query have length ```75```, then you are expected to build an index by specifying ```L = 75```. 
* ```-Lmax <int>``` **Optional (but recommended)**. The maximum length of a unique or doubly-unique substring to be considered in CAMMiQ index. Default value is ```Lmax = 50``` or ```Lmax = 0.5 * L```.
* ```-h <int>|<int1 int2>``` **Optional**. Length of the common prefixes of the unique or doubly-unique substrings to be hashed. Default value is ```h = k``` or ```h = 26```.
  * Note a: The value of ```h``` is required to be *less than or equal to* ```k```. 
  * Note b: ```-h``` parameter can take in one or two integer values. With one input value ```h0```, CAMMiQ will set both hash lengths (for the collection of unique substrings and for the collection of doubly-unique substrings) ```h0```; with two input values ```h1``` and ```h2```, CAMMiQ will set the hash length for unique substrings ```h1``` and the hash length for doubly-unique substrings ```h2```.
* ```-t <int>``` **Optional**. Number of threads used during CAMMiQ's index construction. Note that CAMMiQ uses OpenMP during its index construction, which, by default, is 'auto-threaded' (i.e. attempts to use all available CPUs on a computer).

Example usage:
```
./cammiq --build --doubly_unique -k 26 -L 100 -Lmax 50 -f genome_map1.txt -D /data/fasta_dir/ -t 32 
/* h should not exceed k */
./cammiq --build --doubly_unique -k 26 -L 100 -Lmax 50 -h 25 -f genome_map2.out -D /data/fasta_dir/ -t 64
./cammiq --build --both -k 21 -L 75 -Lmax 75 -h 21 -f bacteria1.fa bacteria2.fa bacteria3.fa
```

#### How do I query the collection of (metagenomic) reads?
Similarly, you'll need to run ```./cammiq --query [options] [parameters]``` from command line, where ```[options]``` include
  * ```--read_cnts``` **Optional**. If ```--read_cnts``` is specified, then CAMMiQ will not produce its standard output (see below ```-o``` option). Instead, CAMMiQ outputs a non-negative matrix where each row corresponds to a query (fastq file); each column corresponds to an **NCBI taxonomic ID** (attention: not a genome ID!); each entry gives the number of reads in a given query that CAMMiQ assigned uniquely to the corresponding taxon.
  
  Here is an example output when CAMMiQ finds ```--read_cnts``` in a command line:
  ```
    7 9 11 14 19
  query_1.fq 0 0 15 100 10000
  query_2.fq 20 0 125 1800 0
  query_3.fq 0 0 0 0 0
  ......
  ```
  
  Note that the output file name can be specified with ```-o``` parameter. 
  
  * ```-doubly_unique``` **Optional**. Only valid when ```--read_cnts``` is specified. CAMMiQ will resolve the ambiguous read counts brought by doubly-unique substrings, and assign each of those reads that only contain doubly-unique substrings from two distinct taxa to one specific taxon. 

and ```[parameters]``` include the following list of (possibly mandatory) parameters.
* ```-f <MAP_FILE>``` **Mandatory**. You should use the same ```<MAP_FILE>``` when building the index in your queries. **Attention:** CAMMiQ is not in charge of verifying the format or correctness (meaning that you use exactly the same file for building the index and querying) of a ```<MAP_FILE>```. When your input ```<MAP_FILE>``` for querying is different from what you used for building the index, some potential "*undefined behavior*" could happen when running CAMMiQ.
* ```-q (-Q) <QUERY_FILE(S)>``` **Mandatory**. ```<QUERY_FILE(S)>``` can be either a list of fastq files, or a directory containing the list of fastq files in your query. A capitalized ```-Q``` indicate the input is a directory, that is,  
  * ```-q``` CAMMiQ takes the list of fastq files.
  * ```-Q``` CAMMiQ takes a directory which contains the list of fastq files.
* ```-i <INDEX_FILES>``` **Mandatory**. As discussed in [What is CAMMiQ index composed of?](####-What-is-CAMMiQ-index-composed-of?), ```<INDEX_FILES>``` include ```*.bin1```, ```*.bin2```, ```genome_lengths.out```,  ```unique_lmer_count_u.out```, ```unique_lmer_count_d.out```. ```-i``` parameter **should take** at least ```*.bin1``` and ```*.bin2```; the default location for CAMMiQ to find those meta-information (```genome_lengths.out```,  ```unique_lmer_count_u.out```, ```unique_lmer_count_d.out```) files is ```./```; if they are not stored in the current directory ```./```, you'll need to specify them explicity.

  One exception is when ```--read_cnts``` is specified. In this case you don't have to input ```unique_lmer_count_u.out``` and ```unique_lmer_count_d.out```; even if you give them, CAMMiQ will *ignore* these files.

* ```-o <OUTPUT_FILE>``` **Optional (but strongly recommended)**. CAMMiQ's standard output file, with ```<OUTPUT_FILE>``` specifying the file name. You may want a different
  
  Here is an example format of ```<OUTPUT_FILE>```
  ```
  sample.fastq
  TAXID	ABUNDANCE	NAME
  1795	0.052824	Mycolicibacterium neoaurum VKM Ac-1815D
  2105	0.045993	Mycoplasma leachii PG50
  547143	0.051095	Hydrogenobaculum sp. 3684
  547145	0.048998	Hydrogenobaculum sp. SHO
  547146	0.042636	Hydrogenobaculum sp. SN
  1476577	0.044080	Candidatus Saccharibacteria oral taxon TM7x
  ......
  ```

* ```-e <int>``` **Optional (but strongly recommended)**. Estimated sequencing error (substitution) rate. Default value is ```e = 0.0``` (i.e., no sequencing errors); for typical Illumina reads, you can try something like ```-e 0.008``` or ```-e 0.01```; however, the more accurate your estimation of error rate, the more accurate CAMMiQ   
* ```-h <int>|<int1 int2>``` **Optional (but need special attention)**. Same as what is described above in the ```--build``` option, ```-h``` can take in one or two integer values. You don't need to speicfy ```-h``` since it's encoded in the corredponding index (```*.bin1``` and ```*.bin2```) files - but if something is found here and not matching the value encoded in the index files, CAMMiQ will report an error.
* ```-t <int>``` **Optional**. Number of threads used during CAMMiQ queries. Same as what is described above in the ```--build``` option. Note that CAMMiQ query has two phases: for the first phase of assigning reads to one or two genomes, the default number of threads is 1; however, the second phase will use IBM ILOG CPLEX Optimization Studio to produce its metagenomic quantification results (or read count results if ```--read_cnts``` is specified), which, by default, is again 'auto-threaded' (i.e. attempts to use all available CPUs on a computer). If your Linux mahcine has more than 32 CPUs, then by default 32 threads will be used.

#### How do I query single cell RNA-seq reads?

### What is the expected computational cost of CAMMiQ?
