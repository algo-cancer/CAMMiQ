### What is CAMMiQ?

CAMMiQ is a software tool for microbial identification and quantification. Specifically, it builds a compact index to manage a database of microbial genomes (in particular, bacterial/viral genomes in NCBI [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) database) and answers queries of the following form: given a set of high throughput sequencing (HTS) reads obtained from a mixture of genomes in the database, identify the corresponding genomes (in the database) and compute their relative abundance.

**_To cite CAMMiQ, please mention_ https://www.biorxiv.org/content/10.1101/2020.06.12.149245v2** 

### How to install CAMMiQ?
**Dependencies**: Though several scripts (to compile the code and download the database) are written using the Bash shell and python, our core programs to build the index and query the set of reads are written in C++11, and need to be compiled using a relatively recent version (The current build of CAMMiQ was tested with version 5.2.0) of gcc that will support C++11. Multithreading is handled using pthread and [OpenMP](https://en.wikipedia.org/wiki/OpenMP).

In addition, you'll need the following components to compile the sources:
* https://github.com/jlabeit/parallel-divsufsort which constructs suffix arrays in a parallelized and lightweight fashion.
* https://github.com/martinus/robin-hood-hashing which is a faster and more memory efficient alternative of STL unordered_map.
* [IBM ILOG CPLEX Optimization Studio](https://www.ibm.com/products/ilog-cplex-optimization-studio) or [Gurobi Optimizer](https://www.gurobi.com/products/gurobi-optimizer/) - CAMMiQ requires either c++ interface for solving a (mixed) integer linear program (ILP) to figure out the most likely composition of genomes in the mixture based on the distribution of querying reads. 

To install CAMMiQ, just clone our repository and run the script install_CAMMiQ.sh we have prepared - it will automatically download the related repos above and compile their codes. The only exception is, you'll have to first download and install the latest version of IBM ILOG CPLEX Optimization Studio (e.g., 12.9.0) or Gurobi Optimizer (e.g., 9.0.0) by yourself.  
```
git clone https://github.com/algo-cancer/CAMMiQ
./install_CAMMiQ.sh --cplex-dir <CPLEX_DIR>
```
or
```
./install_CAMMiQ.sh --gurobi-dir <GUROBI_DIR> --gurobi-version <GUROBI_VERSION>
```
where ```<CPLEX_DIR>``` should be replaced with the **absolute** path to your **CPLEX_Studio** if you use CPLEX; if you use Gurobi, ```<GUROBI_DIR>``` should be replaced with the **absolute** path to **linux64** (under the Gurobi directory) and additionally you'll need to tell the compiler the version of your Gurobi Optimizer (usually both ```a.b``` and ```a.b.c``` formats should work).

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
  * The index file names can be modified through ```-i``` option (see below).

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
  
  As real examples of map files, the 4 index datasets presented in the paper can be found in https://drive.google.com/drive/u/0/folders/1iL0pZ3Jb3QS6AOwO_IGLS2ZDFQJEHpn0.

  *Note.* It is in fact **recommended** that genome IDs grow from 1 to the number of genomes to be indexed. This feature was designed to minimize the memory usage when constructing CAMMiQ's index. See "[What is the expected computational cost of CAMMiQ?](#-What-is-the-expected-computational-cost-of-CAMMiQ?)" for more details.

  *A final note about ```<MAP_FILE>```.* In general, genome IDs in two different lines are *expected* to be different. However, sometimes genome IDs in two distinct lines *can be the same*; in other words, two fasta files can have *the same* genome ID, which means that CAMMiQ will treat them as contigs from the same genome. This functionality turns out to be useful **in some special cases**, e.g., in the ```--read_cnts``` style query, you may want to count the number of reads originated from each genus, yet the corresponding information of each fasta file is given at species/strain level - in this case you may want to use the same genome ID for the genomes belonging to the same species/strain to build the index. You need to be extremely careful about the query scenario when you have two identical genome IDs for different genomes (fasta files) just to avoid any **misuse** of CAMMiQ. 

* ```-D <FASTA_DIR>``` **Mandatory**. ```<FASTA_DIR>``` should contain the list of (fasta) file names given in ```<MAP_FILE>```.
* ```-k <int>``` **Optional**. The minimum length of a unique or doubly-unique substring to be considered in CAMMiQ index. Default value is ```k = 26```.
* ```-L <int>``` **Optional (but strongly recommended)**. Potential read length in a query supported by CAMMiQ index. Default value is ```L = 100```, which fits best for reads with length ```100```; if, for instance, the reads in your query have length ```75```, then you are expected to build an index by specifying ```L = 75```. 
* ```-Lmax <int>``` **Optional (but recommended)**. The maximum length of a unique or doubly-unique substring to be considered in CAMMiQ index. Default value is ```Lmax = 50``` or ```Lmax = 0.5 * L```.
* ```-h <int>|<int1 int2>``` **Optional**. Length of the common prefixes of the unique or doubly-unique substrings to be hashed. Default value is ```h = k``` or ```h = 26```.
  * Note a: The value of ```h``` is required to be *less than or equal to* ```k```. 
  * Note b: ```-h``` parameter can take in one or two integer values. With one input value ```h0```, CAMMiQ will set both hash lengths (for the collection of unique substrings and for the collection of doubly-unique substrings) ```h0```; with two input values ```h1``` and ```h2```, CAMMiQ will set the hash length for unique substrings ```h1``` and the hash length for doubly-unique substrings ```h2```.
* ```-i <INDEX_FILES>``` **Optional**. ```<INDEX_FILES>``` should include one or two file names if ```-i``` is found in ```--build``` option. With ```--unique``` option, the ```*.bin1``` index will be set to a single file name specified here; with ```--doubly_unique``` option the ```*.bin2``` index will be set to a single file name specified here; with ```--both``` option, the name of ```*.bin1``` and ```*.bin2``` file will be reset respectively (See [What is CAMMiQ index composed of?](#-What-is-CAMMiQ-index-composed-of?)). The default file names (if ```-i``` is not specified) are ```index_u.bin1``` and ```index_d.bin2```. The meta-information file names (i.e., ```genome_lengths.out```, ```unique_lmer_count_u.out``` and ```unique_lmer_count_d.out```) cannot be changed; the default location (directory) to store these files is ```./```; however, if the directory of index file names is not ```./```, then the meta-information files will be written to the same directory as that is given in the index file names.
* ```-t <int>``` **Optional**. Number of threads used during CAMMiQ's index construction. Note that CAMMiQ uses OpenMP during its index construction, which, by default, is 'auto-threaded' (i.e. attempts to use all available CPUs on a computer).

**Example usage**:
```
cammiq --build --doubly_unique -k 26 -L 100 -Lmax 50 -f genome_map1.txt -D /data/fasta_dir/ -t 32 
```
The above command line instruction builds doubly-unique substrings (encoded in ```index_d.bin2```) from the genomes listed in ```genome_map1.txt``` stored under the directory ```/data/fasta_dir/```; with substring lengths ranging from [26, 50]; with 32 threads; and to support query reads of length (roughly) 100. Note that -h was not specified so CAMMiQ takes its default value 26.
```
cammiq --build -k 26 -L 120 -Lmax 75 -h 25 -f genome_map2.out -D /data/fasta_dir/ -i /data/Indices/index_u_test.bin1 -t 64
cammiq --build --unique -k 26 -L 120 -Lmax 75 -h 25 -f genome_map2.out -D /data/fasta_dir/ -i /data/Indices/index_u_test.bin1 -t 64
```
Either of the above two command line instructions will build unique substrings (encoded in ```/data/Indices/index_u_test.bin1``` and meta-information files written to ```/data/Indices/```) from the genomes listed in ```genome_map2.txt``` stored again under the directory ```/data/fasta_dir/```; with substring lengths ranging from [26, 75]; with 64 threads; and to support query reads of length (roughly) 120. Note that -h was set to 25, which does not exceed the maximum allowed value ```k=26```.
```
cammiq --build --both -k 21 -L 75 -Lmax 75 -h 21 -f bacteria1.fa bacteria2.fa bacteria3.fa
```
The above command line instruction builds both unique and doubly substrings from the 3 genomes ```bacteria1.fa```, ```bacteria2.fa``` and ```bacteria3.fa```; with substring lengths ranging from [21, 75]; to support query reads of length (roughly) 75. Usually there are much more reference genomes; as a result, listing them all as ```-f``` parameters is not advocated, though supported by CAMMiQ. And since ```-t``` was not specified, CAMMiQ will attempt to use the maximum number of threads according to openmp. 

#### Is there a default database supported by CAMMiQ?
Yes. Similar to many other products, CAMMiQ supports selected complete genomes from **the latest version of** NCBI's RefSeq database, which can be easily downloaded through the python script ```CAMMiQ-download```. The command line options of ```CAMMiQ-download``` include

* ```--dir <DATABASE_PATH>``` The directory to maintain the downloaded genomes. Default path is the current directory ```./```.
* ```--taxa bacteria|viral|archaea|all``` The selected collection of complete genomes from RefSeq to be downloaded. Currently there are four options supported by CAMMiQ: users can choose to download the list (determined by ```--sample``` option) of complete bacterial, viral, archaeal genomes or all three categories from the *latest* version of RefSeq database. 
* ```--sample none|taxid|species``` The option ```none``` keeps all genomes; ```taxid``` keeps one representative genome for each taxonomic ID; ```species``` keeps one representative genome for each species level taxonomic ID (these taxonomic IDs are based on the latest assembly summary files maintained by RefSeq). Note that if ```taxid``` or ```species``` is specified, genomes marked as "representative genome" or "reference genome" (again in the assembly summary files) have the priority to be kept; if no such genomes are found with in each (species level) taxID, then a random genome from this taxID will be kept.
* ```--unzip``` This option decompresses the ```*.gz``` files from RefSeq into ```*.fna``` files.
* ```--quiet``` Disable the output produced by wget to trace the download progress.

#### How do I query the collection of (metagenomic) reads?
Similarly, you'll need to run ```./cammiq --query [options] [parameters]``` from command line, where ```[options]``` include
  * ```--read_cnts``` **Optional**. If ```--read_cnts``` is specified, then CAMMiQ will not produce its standard output (see below ```-o``` option). Instead, CAMMiQ outputs a non-negative (tab-delimited) matrix where each row corresponds to a query (fastq file); each column corresponds to an **NCBI taxonomic ID** (attention: not a genome ID!); each entry gives the number of reads in a given query that CAMMiQ assigned uniquely to the corresponding taxon. ```--read_cnts``` queries correspond to the simplest "Type I" queries in the paper.
  
    Here is an example output when CAMMiQ finds ```--read_cnts``` in a command line:
  
    ```
      7 9 11 14 19
    query_1.fq 0 0 15 100 10000
    query_2.fq 20 0 125 1800 0
    query_3.fq 0 0 0 0 0
    ......
    ```
  
    Note that the output file name can be specified with ```-o``` parameter. 
  
  * ```--doubly_unique``` **Optional**. Only valid when ```--read_cnts``` is specified. CAMMiQ will resolve the ambiguous read counts brought by doubly-unique substrings, and assign each of those reads that only contain doubly-unique substrings from two distinct taxa to one specific taxon. ```--doubly_unique``` queries correspond to "Type II" queries in the paper.

and ```[parameters]``` include the following list of (possibly mandatory) parameters.
* ```-f <MAP_FILE>``` **Mandatory**. You should use the same ```<MAP_FILE>``` when building the index in your queries. **Attention:** CAMMiQ is not in charge of verifying the format or correctness (meaning that you use exactly the same file for building the index and querying) of a ```<MAP_FILE>```. When your input ```<MAP_FILE>``` for querying is different from what you used for building the index, some potential "*undefined behavior*" could happen when running CAMMiQ.
* ```-q (-Q) <QUERY_FILE(S)>``` **Mandatory**. ```<QUERY_FILE(S)>``` can be either a list of fastq files, or a directory containing the list of fastq files in your query. A capitalized ```-Q``` indicate the input is a directory, that is,  
  * ```-q``` CAMMiQ takes the list of fastq files.
  * ```-Q``` CAMMiQ takes a directory which contains the list of fastq files.
* ```-i <INDEX_FILES>``` **Optional**. As per ```cammiq --build``` option, the default index file names (if ```-i``` is not specified) are ```index_u.bin1``` and ```index_d.bin2```. If index files are stored in a different loaction (directory), use ```-i``` to specify the directory and file names: as discussed in [What is CAMMiQ index composed of?](#-What-is-CAMMiQ-index-composed-of?), ```<INDEX_FILES>``` **should include** ```*.bin1``` and ```*.bin2``` plus additional meta-information; the default location for CAMMiQ to find those meta-information (```genome_lengths.out```,  ```unique_lmer_count_u.out```, ```unique_lmer_count_d.out```) files is the same as that provided in ```<INDEX_FILES>```.

  One exception is when ```--read_cnts``` is specified. In this case CAMMiQ will automatically *ignore* the meta-information files.

* ```-o <OUTPUT_FILE>``` **Optional (but strongly recommended)**. CAMMiQ's standard output file, with ```<OUTPUT_FILE>``` specifying the file name. Default output file name is ```quantification_results.out```. ```<OUTPUT_FILE>``` presents the (normalized) abundances resulted from each query fastq file. Each line (except the header lines) from left to right contains 3 tab-delimited fields: NCBI taxonomic ID, Abundance and Organism names. Two queries (fastqs) are separated by a blank line.
  
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
  
  **Attention.** Unlike some other commonly tools (which simply represent *read counts* as genome aubndances), CAMMiQ represents *sequencing depth* for each genome as its abundance. Consider the following example, you have 100 reads from genome (species) A and 100 reads from genome (species) B, but A has twice larger genome length than B - then CAMMiQ will output the abundance of A as 0.333 and the abundance of B as 0.667.

* ```-e <float>``` **Optional (but strongly recommended)**. Estimated average sequencing error (substitution) rate. Default value is ```e = 0.0``` (i.e., no sequencing errors); for typical Illumina reads, you can try something like ```-e 0.008``` or ```-e 0.01```. Of course, the more accurate you estimate this error rate in your queries, the more accurate quantification result (i.e. aubndances of the input genomes) CAMMiQ can produce.   
* ```-h <int>|<int1 int2>``` **Optional (but need special attention)**. Same as what is described above in the ```--build``` option, ```-h``` can take in one or two integer values. You don't need to speicify ```-h``` in ```--query``` option since it's encoded in the corredponding index (```*.bin1``` and ```*.bin2```) files - but if something is found here and not matching the value encoded in the index files, CAMMiQ will raise an error.
* ```-t <int>``` **Optional**. Number of threads used during CAMMiQ queries. Same as what is described above in the ```--build``` option. Note that CAMMiQ query has two phases: in the first phase of assigning reads to one or two genomes, the default number of threads is 1; however, in the second phase CAMMiQ will use IBM ILOG CPLEX Optimization Studio to produce metagenomic quantification results (or read count results if ```--read_cnts``` is specified), which, by default, is again 'auto-threaded' (i.e. attempts to use all available CPUs on a computer). If your Linux mahcine has more than 32 CPUs, then by default 32 threads will be used.

Additionally, you can also modify the following list of parameters.

* ```--enable_ilp_display```, which allows the ILP solvers (IBM ILOG CPLEX Optimization Studio or Gurobi Optimizer) to print their debug information.
* ```--read_length_filter <int>``` CAMMiQ will discard reads shorter than the specified value in ```<QUERY_FILE(S)>```. By default, CAMMiQ will keep all reads in the queries. 
However, users need to pay attention to the possible existence of reads shorter than the hash length specified in ```-h``` (or those built in CAMMiQ indices); in such case CAMMiQ could raise an error. 
* ```--read_cnt_thres <int>``` This is a firm threshold to claim the existence of a certain genome: if it contains sufficient number of unique substrings (or unique *L*-mers, specified by ```--easy_to_identify_thres``` below), then CAMMiQ will apply both the integer value specified here and the floating point value specified in ```--ilp_alpha``` (see below) - those genomes having less than these numbers of reads assigned to will not be output by CAMMiQ; if not (i.e., for those genomes that do not pass ```--easy_to_identify_thres```) then CAMMiQ will only apply the value given in ```--ilp_alpha``` to filter out non-existing genomes. The default value of ```--read_cnt_thres``` is 100.
* ```--easy_to_identify_thres <int>``` When a genome contains more unique *L*-mers than the value specified here, it is considered "easy to identify" and the firm ```--read_cnt_thres``` filter will be applied (see above). The default value is 10000.
* ```--ilp_epsilon <float>``` Default value is 0.01. If CAMMiQ does not give any *feasible* solution for its quantification results, then users should attempt to increase this value for generating at least one (more relaxed) feasible solution.
* ```--ilp_alpha <float>``` **Important**. CAMMiQ's relative read count filter (Intuitively, it can be considered as the "resolution" of CAMMiQ). If a genome has less than ```ilp_alpha * number of unique/doubly-unique L mers``` reads being assigned to, then it will be filtered out in the quantification process. The default value of ```--ilp_alpha``` is 0.0001. See CAMMiQ paper https://www.biorxiv.org/content/10.1101/2020.06.12.149245v2 for further details.
* ```--max_depth <float>``` Maximum allowed coverage per genome. The default value is 100.0.
* ```--unique_read_cnt_thres <int>``` A special parameter for ```--read_cnts (--doubly-unique)``` queries; it's mainly designed for single cell queries. In this case when a genome has sufficient number of reads assigned to (i.e., pass the threshold), it must be included in CAMMiQ's identification output.
* ```--doubly_unique_read_cnt_thres <int>``` A special parameter for ```--read_cnts (--doubly-unique)``` queries; it's also designed for single cell queries. When a genome has less than ```--unique_read_cnt_thres``` unique reads and less than ```--doubly_unique_read_cnt_thres``` doubly-unique reads (that is, these reads are ambiguously assigned to the corresponding genome), it must be excluded/filtered out in CAMMiQ's identification output.

**Example usage**:
```
cammiq --query --read_cnts -h 26 26 -f genome_map1.out -i index_u.bin1 index_d.bin2 -Q /data/fastqs/ --read_length_filter 70 -o cammiq_identification.txt
cammiq --query --read_cnts --doubly_unique -h 26 26 -f genome_map1.out -i index_u.bin1 index_d.bin2 -Q /data/fastqs/ -o cammiq_identification.txt --unique_read_cnt_thres 20 --doubly_unique_read_cnt_thres 10 -t 16
```
The first instruction queries against ```index_u.bin1``` and ```index_d.bin2``` the reads from all fastq files under the directory ```/data/fastqs/```; with each fastq forming a distinct query; with read length less than 70 being discarded; with the number of reads uniquely assigned to the genomes in ```genome_map1.out``` (which must correspond to the index files) written into ```cammiq_identification.txt```, for future usage.

The second instruction queries the same set of reads, but this time CAMMiQ considers the reads ambiguously assigned to two genomes and resolves such ambiguity through ILP solvers (IBM ILOG CPLEX Optimization Studio or Gurobi Optimizer), with parameter ```--unique_read_cnt_thres``` set to 20 and ```--doubly_unique_read_cnt_thres``` set to 10. The maximum allowed number of threads for ILP solvers is set to 16.

```
cammiq --query -f genome_map2.out -i ../index_u_test.bin1 ../index_d_test.bin2 -q query1.fastq query2.fastq -o cammiq_quantification.out -t 8
cammiq --query -h 31 -f genome_map3.out -i index_u.bin1 index_d.bin2 -Q /data/fastqs/ -o gurobi_solutions.out -t 8 --ilp_alpha 0.001
cammiq --query -h 25 25 -f genome_map4.out -i index_u.bin1 index_d.bin2 -Q /data/test_fastqs/ -o cplex_solutions.out -t 8 --enable_ilp_display --ilp_alpha 0.001 --read_cnt_thres 50 --ilp_epsilon 0.1
```
The third instruction queries against ```index_u_test.bin1``` and ```index_d_test.bin2``` (both stored in the parent level of directory ```../```; with meta-information also stored under ```../```) the reads from two given fastq files ```query1.fastq``` and ```query2.fastq```; with each fastq forming a distinct query; with the corresponding genomes (of course, within the list of ```genome_map2.out```) and their relative bundances written into ```cammiq_quantification.out```; with 8 threads both for querying the index and using ILP to quantify the abundances.

The fourth instruction queries against ```index_u.bin1``` and ```index_d.bin2``` the reads from all fastq files under the above directory ```/data/fastqs/```; with each fastq forming a distinct query; with the corresponding genomes (within the list of ```genome_map3.out```) and their relative bundances written into ```gurobi_solutions.out```; with again 8 threads both for querying the index and using ILP to quantify the abundances; with parameter ```--ilp_alpha``` set to ```0.001``` (this correspond to part of our results shown in the paper). Note that with ```-h 31``` set both hash lengths to ```31``` - in this case if the hash lengths stored in the indices do not match, CAMMiQ will report an error.

The fifth instruction queries against ```index_u.bin1``` and ```index_d.bin2``` the same set of fastqs; with each fastq forming a distinct query; with the corresponding genomes (within the list of ```genome_map4.out```) and their relative bundances written into ```cplex_solutions.out```; with CPLEX (ILP solver) debug information enabled on the screen printing output; and with multiple parameters tuned. Same attention needed to be paid on specified hash lengths.

#### How do I query single cell RNA-seq reads?
Since RNA-seq reads usually do not follow (roughly) uniform distribution across the entire genome, users are suggested to apply ```--read_cnts``` and/or ```--doubly_unique``` queries in command line options, which will output the number of reads identified though CAMMiQ indices which are assigned to each genome in the index database. Multiple parameters including e.g., ```--read_length_filter```; ```--(doubly_)unique_read_cnt_thres``` can be tuned to control CAMMiQ's output. 

### What is the expected computational cost of CAMMiQ?
Since both index construction and querying of CAMMiQ run in a reasonable amount time (usually a few hours with multiple threads for index construction and minutes for querying a few million reads; additionally CAMMiQ may require a few minutes to load the indices into main memory when you query the first fastq file), the main bottleneck comes from main memory usage. Here is a quick estimation of the memory needed for index construction: suppose your total size of input genomes is ```N```, then CAMMiQ will need to reserve up to ```37 * N``` Bytes (plus the memory for maintaining the index structures) of RAM in total during its index construction step (plus certain amount of disk space). For querying, suppose the total index size for ```*.bin1/*.bin1.aux``` and ```*.bin2/*.bin2.aux``` files is ```IDX_SIZE```, and the total size of query files is ```QRY_SIZE```, then CAMMiQ will need roughly ```8 * IDX_SIZE + QRY_SIZE``` Bytes of RAM in total during its querying step, plus the memory required by the ILP solvers.

