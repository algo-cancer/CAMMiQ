### What is CAMMiQ?

CAMMiQ is a software tool for microbial identification and quantification. Specifically, it builds a compact index to manage a database of microbial genomes (in particular, bacterial/viral genomes in NCBI [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) database) and answers queries of the following form: given a set of high throughput sequencing (HTS) reads obtained from a mixture of genomes in the database, identify the corresponding genomes (in the database) and compute their relative abundance.

### How to install CAMMiQ?
Dependencies: Though several scripts (to compile the code and download the database) are written using the Bash shell and python, our core programs to build the index and query the set of reads are written in C++11, and need to be compiled using a relatively recent version (The current build of CAMMiQ was tested with version 5.2.0) of gcc that will support C++11. Multithreading is handled using pthread and [OpenMP](https://en.wikipedia.org/wiki/OpenMP).

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
To begin using CAMMiQ, you will first need to install it, and then either download or create a database.
