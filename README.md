# What is CAMMiQ?

CAMMiQ is a software tool for microbial identification and quantification. Specifically, it builds an index to manage a database of microbial genomes and 
answers queries of the following form: given a set of HTS reads obtained from a mixture of genomes in the database,
identify the corresponding genomes in the database and compute their relative abundance.

# How to install CAMMiQ?
Dependencies: Though several scripts are written using the Bash shell and python, core programs needed to build the database and run the are written in C++11, and need to be compiled using a somewhat recent version of g++ that will support C++11. Multithreading is handled using OpenMP. Downloads of NCBI data are performed by wget and rsync.

To install CAMMiQ, just clone our repository and run the script install_CAMMiQ.sh we have prepared :

# How do I use CAMMiQ?
To begin using CAMMiQ, you will first need to install it, and then either download or create a database.
