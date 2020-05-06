#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include <cstdlib>
//#include <climits>

#include "util.hpp"
//#include "gsa.hpp"
#include "build.hpp"

bool validFile(const char* _file) {
        FILE * fd = fopen(_file, "r");
        if (fd == NULL)
		return false;
        fclose(fd);
        return true;
}

void printUsage() {
	fprintf(stdout, "\n");
	fprintf(stdout, "Fast metagenomic species abundance estimation.\n\n");
	fprintf(stdout, "Usage: ./t -k <integer> -L <integer> -f <fasta files> -d <directory>\n\n");
	fprintf(stdout, "Definitions of parameters: \n\n");
	fprintf(stdout, "-k <integer>,\t minimum substring length: integer, >= 5 and <= read length. The default value is 20.\n");
	fprintf(stdout, "-L <integer>,\t read length: integer, the default value is 100.\n");
	fprintf(stdout, "-f <string>,\t fasta file names separated by space.\n");
	fprintf(stdout, "-d <string>,\t directoty containing the fasta files. \n");
	fprintf(stdout, "-t <integer>,\t number of threads. \n\n");
	fprintf(stdout, "--help,\t to print user options.\n");
	fprintf(stdout, "--version,\t to print the version information.\n\n");
	fprintf(stdout, "Version: 0.0; Contact: kzhu@indiana.edu.\n");
}

int main(int argc, char** argv) {
	if (argc == 2) {
		std::string val(argv[1]);
		if (val == "--help" || val == "--HELP" ) {
			printUsage();
			return 0;
		}
	}
	if (argc <= 6) {
		fprintf(stderr, "At least 3 parameters are required.\n");
		printUsage();
		exit(EXIT_FAILURE);
	}

	/* Check the parameters. */
	int K = 20, L = 100, t = 1;
	std::string fa_name = "", fa_dir = "", fm_name = "../Human_Gut/species_map_gut.out";
	//std::string fa_name = "", fa_dir = "", fm_name = "./s1.out";
	std::vector<std::string> fa_names;
	FastaReader *main_fr = NULL;
	std::string mode;
	for (int i = 1; i < argc; i++) {
		std::string val(argv[i]);
		if (val == "-k") {
			if (++i >= argc) {
				fprintf(stderr, "Please specify the minimum unique substring length.\n"); 
				exit(EXIT_FAILURE);
			}
			K = atoi(argv[i]); 
			if (K <= 4 || K > maxK) {
				fprintf(stderr, "The mus length should be in range [5, %d].\n", maxK); 
				exit(EXIT_FAILURE);
			}
			continue;
		}
		if (val == "-L") {
			if (++i >= argc) {
				fprintf(stderr, "Please specify the read length / maximum unique substring length.\n"); 
				exit(EXIT_FAILURE);
			}
			L = atoi(argv[i]); 
			if (L < K || L > maxL) {
				fprintf(stderr, "The max us length should be in range [%d, %d].\n", K, maxL); 
				exit(EXIT_FAILURE);
			}
			continue;
		}
		if (val == "-m") {
			if (++i >= argc) {
				//fprintf(stderr, "Please specify.\n"); 
				exit(EXIT_FAILURE);
			}
			mode = argv[i];
			continue;
		}
		if (val == "-t") {
			if (++i >= argc) {
				fprintf(stderr, "Please specify the worker threads number.\n"); 
				exit(EXIT_FAILURE);
			}
			t = atoi(argv[i]); 
			continue;
		}
		if (val == "-f") {
			if (++i >= argc) {
				fprintf(stderr, "Please specify the fasta file names.\n");
				exit(EXIT_FAILURE);
			}
			while (i < argc && argv[i][0] != '-') {
				fa_name = argv[i++];
				if (validFile(fa_name.c_str()))
					fa_names.push_back(fa_name);
				else
					fprintf(stderr, "Failed to find input file %s.\n", fa_name.c_str());
			}
			if (fa_names.empty()) {
				fprintf(stderr, "Please specify valid fasta file names.\n");
				exit(EXIT_FAILURE);	
			}
			i--;
			continue;
		}
		if (val == "-d") {
			if (++i >= argc) {
				fprintf(stderr, "Please specify the directory containing fasta files.\n"); 
				exit(EXIT_FAILURE);
			}
			fa_dir = argv[i];
			if (!validFile(fa_dir.c_str())) {
				fprintf(stderr, "Failed to find input directory %s.\n", fa_dir.c_str());
				exit(EXIT_FAILURE);
			}
			continue;
		}
		fprintf(stderr, "Failed to recognize option: %s. \n", val.c_str());
		exit(EXIT_FAILURE);
	}

	if (!fa_names.empty()) {
		if (fa_dir.length() > 0)
			fprintf(stderr, "Ignoring the input directory.\n");
		main_fr = new FastaReader(L, K, t, fa_names);
		//main_fr->readAllFasta();
	} else {
		if (fa_dir.length() > 0) {
			/*std::string test = "GGACCCCTAATTA";
			uint8_t *test_ = new uint8_t[26];
			for (size_t i = 0; i < test.length(); i++)
				test_[i] = test[i] + 165;
			for (size_t i = 0; i < test.length(); i++)
				test_[test.length() + i] = test[i] + 165;
			Hash *hash_test = new Hash(13);
			uint32_t bucket = hash_test->computeHashVal(test_);
			hash_test->insert32_d(bucket, test_, 13, 1);
			
			hash_test->insert32_d(bucket, test_, 13, 1);*/
			
			main_fr = new FastaReader(L, K, t);
			//main_fr->prepFasta(fa_dir);
			main_fr->readFnMap(fa_dir, fm_name);
			main_fr->readAllFasta();
			
			main_fr->allocSuffixArray();
			main_fr->setHashLength(11);
			if (mode.compare("stat") == 0)
				main_fr->computeIndex(0);
			if (mode.compare("non_overlap") == 0)
				main_fr->computeIndex(1);
			if (mode.compare("minimum") == 0)
				main_fr->computeIndex(2);
			if (mode.compare("doubly_unique") == 0)
				main_fr->computeIndex(4);
			if (mode.compare("stat_clark") == 0)
				main_fr->computeIndex(3);}
			//main_fr->analIndex(K, L);
			//main_fr->outputRefLength();
	}
	if (main_fr)
		delete main_fr;

	return 0;
}
