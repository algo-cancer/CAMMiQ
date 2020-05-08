#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include <cstdlib>

#include "util.hpp"
#include "build.hpp"
#include "query.hpp"

bool validFile(const char* _file) {
        FILE * fd = fopen(_file, "r");
        if (fd == NULL)
		return false;
        fclose(fd);
        return true;
}

void printUsage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "CAMMiQ: Metagenomic microbial abundance quantification.\n\n");
	fprintf(stderr, "Usage: ./cammiq --<option> parameters\n\n");
	fprintf(stderr, "Definitions of parameters: \n\n");
	fprintf(stderr, "option = build | query .\n");
	fprintf(stderr, "-k <integer>,\t minimum substring length: integer, >= 5 and <= read length. Default value is 26.\n");
	fprintf(stderr, "-L <integer>,\t read length: integer, default value is 100.\n");
	fprintf(stderr, "-h <integer>,\t hash length: integer, default value is 26.\n");
	fprintf(stderr, "-f <string>,\t file names separated by space.\n");
	fprintf(stderr, "-d <string>,\t directoty containing the fasta or fastq files. \n");
	fprintf(stderr, "-i <string>,\t indexing options: one in unique | doubly_unique | both.\n");
	fprintf(stderr, "-t <integer>,\t number of threads. \n\n");
	fprintf(stderr, "--help,\t to print user options.\n");
	//fprintf(stderr, "--version,\t to print the version information.\n\n");
	//fprintf(stderr, "Version: 0.0; Contact: kzhu@indiana.edu.\n");
}

int main(int argc, char** argv) {
	if (argc == 2) {
		std::string val(argv[1]);
		if (val == "--help" || val == "--HELP" ) {
			printUsage();
			return 0;
		}
		exit(EXIT_FAILURE);
	}

	/* Check the parameters. */
	int K = 26, L = 100, t = 1, h = -1, h1 = -1, h2 = -1;
	std::string fa_name = "", fa_dir = "", fm_name = "", fi_name = "", fq_name = "", fq_dir = "";
	std::vector<std::string> fa_names;
	std::vector<std::string> fq_names;
	FastaReader *main_fr = NULL;
	FqReader *main_fqr = NULL;
	int mode = 0;
	std::string idx_option;
	std::string output;
	float erate = 0.01;
	int u = 0, hash_option = 32;
	std::string i1fn, i2fn;
	for (int i = 1; i < argc; i++) {
		std::string val(argv[i]);
		if (val == "--build") {
			mode = 0;
			continue;
		}
		if (val == "--query") {
			mode = 1;
			continue;
		}
		if (val == "-k") {
			if (mode > 0) {
				fprintf(stderr, "Parameter k is only valid in mode BUILD.\n"); 
				exit(EXIT_FAILURE);
			}
			if (++i >= argc) {
				fprintf(stderr, "Please specify the minimum substring length.\n"); 
				exit(EXIT_FAILURE);
			}
			K = atoi(argv[i]); 
			if (K <= 4 || K > maxK) {
				fprintf(stderr, "The min substring length should be in range [5, %d].\n", maxK); 
				exit(EXIT_FAILURE);
			}
			continue;
		}
		if (val == "-L") {
			if (mode > 0) {
				fprintf(stderr, "Parameter L is only valid in mode BUILD.\n"); 
				exit(EXIT_FAILURE);
			}
			if (++i >= argc) {
				fprintf(stderr, "Please specify the read length / maximum unique substring length.\n"); 
				exit(EXIT_FAILURE);
			}
			L = atoi(argv[i]);
			if (L < K || L > maxL) {
				fprintf(stderr, "The max substring length should be in range [%d, %d].\n", K, maxL); 
				exit(EXIT_FAILURE);
			}
			continue;
		}
		if (val == "-i") {
			if (++i >= argc) {
				fprintf(stderr, "Please specify indexing options.\n"); 
				exit(EXIT_FAILURE);
			}
			idx_option = argv[i];
			continue;
		}
		if (val == "-o") {
			output = argv[++i];
			continue;
		}
		if (val == "-e") {
			erate = std::stof(argv[++i]);
			continue;
		}

		if (val == "-h") {
			if (++i >= argc) {
				fprintf(stderr, "Please specify hash length as an integer.\n"); 
				exit(EXIT_FAILURE);
			}
			if (i + 1 < argc && argv[i + 1][0] != '-') {
				if (mode == 0) {
					fprintf(stderr, "Multiple values of h is only valid in mode QUERY.\n"); 
					exit(EXIT_FAILURE);
				}
				h1 = atoi(argv[i++]); 
				h2 = atoi(argv[i]);
			} else
				h = atoi(argv[i]);
			if ((h != -1) && (h <= 4 || h >= 32)) {
				fprintf(stderr, "The hash length should be in range [5, 31].\n"); 
				exit(EXIT_FAILURE);
			}
			if ((h1 != -1) && (h1 <= 4 || h1 >= 32)) {
				fprintf(stderr, "The hash length should be in range [5, 31].\n"); 
				exit(EXIT_FAILURE);
			}
			if ((h2 != -1) && (h2 <= 4 || h2 >= 32)) {
				fprintf(stderr, "The hash length should be in range [5, 31].\n"); 
				exit(EXIT_FAILURE);
			}
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
				fprintf(stderr, "Please specify file names.\n");
				exit(EXIT_FAILURE);
			}
			switch (mode) {
				case 1:
					while (i < argc && argv[i][0] != '-') {
						std::string filename = argv[i++];
						//fprintf(stderr, "%s", filename.c_str());
						std::string ext = filename.substr(filename.find_last_of(".") + 1);
						if (ext == "fq" || ext == "fastq") {
							fq_name = filename;
							if (validFile(fq_name.c_str()))
								fq_names.push_back(fq_name);
							else {
								fprintf(stderr, "Failed to find input file %s.\n", fq_name.c_str());
								exit(EXIT_FAILURE);
							}
						}
						if (ext == "out" || ext == "map")
							fm_name = filename;
						if (ext == "idx" || ext == "bin")
							fi_name = filename;
						if (ext == "idx1" || ext == "bin1")
							i1fn = filename;
						if (ext == "idx2" || ext == "bin2")
							i2fn = filename;
					}
					i--;
					break;
				default:
					while (i < argc && argv[i][0] != '-') {
						std::string filename = argv[i++];
						std::string ext = filename.substr(filename.find_last_of(".") + 1);
						if (ext == "fna" || ext == "fasta" || ext == "ffn") {
							fa_name = filename;
							if (validFile(fa_name.c_str()))
								fa_names.push_back(fa_name);
							else {
								fprintf(stderr, "Failed to find input file %s.\n", fa_name.c_str());
								exit(EXIT_FAILURE);
							}
						}
						if (ext == "out" || ext == "map") {
							fm_name = filename;
						}
					}
					i--;
					break;
			}
			continue;
		}
		if (val == "-d") {
			switch (mode) {
				case 1: //fastq directory
					if (++i >= argc) {
						fprintf(stderr, "Please specify the directory containing fastq files.\n"); 
						exit(EXIT_FAILURE);
					}
					fq_dir = argv[i];
					if (!validFile(fq_dir.c_str())) {
						fprintf(stderr, "Failed to find input directory %s.\n", fq_dir.c_str());
						exit(EXIT_FAILURE);
					}
					break;
				default: // fasta directory
					if (++i >= argc) {
						fprintf(stderr, "Please specify the directory containing fasta files.\n"); 
						exit(EXIT_FAILURE);
					}
					fa_dir = argv[i];
					if (!validFile(fa_dir.c_str())) {
						fprintf(stderr, "Failed to find input directory %s.\n", fa_dir.c_str());
						exit(EXIT_FAILURE);
					}
					break;
			}
			continue;
		}
		fprintf(stderr, "Failed to recognize option: %s. \n", val.c_str());
		exit(EXIT_FAILURE);
	}

	if (h == -1 && h1 == -1 && h2 == -1) {
		fprintf(stderr, "Please specify a valid hash length.\n");
		exit(EXIT_FAILURE);
	}
	
	switch (mode) {
		case 0:
			if (!fa_names.empty()) {
				if (fa_dir.length() > 0)
					fprintf(stderr, "Ignoring the input directory.\n");
				main_fr = new FastaReader(L, K, t, fa_names);
			} else
				main_fr = new FastaReader(L, K, t);

			main_fr->readFnMap(fa_dir, fm_name);
			main_fr->readAllFasta();
			main_fr->allocSuffixArray(0);

			if (idx_option.compare("unique") == 0) {
				main_fr->allocSuffixArray(1);
				main_fr->setHashLength(h);
				main_fr->computeIndex(1);
			}
			if (idx_option.compare("doubly_unique") == 0) {
				main_fr->allocSuffixArray(2);
				main_fr->setHashLength(h);
				main_fr->computeIndex(2);
			}
			if (idx_option.compare("both") == 0) {
				main_fr->allocSuffixArray(1);
				main_fr->setHashLength(h1);
				main_fr->computeIndex(1);
				main_fr->allocSuffixArray(2);
				main_fr->setHashLength(h2);
				main_fr->computeIndex(2);
			}
			if (main_fr)
				delete main_fr;
			break;
		case 1:
			//fprintf(stderr, "fm name: %s.\n", fm_name.c_str());
			if (fi_name != "") {
				u = (idx_option.compare("doubly_unique") == 0) ? 1 : 0;
				if (h == -1) {
					fprintf(stderr, "Please specify a valid hash length.\n");
					exit(EXIT_FAILURE);
				}
				hash_option = (h >= 16) ? 64 : 32;
				main_fqr = new FqReader(u, hash_option, h, fi_name, fm_name);
			} else {
				if (h1 == -1 || h2 == -1) {
					fprintf(stderr, "Please specify a valid hash length.\n");
					exit(EXIT_FAILURE);
				}
				main_fqr = new FqReader(32, h1, i1fn, h2, i2fn, fm_name, output, erate);
				//assert(h1 == h2);
				//main_fqr = new FqReader(32, h1, i1fn, i2fn, fm_name, output, erate);
			}
			//main_fqr->loadIdx();
			main_fqr->loadIdx_p();
			main_fqr->loadSmap();
			if (!fq_names.empty())
				main_fqr->queryFastq_p(fq_names);
				//main_fqr->queryFastq(fq_names);
			else {
				main_fqr->getFqList(fq_dir);
				main_fqr->queryAllFastq();
			}
			break;
		case 2:
			break;
	}

	return 0;
}
