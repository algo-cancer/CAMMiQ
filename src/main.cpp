#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include <cstdlib>
//#include <cassert>
#include <omp.h>

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
	fprintf(stderr, "Usage: ./cammiq --<options> (--<options_for_build> or --<options_for_query>) parameters\n\n");
	fprintf(stderr, "Definitions of parameters: \n\n");
	fprintf(stderr, "options = build | query.\n");
	fprintf(stderr, "options_for_build = unique | doubly_unique | both.\n");
	fprintf(stderr, "options_for_query = read_cnts | doubly_unique.\n");
	fprintf(stderr, "-k <integer>,\t minimum substring length: integer, >= 5 and <= read length. Default value is 26.\n");
	fprintf(stderr, "-Lmax <integer>,\t maximum substring length: integer, > k and <= read length. Default value is 50.\n");
	fprintf(stderr, "-L <integer>,\t read length: integer, default value is 100.\n");
	fprintf(stderr, "-h <integers>,\t hash length(s): integer, default value is 26.\n");
	fprintf(stderr, "-f <strings>,\t map file name and possibly other file names separated by space.\n");
	fprintf(stderr, "-D <string>,\t directoty containing the fasta files.\n");
	fprintf(stderr, "-Q <string>,\t directoty containing the fastq files.\n");
	fprintf(stderr, "-q <strings>,\t querying file names separated by space.\n");
	fprintf(stderr, "-i <strings>,\t indexing file names separated by space.\n");
	fprintf(stderr, "-o <string>,\t output file name.\n");
	fprintf(stderr, "-e <float>,\t expected sequencing error probability in queries.\n");
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
	int K = -1, Lmax = -1, L = -1, t = 1, h = -1, h1 = -1, h2 = -1;
	std::string fa_name = "", fa_dir = "", fm_name = "", fi_name1 = "", fi_name2 = "", fq_name = "", fq_dir = "";
	std::vector<std::string> fa_names;
	std::vector<std::string> fq_names;
	FastaReader *main_fr = NULL;
	FqReader *main_fqr = NULL;
	int mode = -1, id_mode = 0;
	bool debug_info = 0;
	std::string idx_option;
	std::string output;
	float erate = 0.01;
	std::vector<double> fine_parameters;
	//int u = 0, hash_option = 32;
	//std::string i1fn, i2fn;
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
		if (val == "--unique") {
			if (mode > 0) {
				fprintf(stderr, "Option --unique is only valid in mode BUILD.\n"); 
				exit(EXIT_FAILURE);
			}
			idx_option = "unique";
			continue;
		}
		if (val == "--doubly_unique") {
			if (mode > 0) {
				fprintf(stderr, "Option --unique is only valid in mode BUILD.\n"); 
				exit(EXIT_FAILURE);
			}
			idx_option = "doubly_unique";
			continue;
		}
		if (val == "--both") {
			if (mode > 0) {
				fprintf(stderr, "Option --unique is only valid in mode BUILD.\n"); 
				exit(EXIT_FAILURE);
			}
			idx_option = "both";
			continue;
		}
		if (val == "--read_cnts") {
			if (mode <= 0) {
				fprintf(stderr, "Option --read_cnts is only valid in mode QUERY.\n"); 
				exit(EXIT_FAILURE);
			}
			id_mode = 1;
			continue;
		}
		if (val == "--enable_ilp_display") {
			if (mode <= 0) {
				fprintf(stderr, "Option --enable_ilp_display is only valid in mode QUERY.\n"); 
				exit(EXIT_FAILURE);
			}
			//i++;
			debug_info = 1;
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
			if (K < 0) {
				fprintf(stderr, "Please first specify a valid value of k.\n"); 
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
		if (val == "-Lmax") {
			if (mode > 0) {
				fprintf(stderr, "Parameter Lmax is only valid in mode BUILD.\n"); 
				exit(EXIT_FAILURE);
			}
			if (L < 0) {
				fprintf(stderr, "Please first specify a valid value of L.\n"); 
				exit(EXIT_FAILURE);
			}
			if (++i >= argc) {
				fprintf(stderr, "Please specify the maximum substring length.\n"); 
				exit(EXIT_FAILURE);
			}
			Lmax = atoi(argv[i]); 
			if (Lmax < K + 1 || Lmax > L) {
				fprintf(stderr, "The max substring length should be in range [%d, %d].\n", K + 1, L); 
				exit(EXIT_FAILURE);
			}
			continue;
		}
		if (val == "-i") {
			if (mode <= 0) {
				fprintf(stderr, "Parameter i is only valid in mode QUERY.\n"); 
				exit(EXIT_FAILURE);
			}
			if (++i >= argc) {
				fprintf(stderr, "Please specify index file names.\n"); 
				exit(EXIT_FAILURE);
			}
			while (i < argc && argv[i][0] != '-') {
				std::string filename = argv[i++];
				std::string ext = filename.substr(filename.find_last_of(".") + 1);
				//if (ext == "idx" || ext == "bin")
				//	fi_name = filename;
				if (ext == "idx1" || ext == "bin1")
					fi_name1 = filename;
				if (ext == "idx2" || ext == "bin2")
					fi_name2 = filename;
			}
			i--;
			continue;
		}
		if (val == "-o") {
			if (mode <= 0) {
				fprintf(stderr, "Parameter o is only valid in mode QUERY.\n"); 
				exit(EXIT_FAILURE);
			}
			if (++i >= argc) {
				fprintf(stderr, "Please specify an output file name.\n"); 
				exit(EXIT_FAILURE);
			}
			output = argv[i];
			continue;
		}
		if (val == "-e") {
			if (mode <= 0) {
				fprintf(stderr, "Parameter e is only valid in mode QUERY.\n"); 
				exit(EXIT_FAILURE);
			}
			if (++i >= argc) {
				fprintf(stderr, "Please specify the expected sequencing error rate.\n"); 
				exit(EXIT_FAILURE);
			}
			erate = std::stof(argv[i]);
			if (erate < 0.0 || erate > 0.2) {
				fprintf(stderr, "The error rate should be in range [0, 0.2].\n"); 
				exit(EXIT_FAILURE);
			}
			continue;
		}
		if (val == "-h") {
			if (++i >= argc) {
				fprintf(stderr, "Please specify hash length as an integer.\n"); 
				exit(EXIT_FAILURE);
			}
			if (i + 1 < argc && argv[i + 1][0] != '-') {
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
			if ((mode == 0) && ((h > K) || (h1 > K) || (h2 > K))) {
				fprintf(stderr, "The hash length should be less than or equal to k.\n"); 
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
			omp_set_num_threads(atoi(argv[i])); 
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
						std::string ext = filename.substr(filename.find_last_of(".") + 1);
						if (ext == "out" || ext == "map")
							fm_name = filename;
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
		if (val == "-q") {
			if (++i >= argc) {
				fprintf(stderr, "Please specify query file names.\n");
				exit(EXIT_FAILURE);
			}
			while (i < argc && argv[i][0] != '-') {
				std::string filename = argv[i++];
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
			}
			i--;
			continue;
		}
		if (val == "-Q") {
			if (++i >= argc) {
				fprintf(stderr, "Please specify the directory containing fastq files.\n"); 
				exit(EXIT_FAILURE);
			}
			fq_dir = argv[i];
			if (!validFile(fq_dir.c_str())) {
				fprintf(stderr, "Failed to find input directory %s.\n", fq_dir.c_str());
				exit(EXIT_FAILURE);
			}
			continue;
		}
		if (val == "-D") {
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
	
	if (h == -1 && h1 == -1 && h2 == -1) {
		fprintf(stderr, "Warning: missing parameter h, set to default h = 26.\n");
		h = 26;
		h1 = 26;
		h2 = 26;
	}
	
	switch (mode) {
		case 0:
			//NEED TO CHECK map file!
			if (K == -1) {
				fprintf(stderr, "Warning: missing parameter k, set to default k = 26.\n");
				K = 26;
			}
			if (L == -1) {
				fprintf(stderr, "Warning: missing parameter L, set to default L = 100.\n");
				L = 100;
			}
			if (Lmax == -1) {
				fprintf(stderr, "Warning: missing parameter Lmax, set to default Lmax = 50.\n");
				Lmax = 50;
			}
			if (!fa_names.empty()) {
				if (fa_dir.length() > 0)
					fprintf(stderr, "Ignoring the input directory.\n");
				main_fr = new FastaReader(L, Lmax, K, t, fa_names);
			} else
				main_fr = new FastaReader(L, Lmax, K, t);

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
				main_fr->allocSuffixArray(3);
				if (h == -1)
					main_fr->setHashLength(h1);
				else
					main_fr->setHashLength(h);
				main_fr->computeIndex(1);
				main_fr->allocSuffixArray(2);
				if (h == -1)
					main_fr->setHashLength(h2);
				else
					main_fr->setHashLength(h);
				main_fr->computeIndex(3);
			}
			if (main_fr)
				delete main_fr;
			break;
		case 1:
			if (h == -1) {
				if (h1 == -1 || h2 == -1) {
					fprintf(stderr, "Warning: Hash length not specified, using that encoded in the index.\n");
					main_fqr = new FqReader(fi_name1, fi_name2, fm_name, output, erate, debug_info);
				}
				main_fqr = new FqReader(h1, fi_name1, h2, fi_name2, fm_name, output, erate, debug_info);
			} else {
				main_fqr = new FqReader(h, fi_name1, fi_name2, fm_name, output, erate, debug_info);
			}
			main_fqr->loadIdx_p();
			main_fqr->loadSmap();
			main_fqr->nthreads = t;
			if (!fq_names.empty()) {
				if (id_mode == 0)
					main_fqr->queryFastq_p(fq_names, 0);
				else
					main_fqr->queryFastq_sc(fq_names, 70);
			} else {
				if (fq_dir != "") {
					if (id_mode == 0)
						main_fqr->queryFastq_p(fq_dir, 0);
					else
						main_fqr->queryFastq_sc(fq_dir, 70);
				} else {
					fprintf(stderr, "Please specify at least one query file or directory.\n");
					exit(EXIT_FAILURE);
				}
			}
			break;
		default:
			break;
	}

	return 0;
}
