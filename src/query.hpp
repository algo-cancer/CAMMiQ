#ifndef FQREADER_HPP_
#define FQREADER_HPP_

#include <cassert>
#include <pthread.h>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>

#include "hashtrie.hpp"

struct Genome {
	public:
		uint64_t read_cnts_u;
		uint64_t read_cnts_d;
		uint32_t glength;
		uint32_t nus;
		uint32_t nds;
		uint32_t taxID;
		std::string name;

		Genome(uint32_t, std::string&);
		~Genome() {};
};

class FqReader {
	private:
		/* Fastq file names. */
		std::vector<std::string> qfilenames;
		std::string current_filename;

		/* Reads. */
		size_t max_rl = 256;
		std::vector<std::vector<uint8_t*>> reads;
		std::vector<std::vector<uint8_t>> rlengths;
		std::vector<size_t> tlengths;
		
		/* Number of "conflict" (having more than one refID) reads. */
		size_t nconf = 0;
		
		/* Number of "unlabeled" (cannot assign one refID) reads. */
		size_t nundet = 0;
		
		/* Meta informations. */
		std::string MAPFILE;
		std::vector<Genome*> genomes;
		int *exist = NULL;
		std::map<std::pair<uint32_t, uint32_t>, uint64_t> read_cnts_b;

		/* Hash parameters. */
		uint32_t hash_len_u = -1;
		uint32_t hash_len_d = -1;
		std::string IDXFILEU;
		std::string IDXFILED;
		Hash *ht_u = NULL;
		Hash *ht_d = NULL;
		float erate_;
		std::string OUTPUTFILE;

		bool debug_display = 0;

	public:
		int nthreads = 1;
		
		FqReader(std::string&, std::string&, std::string&, std::string&, float, bool);
		FqReader(uint32_t, std::string&, std::string&, std::string&, std::string&, float, bool);
		FqReader(uint32_t, std::string&, uint32_t, std::string&, std::string&, std::string&, float, bool);
		~FqReader();

		void clearReads();

		void loadIdx();
		void loadIdx_p();
		void* loadIdx_u_() {
			ht_u->loadIdx64_p(IDXFILEU);
			if (hash_len_u > 0)
				assert(ht_u->getHashLength() == hash_len_u);
			else
				hash_len_u = ht_u->getHashLength();
			return 0;
		}
		void* loadIdx_d_() {
			ht_d->loadIdx64_p(IDXFILED);
			if (hash_len_d > 0)
				assert(ht_d->getHashLength() == hash_len_d);
			else
				hash_len_d = ht_d->getHashLength();
			return 0;
		}
		static void *loadIdx_u(void *obj) {
			return ((FqReader*) obj)->loadIdx_u_();
		}
		static void *loadIdx_d(void *obj) {
			return ((FqReader*) obj)->loadIdx_d_();
		}

		void loadSmap();
		void loadGenomeLength();

		/* FASTQ Operations. */
		void getFqList(std::string&);
		void readFastq(std::string&, size_t);
		void readFastq(std::string&, size_t, size_t);
		void readallFastq();
		void readallFastq(size_t);
		void prepallFastq();
		//void prepallFastq_sc();

		void getFqnameWithoutDir(size_t);

		void query64_p(size_t);
		void query64mt_p(size_t);
		void query64_sc(size_t);

		void queryFastq_p(std::string&, size_t, std::vector<double>&);
		void queryFastq_p(std::vector<std::string>&, size_t, std::vector<double>&);
		void queryFastq_sc(int, std::string&, size_t, std::vector<double>&);
		void queryFastq_sc(int, std::vector<std::string>&, size_t, std::vector<double>&);

		void getRC(uint8_t*, uint8_t*, size_t);

		void runILP_cplex(size_t, int, uint32_t, double, double, double, double);
		void runILP_gurobi(size_t, int, uint32_t, double, double, double, double);
		void runILPsc_cplex(size_t, uint32_t, uint32_t);
		void runILPsc_gurobi(size_t, uint32_t, uint32_t);

		void outputUniqueCnts(size_t);

		void resetCounters();
		void resetCounters_sc();
		
		static int symbolIdx[256];

		static int rcIdx[128];

		static char alphabet[4];
};

#endif
