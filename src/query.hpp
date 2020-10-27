#ifndef FQREADER_HPP_
#define FQREADER_HPP_

#include <pthread.h>
#include <vector>
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

		/* Reads. */
		size_t max_rl = 256;
		std::vector<uint8_t*> reads;
		
		/* Number of "conflict" (having more than one refID) reads. */
		size_t nconf = 0;
		
		/* Number of "unlabeled" (cannot assign one refID) reads. */
		size_t nundet = 0;
		
		/* Number of reads */
		std::string MAPFILE;
		std::vector<Genome*> genomes;
		int *exist = NULL;

		/* Hash parameters. */
		int doubly_unique = 0;
		int hash_option = 32;
		uint32_t hash_len_u = 12;
		uint32_t hash_len_d = 12;
		std::string IDXFILEU;
		std::string IDXFILED;
		Hash *ht_u = NULL;
		Hash *ht_d = NULL;
		float erate_;
		std::string OUTPUTFILE;

	public:
		int tid_ = 0;
		pthread_mutex_t thread_lock;
		
		FqReader(int, int, uint32_t, std::string&, std::string&);
		FqReader(int, uint32_t, std::string&, uint32_t, std::string&, std::string&);
		FqReader(int, uint32_t, std::string&, uint32_t, std::string&, std::string&, std::string&, float);
		FqReader(int, uint32_t, std::string&, std::string&, std::string&, std::string&, float);
		~FqReader();

		void resetReads();

		void loadIdx();
		void loadIdx_p();
		void* loadIdx_p__();
		static void *loadIdx_p_(void *obj) {
			return ((FqReader*) obj)->loadIdx_p__();
		}

		void loadSmap();
		void loadGenomeLength();

		void getFqList(std::string&);

		void readFastq(std::string&);

		//void query32_s(size_t);
		//void query64_s(size_t);
		// void query32_p(size_t);
		void query64_p(size_t);

		//void queryFastq(std::vector<std::string>&);
		void queryFastq_p(std::vector<std::string>&);

		//void queryAllFastq();

		void getRC(uint8_t*, uint8_t*, size_t);

		void runILP_p(int, int, uint32_t, double, double, double, double);
		
		static int symbolIdx[256];

		static int rcIdx[128];
};

#endif
