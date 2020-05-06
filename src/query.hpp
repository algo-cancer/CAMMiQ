#ifndef FQREADER_HPP_
#define FQREADER_HPP_

#include <vector>
#include <unordered_set>
#include <unordered_map>

#include "hashtrie.hpp"

class FqReader {
	private:
		/* Fastq file names. */
		std::vector<std::string> qfilenames;

		/* Reads. */
		size_t max_rl = 0;
		std::vector<uint8_t*> reads;
		//std::vector<uint32_t> true_labels;
		
		/* Number of "conflict" (having more than one refID) reads. */
		size_t nconf = 0;
		
		/* Number of "unlabeled" (cannot assign one refID) reads. */
		size_t nundet = 0;
		
		/* Number of reads */
		std::string MAPFILE;
		std::unordered_map<uint32_t, size_t> species_order;
		std::unordered_map<uint32_t, uint64_t> read_cnts_u;
		std::unordered_map<uint32_t, uint64_t> read_cnts_d;
		std::unordered_map<uint32_t, uint32_t> genome_lengths;
		std::unordered_map<uint32_t, uint32_t> species_nus;
		std::unordered_map<uint32_t, uint32_t> species_nds;
		//std::unordered_map<uint32_t, uint32_t> species_tus;
		//std::unordered_map<uint32_t, uint32_t> species_tds;
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
		FqReader(int, int, uint32_t, std::string&, std::string&);
		FqReader(int, uint32_t, std::string&, uint32_t, std::string&, std::string&);
		FqReader(int, uint32_t, std::string&, uint32_t, std::string&, std::string&, std::string&, float);
		FqReader(int, uint32_t, std::string&, std::string&, std::string&, std::string&, float);
		~FqReader();

		void resetReads();

		void loadIdx();
		void loadIdx_p();
		//void loadIdx_g();

		void loadSmap();
		void loadGenomeLength();

		void getFqList(std::string&);

		void readFastq(std::string&);

		void query32();
		void query64();
		//void quickQuery()

		void query32_d();
		void query64_d();

		//void query32_m();
		//void query64_m();

		//std::unordered_set<pleafNode*>& query32_p();
		void query32_p();
		void query64_p();
		//void query32_a();
		//void reformatIdx32_p();

		void queryFastq(std::vector<std::string>&);
		void queryFastq_p(std::vector<std::string>&);

		void queryAllFastq();

		void getRC(uint8_t*, uint8_t*, size_t);

		//void runILP_a(int, int, uint32_t, double);
		void runILP_p(int, int, uint32_t, double, double, double, double);
		
		static int symbolIdx[256];

		static int rcIdx[128];
};

#endif
