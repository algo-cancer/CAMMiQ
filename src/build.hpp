#ifndef BUILD_HPP 
#define BUILD_HPP

#include <vector> 
#include <string>
#include <unordered_map>
#include <pthread.h>
#include "hashtrie.hpp"
#include "gsa.hpp"

class FastaReader {
	private:
		/* Maximum and minimum unique substring length. */
		uint32_t L = 100;
		uint32_t K_ = 30;
		int minuL = 20;
		uint32_t maxuL = 50;
		int num_threads = 1;
		int current_tid = 0;
		pthread_mutex_t thread_lock;
		static int mode_;

		/* Current genome sequence length & num genomes & num contigs. */
		uint64_t N_ = 0;
		uint32_t M_ = 0;
		uint32_t C_ = 0;
		int max_rid = 0;

		/* The input reference sequences. */
		std::vector<uint8_t*> seqs;
		uint8_t* seqs_ = NULL;
		std::string rc_contig = "";
		size_t cur = 0, pos = 0;
		std::vector<uint64_t> contig_pos;
		std::vector<uint64_t> ref_pos;
		std::vector<uint32_t> refID;
		std::vector<bool> exist_unique_;
		uint16_t *mu_index = NULL;

		/* */
		SuffixArray* sa = NULL;
		uint8_t *occ = NULL;
		uint8_t *occ2 = NULL;
		uint32_t *GSA2 = NULL;
		uint16_t *GSA2_ = NULL;
		std::unordered_map<std::string, uint32_t> filenames;

		/* */
		pthread_spinlock_t* hasht_access;
		uint32_t HASH_LEN_ = 12;
		Hash *hasht = NULL;

	public: 
		const size_t block_size_ = 2097152;
		const uint8_t base_offset = 165;

		/* Statistics. */
		std::vector<uint64_t> minusl;
		std::vector<uint64_t> maxusl;
		std::vector<double> sumusl;
		std::vector<uint32_t> uscount;
		std::vector<uint32_t> uLmcount;
		std::vector<uint32_t> uLmcount_clark; 

		/* Constructors. */
		FastaReader(int, int, int, int);
		FastaReader(int, int, int, int, std::vector<std::string>&); 
		~FastaReader() {
			for (auto seq : seqs)
				if (seq != NULL)
					delete []seq;
			if (seqs_ != NULL)
				delete [] seqs_;
			if (sa != NULL)
				delete sa;
			if (hasht != NULL)
				delete hasht;
			if (hasht_access != NULL)
				delete hasht_access;
		}

		/* Read in fasta files. */
		void allocSeq();
		bool insertBase(std::string&);
		bool insertBaseRC(std::string&);
		bool insertContig();
		bool insertRC();
		void prepFasta(std::string&);
		void readAllFasta();
		void readFasta(std::string&);
		void readFnMap(std::string&, std::string&);

		/* Compute minimum unique substrings. */
		void* computeIndexmin();
		void* computeIndexmin_d();
		void* computeIndexmin_d_();
		static void *computeIndex_t(void *obj) {
			switch (mode_) {
				case 1:
					return ((FastaReader*) obj)->computeIndexmin();
				case 2:
				case 3:
					return ((FastaReader*) obj)->computeIndexmin_d();
				default:
					return ((FastaReader*) obj)->computeIndexmin_d_();
			}
		}

		void allocSuffixArray(int);

		void getMinUniqueSubstring();

		/* Build Index. */
		void computeIndex(int);
		//void analIndex(int, int);
		//void outputRefLength();
		void insert32(uint64_t, uint32_t, uint32_t, uint8_t);
		void insert32_d(uint64_t, uint32_t, uint32_t, uint32_t, uint8_t, uint8_t);
		//void try_increase_cnt_32(uint64_t, uint32_t, uint32_t);
		void insert64(uint64_t, uint32_t, uint32_t, uint8_t);
		void insert64_d(uint64_t, uint32_t, uint32_t, uint32_t, uint8_t, uint8_t);

		void setHashLength(uint32_t);

		//void try_printing_index_32(uint64_t, uint32_t);
}; 

#endif 
