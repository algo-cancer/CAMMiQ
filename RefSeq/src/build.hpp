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
		int minL = 20;
		int num_threads = 1;
		int current_tid = 0;
		pthread_mutex_t thread_lock;
		static int mode_;

		/* Current genome sequence length & num genomes & num contigs. */
		uint64_t N_ = 0;
		uint32_t M_ = 0;
		uint32_t C_ = 0;

		/* The input reference sequences. */
		std::vector<uint8_t*> seqs;
		uint8_t* seqs_ = NULL;
		size_t cur = 0, pos = 0;
		std::vector<uint64_t> contig_pos;
		std::vector<uint64_t> ref_pos;
		std::vector<uint32_t> refID;
		std::vector<bool> exist_unique_;
		uint64_t *mu_index = NULL;

		/* */
		SuffixArray* sa = NULL;
		std::unordered_map<std::string, uint32_t> filenames;

		/* */
		int num_locks = 256;
		std::vector<pthread_spinlock_t*> hasht_access;
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
		FastaReader(int, int, int);
		FastaReader(int, int, int, std::vector<std::string>&); 
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
			for (auto spinlock : hasht_access)
				if (spinlock != NULL)
					delete spinlock;
		}

		/* Read in fasta files. */
		void allocSeq();
		bool insertBase(std::string&);
		bool insertContig();
		void prepFasta(std::string&);
		void readAllFasta();
		void readFasta(std::string&);
		void readFnMap(std::string&, std::string&);

		/* Compute minimum unique substrings. */
		void* computeIndexstats();
		void* computeIndexno();
		void* computeIndexmin();
		void* computeIndexmin_d();
		void* analIndexCLARK();
		static void *computeIndex_t(void *obj) {
			switch (mode_) {
				case 0:
					return ((FastaReader*) obj)->computeIndexstats();
				case 1:
					return ((FastaReader*) obj)->computeIndexno();
				case 2:
					return ((FastaReader*) obj)->computeIndexmin();
				case 4:
					return ((FastaReader*) obj)->computeIndexmin_d();
				default:
					return ((FastaReader*) obj)->analIndexCLARK();
			}
		}
		void allocSuffixArray();

		void getMinUniqueSubstring();

		/* Build Index. */
		void computeIndex(int);
		//void analIndex(int, int);
		//void outputRefLength();
		void insert32(uint64_t, uint32_t, uint32_t);
		void insert32_d(uint64_t, uint32_t, uint32_t);
		void insert64(uint64_t, uint32_t, uint32_t);

		void setNumLocks(int);
		void setHashLength(uint32_t);
}; 

#endif 
