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

		// Modified 10/28
		uint32_t ext_len = 3;

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
		
		//new: keep the positions
		//uint32_t ci_inserted = 0;
		//uint64_t last_inserted_i = 0;
		//uint64_t last_inserted_l = 0;
		//std::unordered_map<uint32_t, uint32_t> uc_cnts;

		/* The input reference sequences. */
		std::vector<uint8_t*> seqs;
		uint8_t* seqs_ = NULL;
		std::string rc_contig = "";
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

		//std::vector<uint32_t> least_50_sp = {573062, 1898044, 399726, 768493, 1476577, 547144, 550542, 547145, 2105, 1914410, 544556, 1703967, 1423, 1891098, 2016519, 2572089, 1703962, 1703965, 1962264, 29323, 28152, 2567883, 2567881, 1703964, 1837130, 1949203, 1967781, 2565557, 1806506, 79881, 1570328, 475375, 547146, 1703961, 1148, 496866, 547143, 1795, 2565560, 2565561, 2116702, 768490, 83617, 1703968, 1490052, 2144, 1703966, 1468409, 1740162, 2016518}; 

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
		bool insertBaseRC(std::string&);
		bool insertContig();
		bool insertRC();
		void prepFasta(std::string&);
		void readAllFasta();
		void readFasta(std::string&);
		void readFnMap(std::string&, std::string&);

		/* Compute minimum unique substrings. */
		void* computeIndexstats();
		void* computeIndexdense();
		void* computeIndexdense_d();
		void* computeIndexmin();
		//void* computeIndexsum();
		//void* computeIndexmin_ext();
		void* computeIndexmin_d();
		void* analIndexCLARK();
		static void *computeIndex_t(void *obj) {
			switch (mode_) {
				case 0:
					return ((FastaReader*) obj)->computeIndexstats();
				case 1:
					return ((FastaReader*) obj)->computeIndexdense();
				case 11:
					return ((FastaReader*) obj)->computeIndexdense_d();
				case 2:
					return ((FastaReader*) obj)->computeIndexmin();
				case 4:
					return ((FastaReader*) obj)->computeIndexmin_d();
				case 5:
					//return ((FastaReader*) obj)->computeIndexsum();
				default:
					return ((FastaReader*) obj)->analIndexCLARK();
			}
		}

		void allocSuffixArray(bool);

		void getMinUniqueSubstring();

		/* Build Index. */
		void computeIndex(int);
		//void analIndex(int, int);
		//void outputRefLength();
		void insert32(uint64_t, uint32_t, uint32_t, uint8_t*);
		void insert32_d(uint64_t, uint32_t, uint32_t, uint32_t*, uint8_t*, uint8_t*);
		//void try_increase_cnt_32(uint64_t, uint32_t, uint32_t);
		void insert64(uint64_t, uint32_t, uint32_t, uint8_t*);
		void insert64_d(uint64_t, uint32_t, uint32_t, uint32_t*, uint8_t*, uint8_t*);

		void setNumLocks(int);
		void setHashLength(uint32_t);

		//void try_printing_index_32(uint64_t, uint32_t);
}; 

#endif 
