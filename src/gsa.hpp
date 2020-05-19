#ifndef GSA_HPP_ 
#define GSA_HPP_

#include <vector>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include "divsufsort.h"

class SuffixArray {
	private:
		/* Buffer of size n + 2. */
		uint64_t *buffer = NULL;

		/* Num bits for REV and genome ID. */
		int gsa_unit_ = 16;

		/* Generalized enhanced suffix array. */
		int64_t *SA = NULL;
		uint64_t *REV = NULL;
		uint32_t *GSA32 = NULL;
		uint16_t *GSA16 = NULL;
		uint16_t *LCP = NULL;

		/* Longest common prefix array to compute. */
		uint16_t *LCP0 = NULL;
		bool lcp_prepared__ = 0;
		
	public:
		uint8_t *occ = NULL;
		uint8_t *occ2 = NULL;
		uint32_t *GSA32_ = NULL;
		uint16_t *GSA16_ = NULL;
		
		/* Constructors. */
		SuffixArray(uint64_t, int);
		~SuffixArray() {
			if (SA != NULL)
				delete []SA;
			if (REV != NULL)
				delete []REV;
			if (buffer != NULL)
				delete []buffer;
			if (gsa_unit_ == 32 && GSA32_ != NULL)
				delete []GSA32_;
		}
		void allocBuffer(uint64_t);

		/* Compute the suffix array SA of s = s1$s2$...$sm;
		   The suffix array of the string [s, s + n) is stored into SA[0, n]. */
		void computeSuffixArray(uint8_t*, uint64_t, bool, bool);

		/* Compute the reverse SA index.  REVSA[i] gives the index of the suffix
		   starting a i in the SA array. */
 		void computeRevSuffixArray(uint64_t, bool);

		/* Compute the generalized SA index.
		   The id (1 <= id <= m) of each suffix is stored into GSA[0, n]. */
		void fillGnrSuffixArray16(uint16_t*, uint64_t*, std::vector<uint64_t>&, 
						std::vector<uint32_t> &, uint64_t, uint64_t);
		void fillGnrSuffixArray32(uint32_t*, uint64_t*, std::vector<uint64_t>&, 
						std::vector<uint32_t> &, uint64_t, uint64_t);
		void computeGnrSuffixArray16(std::vector<uint64_t>&, std::vector<uint32_t>&, 
						uint64_t, bool, bool);
		void computeGnrSuffixArray32(std::vector<uint64_t>&, std::vector<uint32_t>&, 
						uint64_t, bool, bool);

		/* Compute the longest common prefix between adjacent suffices. */
		void computeAvgLcp(uint8_t*, uint64_t, bool);
		void computeLcpArray(uint8_t*, uint64_t, bool);
		void prepareGnrLCP(uint64_t, bool);
		void reloadLCP(uint64_t, bool);
		void computeGnrLcpArray16(uint64_t, uint16_t, bool);
		void computeGnrLcpArray32(uint64_t, uint16_t, bool);
		void computeGnrLcpArray16_d(uint64_t, uint16_t, uint16_t, bool);
		void computeGnrLcpArray32_d(uint64_t, uint16_t, uint16_t, bool);
		void computeOCC16(uint64_t, bool);
		void computeOCC32(uint64_t, bool);
		void computeOCC16_d(uint64_t, uint16_t, bool);
		void computeOCC32_d(uint64_t, uint16_t, bool);

		/* Compute shortest (doubly-)unique substrings. 
		   The MU array is stored in LCP[0, n] for memory reuse.*/
		void computeMinUnique(uint64_t, bool);
		void computeMinUnique(uint64_t, uint16_t, bool);

		/* IO helpers. */
		void writeArray16(uint16_t*, uint64_t, std::string&);
		void writeArray32(uint32_t*, uint64_t, std::string&);
		void writeArray64(uint64_t*, uint64_t, std::string&);
		void readArray16(uint16_t*, uint64_t, std::string&);
		void readArray32(uint32_t*, uint64_t, std::string&);
		void readArray64(uint64_t*, uint64_t, std::string&);
		
		/* Interface. */
		uint16_t* run(uint8_t*, std::vector<uint64_t>&, std::vector<uint32_t>&, 
				uint64_t, int, uint16_t, uint16_t, bool, bool, bool);
};

#endif
