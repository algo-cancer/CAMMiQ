#ifndef GSA_HPP_ 
#define GSA_HPP_

#include <vector>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include "divsufsort.h"

class SuffixArray {
	private:
		int64_t *SA = NULL;
		uint64_t *REV = NULL;

		uint32_t *GSA = NULL;
		uint8_t *LCP = NULL;
		uint8_t *LCP0 = NULL;
		
	public:
		uint8_t *occ = NULL;
		uint8_t *occ2 = NULL;
		uint32_t *GSA2 = NULL;
		
		/* Constructors. */
		SuffixArray(uint64_t);
		~SuffixArray() {
			if (SA != NULL)
				delete []SA;
			if (REV != NULL)
				delete []REV;
			if (GSA != NULL)
				delete []GSA;
			if (LCP != NULL)
				delete []LCP;
			if (LCP0 != NULL)
				delete []LCP0;
		}

		/* Compute the suffix array SA of s = s1$s2$...$sm;
		   The suffix array of the string [s, s + n) is stored into SA[0, n]. */
		void computeSuffixArray(uint8_t*, uint64_t, bool, bool);

		/* Computes the reverse SA index.  REVSA[i] gives the index of the suffix
		   starting a i in the SA array. */
 		void computeRevSuffixArray(uint64_t, bool);

		/* Computes the generalized SA index.
		   The id (1 <= id <= m) of each suffix is stroed into GSA[0, n]. */
		void computeGnrSuffixArray(std::vector<uint64_t>&, std::vector<uint32_t>&, uint64_t, bool);

		/* Computes the longest common prefix between adjacent suffixes. */
		void computeLcpArray(uint8_t*, uint64_t, bool, bool);

		/**/
		void computeGnrLcpArray(uint64_t, bool);
		void computeGnrLcpArray_d(uint64_t, bool);
	
		/**/
		void computeGnrLcpArray_ext(uint64_t, uint8_t, bool);
		void computeGnrLcpArray_dext(uint64_t, uint8_t, uint8_t, bool);

		/**/
		void computeMinUnique(uint64_t, bool);

		uint64_t* run(uint8_t*, std::vector<uint64_t>&, std::vector<uint32_t>&, uint64_t, int, uint8_t, uint8_t, bool, bool, bool);

		//uint64_t* run_ext_d(uint8_t*, std::vector<uint64_t>&, std::vector<uint32_t>&, uint64_t, int);

		//std::pair<uint64_t*, uint64_t*> run_both(uint8_t*, std::vector<uint64_t>&, std::vector<uint32_t>&, uint64_t);
};

#endif
