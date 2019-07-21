#ifndef GSA_HPP_ 
#define GSA_HPP_

#include <vector>
//#include <ctime>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include "divsufsort.h"

class SuffixArray {
	private:
		int64_t *SA = NULL;
		uint64_t *REV = NULL;

		uint32_t *GSA = NULL;
		uint16_t *LCP = NULL;
		uint16_t *GLCP = NULL;

	public:
		//std::clock_t start;
		//double duration;
		bool debug = false;
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
			if (GLCP != NULL)
				delete []GLCP;
		}

		/* Compute the suffix array SA of s = s1$s2$...$sm;
		   The suffix array of the string [s, s + n) is stored into SA[0, n]. */
		void computeSuffixArray(uint8_t*, uint64_t);

		/* Computes the reverse SA index.  REVSA[i] gives the index of the suffix
		   starting a i in the SA array. */
 		void computeRevSuffixArray(uint64_t);

		/* Computes the generalized SA index.
		   The id (1 <= id <= m) of each suffix is stroed into GSA[0, n]. */
		void computeGnrSuffixArray(std::vector<uint64_t>&, std::vector<uint32_t>&, uint64_t);

		/* Computes the longest common prefix between adjacent suffixes. */
		void computeLcpArray(uint8_t*, uint64_t);

		/**/
		void computeGnrLcpArray(uint64_t);
		void computeGnr2LcpArray(uint64_t);

		/**/
		void computeMinUnique(uint64_t);

		uint64_t* run(uint8_t*, std::vector<uint64_t>&, std::vector<uint32_t>&, uint64_t);

};

#endif
