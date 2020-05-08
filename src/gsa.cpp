#include <cstring>
#include <chrono>
#include "divsufsort.h"
#include "gsa.hpp"

SuffixArray::SuffixArray(uint64_t n) {
	SA = new int64_t[n + 2];
	REV = new uint64_t[n + 2];
	GSA = new uint32_t[n + 2];
	LCP = new uint8_t[n + 2];
	LCP0 = new uint8_t[n + 2];
	GSA2 = new uint32_t[n + 2];
	occ = new uint8_t[n + 2];
	occ2 = new uint8_t[n + 2];
}

void SuffixArray::computeSuffixArray(uint8_t* s, uint64_t n, bool debug, bool debug_sa) {
	auto start = std::chrono::high_resolution_clock::now();
	if (divsufsort(s, SA, n) != 0) {
		fprintf(stderr, "Error occurred in computing suffix array.\n");
		abort();
	} else {
		if (debug) {
			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
					(std::chrono::high_resolution_clock::now() - start).count();
			fprintf(stderr, "Time for computing suffix array: %lu ms.\n", duration);
		}
	}
	if (debug_sa)
		if (sufcheck(s, SA, n, true) != 0) {
			fprintf(stderr, "Sufcheck failed.\n");
			abort();
		}
}

void SuffixArray::computeRevSuffixArray(uint64_t n, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	memset(REV, 0, sizeof(uint64_t) * (n + 2));
	#pragma omp parallel for
	for (uint64_t i = 0; i < n; i++) {
		uint64_t SA_i = SA[i];
		if (REV[SA_i] == 0) {
			REV[SA_i] = i;
		}
		else {
			fprintf(stderr, "i: %lu; SA[i]: %lu; REV: %lu.\n", i, SA_i, REV[SA_i]);
			fprintf(stderr, "Error in suffix array.\n");
			abort();
		}
	}
	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing reverse suffix array: %lu ms.\n", duration);
	}
}

void SuffixArray::computeGnrSuffixArray(std::vector<uint64_t> &spos, 
					std::vector<uint32_t> &sid, uint64_t n, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	//memset(GSA, 0, sizeof(uint32_t) * (n + 2));
	uint64_t j = 0;
	#pragma omp parallel for firstprivate(j)
	for (uint64_t i = 0; i < n; i++) {
		while (j < spos.size() && i >= spos[j])
			j++;
		GSA[REV[i]] = sid[j];
	}
	GSA[n] = 0;
	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing generalized suffix array: %lu ms.\n", duration);
	}
}

void SuffixArray::computeLcpArray(uint8_t* s, uint64_t n, bool debug, bool print_avg) {
	auto start = std::chrono::high_resolution_clock::now();
	uint64_t totalLCP = 0;
	double avgLCP = 0.0;
	uint64_t len = 0;
	//memset(LCP, 0, sizeof(uint16_t) * n);
	if (print_avg == 0) {
		#pragma omp parallel for firstprivate(len)
		for (uint64_t i = 0; i < n; i++) {
			//k: 0 - n - 1; SA: 0 - n - 1
			uint64_t k = REV[i];
			if (k == 0)
				continue;
			uint64_t j = SA[k - 1];
			for(; i + len < n && j + len < n && s[i + len] == s[j + len]; len++);
			LCP[k] = (len >= 0xFFull) ? UINT8_MAX : (uint8_t)(len & 0xFF);
			if (len > 0) len--;
		}
	} else {
		for (uint64_t i = 0; i < n; i++) {
			//k: 0 - n - 1; SA: 0 - n - 1
			uint64_t k = REV[i];
			if (k == 0)
				continue;
			uint64_t j = SA[k - 1];
			for(; i + len < n && j + len < n && s[i + len] == s[j + len]; len++);
			totalLCP += len;
			LCP[k] = (len >= 0xFFull) ? UINT8_MAX : (uint8_t)(len & 0xFF);
			if (len > 0) len--;
		}
	}
	LCP[n] = 0;
	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing LCP array: %lu ms.\n", duration);
	}
	if (print_avg) {
		avgLCP = 1.0 * totalLCP / n;
		fprintf(stderr, "The average LCP is %lf.\n", avgLCP);
	}
}

void SuffixArray::computeGnrLcpArray(uint64_t n, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	//memset(LCP0, 0, sizeof(uint8_t) * (n + 2));
	uint8_t minlcp = UINT8_MAX;
	uint64_t nextd = 0, begin = 0, end = n - 1;
	for (; GSA[end] == GSA[end - 1]; end--);
	for (uint64_t i = begin; i < end; i += (nextd + 1), minlcp = UINT8_MAX, nextd = 0) {
		for (; GSA[i + nextd] == GSA[i + nextd + 1]; nextd++);
		for (int64_t j = nextd; j >= 0; j--) {
			minlcp = std::min(minlcp, LCP[i + j + 1]);
			LCP0[i + j] = minlcp;
		}
	}
	for (uint64_t i = end; i < n; i++)
		LCP0[i] = 0;
	minlcp = UINT8_MAX;
	nextd = 0;
	for (end = 0; GSA[end] == GSA[end + 1]; end++);
	begin = n - 1;
	for (uint64_t i = begin; i > end; i -= (nextd + 1), minlcp = UINT8_MAX, nextd = 0) {
		for (; GSA[i - nextd] == GSA[i - nextd - 1]; nextd++);
		for (int64_t j = nextd; j >= 0; j--) {
			minlcp = std::min(minlcp, LCP[i - j]);
			LCP0[i - j] = std::max(LCP0[i - j], minlcp);
		}
	}
	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing LCP0 array: %lu ms.\n", duration);
	}
}

void SuffixArray::computeGnrLcpArray_ext(uint64_t n, uint8_t el, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	uint8_t minlcp = UINT8_MAX;
	uint64_t nextd = 0, begin = 0, end = n - 1;
	for (; GSA[end] == GSA[end - 1]; end--);

	//#pragma omp parallel for
	for (uint64_t i = begin; i < end; i += (nextd + 1), minlcp = UINT8_MAX, nextd = 0) {
		for (; GSA[i + nextd] == GSA[i + nextd + 1]; nextd++);
		for (int64_t j = nextd; j >= 0; j--) {
			minlcp = std::min(minlcp, LCP[i + j + 1]);
			LCP0[i + j] = std::max(el, minlcp);
		}
	}
	
	for (uint64_t i = end; i < n; i++)
		LCP0[i] = 0;
	minlcp = UINT8_MAX;
	nextd = 0;
	for (end = 0; GSA[end] == GSA[end + 1]; end++);
	begin = n - 1;
	
	//#pragma omp parallel for
	for (uint64_t i = begin; i > end; i -= (nextd + 1), minlcp = UINT8_MAX, nextd = 0) {
		for (; GSA[i - nextd] == GSA[i - nextd - 1]; nextd++);
		for (int64_t j = nextd; j >= 0; j--) {
			minlcp = std::min(minlcp, LCP[i - j]);
			LCP0[i - j] = std::max(LCP0[i - j], minlcp);
		}
	} 

	memset(occ, 1, sizeof(uint8_t) * (n + 2));
	for (uint64_t i = begin; i > end; i -= (nextd + 1)) {
		for (nextd = 0; GSA[i - nextd] == GSA[i - nextd - 1]; nextd++);
		//occ[i] = 1;
		if (nextd > 0) {
			for (uint64_t j = 0; j <= nextd; j++) {
				minlcp = UINT8_MAX;
				for (int64_t j_ = j - 1; j_ >= 0; j_--) {
					minlcp = std::min(minlcp, LCP[i - nextd + j_ + 1]);
					if (minlcp > LCP0[i - nextd + j]) 
						occ[SA[i - nextd + j]]++;
				}
				minlcp = UINT8_MAX;
				for (uint64_t j_ = j + 1; j_ <= nextd; j_++) {
					minlcp = std::min(minlcp, LCP[i - nextd + j_]);
					if (minlcp > LCP0[i - nextd + j])
						occ[SA[i - nextd + j]]++;
				}
			}
		}
	} 

	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing generalized LCP array: %lu ms.\n", duration);
	}
}

void SuffixArray::computeGnrLcpArray_d(uint64_t n, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	uint8_t minlcp = UINT8_MAX;
	uint64_t nextd = 0, begin = 0, end = n - 1;
	for (; GSA[end] == GSA[end - 1]; end--);

	for (uint64_t i = begin; i < end; i += (nextd + 1), minlcp = UINT8_MAX, nextd = 0) {
		for (; GSA[i + nextd] == GSA[i + nextd + 1]; nextd++);
		for (int64_t j = nextd; j >= 0; j--) {
			minlcp = std::min(minlcp, LCP[i + j + 1]);
			LCP0[i + j] = minlcp;
		}
	}
	for (uint64_t i = end; i < n; i++)
		LCP0[i] = 0;
	minlcp = UINT8_MAX;
	nextd = 0;
	begin = n - 1;
	for (end = 0; GSA[end] == GSA[end + 1]; end++);
	for (; GSA[end] == GSA[end + 1]; end++);

	for (uint64_t i = begin; i > end; i -= (nextd + 1), minlcp = UINT8_MAX, nextd = 0) {
		for (; GSA[i - nextd] == GSA[i - nextd - 1]; nextd++);
		for (int64_t j = nextd; j >= 0; j--) {
			minlcp = std::min(minlcp, LCP[i - j]);
			if (LCP0[i - j] < minlcp) {
				uint8_t min2lcp = UINT8_MAX;
				uint64_t i_ = i - nextd - 1;
				for (; GSA[i_] == GSA[i_ - 1]; i_--)	
					 min2lcp = std::min(min2lcp, LCP[i_]);
				min2lcp = std::min(min2lcp, LCP[i_]);
				min2lcp = std::min(min2lcp, minlcp);
				LCP0[i - j] = std::max(LCP0[i - j], min2lcp);
			}
			if (LCP0[i - j] > minlcp) {
				uint8_t min2lcp = UINT8_MAX;
				uint64_t i_ = i;
				for (; GSA[i_] == GSA[i_ + 1] && i_ < n; i_++)
					min2lcp = std::min(min2lcp, LCP[i_ + 1]);
				min2lcp = std::min(min2lcp, LCP[i_ + 1]);
				for (i_ = i_ + 1; GSA[i_] == GSA[i_ + 1] && i_ < n; i_++)
					min2lcp = std::min(min2lcp, LCP[i_ + 1]);
				min2lcp = std::min(min2lcp, LCP[i_ + 1]);
				LCP0[i - j] = std::max(minlcp, min2lcp);
			}
				
			//LCP0[i - j] = std::min(LCP0[i - j], minlcp);
		}
	}
	/*for (uint64_t i = 10080737490; i < 10080737500; i++) {
		fprintf(stderr, "i: %lu; LCP0[i]: %u\n", i, LCP0[i]);
	}*/
	
	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing generalized LCP array: %lu ms.\n", duration);
	}
}

void SuffixArray::computeGnrLcpArray_dext(uint64_t n, uint8_t el, uint8_t ulmax, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	uint8_t minlcp = UINT8_MAX;
	uint64_t nextd = 0, begin = 0, end = n - 1;
	for (; GSA[end] == GSA[end - 1]; end--);
	//fprintf(stderr, "ulmax: %u\n", ulmax);
	memset(GSA2, 0, sizeof(uint32_t) * (n + 2));

	for (uint64_t i = begin; i < end; i += (nextd + 1), minlcp = UINT8_MAX, nextd = 0) {
		for (; GSA[i + nextd] == GSA[i + nextd + 1]; nextd++);
		for (int64_t j = nextd; j >= 0; j--) {
			minlcp = std::min(minlcp, LCP[i + j + 1]);
			LCP0[i + j] = minlcp;
			GSA2[SA[i + j]] = GSA[i + nextd + 1];
		}
	}
	for (uint64_t i = end; i < n; i++)
		LCP0[i] = 0;
	minlcp = UINT8_MAX;
	nextd = 0;
	begin = n - 1;
	for (end = 0; GSA[end] == GSA[end + 1]; end++);
	for (; GSA[end] == GSA[end + 1]; end++);

	for (uint64_t i = begin; i > end; i -= (nextd + 1), minlcp = UINT8_MAX, nextd = 0) {
		for (; GSA[i - nextd] == GSA[i - nextd - 1]; nextd++);
		for (int64_t j = nextd; j >= 0; j--) {
			minlcp = std::min(minlcp, LCP[i - j]);
			if (LCP0[i - j] < minlcp) {
				uint8_t min2lcp = UINT8_MAX;
				uint64_t i_ = i - nextd - 1;
				for (; GSA[i_] == GSA[i_ - 1]; i_--)	
					 min2lcp = std::min(min2lcp, LCP[i_]);
				min2lcp = std::min(min2lcp, LCP[i_]);
				min2lcp = std::min(min2lcp, minlcp);
				//if (i > 10000000 && i < 10000100)
				//	fprintf(stderr, "LCP[%lu] = %u; minlcp = %u; min2lcp = %u\n", i - j, LCP0[i - j], minlcp, min2lcp);
				LCP0[i - j] = std::max(LCP0[i - j], min2lcp);
				LCP0[i - j] = std::max(LCP0[i - j], el);
				GSA2[SA[i - j]] = GSA[i - nextd - 1];
				if (LCP0[i - j] >= minlcp)
					LCP0[i - j] = ulmax + 2;
				//if (i > 10000000 & i < 10000100)
				//	fprintf(stderr, "LCP[%lu] = %u\n", i - j, LCP0[i - j]);
			} else {
				if (LCP0[i - j] > minlcp) {
					uint8_t min2lcp = UINT8_MAX;
					uint64_t i_ = i;
					for (; GSA[i_] == GSA[i_ + 1] && i_ < n; i_++)
						min2lcp = std::min(min2lcp, LCP[i_ + 1]);
					min2lcp = std::min(min2lcp, LCP[i_ + 1]); 
					for (i_ = i_ + 1; GSA[i_] == GSA[i_ + 1] && i_ < n; i_++)
						min2lcp = std::min(min2lcp, LCP[i_ + 1]);
					min2lcp = std::min(min2lcp, LCP[i_ + 1]);
					uint8_t lcp0 = std::max(minlcp, min2lcp);
					lcp0 = std::max(lcp0, el);
					if (lcp0 >= LCP0[i - j])
						LCP0[i - j] = ulmax + 2;
					else
						LCP0[i - j] = lcp0;
				} else
					LCP0[i - j] = ulmax + 2;
			}
		}
	}

	memset(occ, 0, sizeof(uint8_t) * (n + 2));
	memset(occ2, 0, sizeof(uint8_t) * (n + 2));
	for (uint64_t i = begin; i > end; i--) {
		if (LCP0[i] <= ulmax) {
			occ[SA[i]] = 1;
			minlcp = UINT8_MAX;
			for (uint64_t j = 0; (i - j > end) && ((GSA[i - j - 1] == GSA[i]) || (GSA[i - j - 1] == GSA2[SA[i]])); j++) {
				minlcp = std::min(minlcp, LCP[i - j]);
				if (minlcp > LCP0[i]) {
					if (GSA[i - j - 1] == GSA[i])
						occ[SA[i]]++;
					if (GSA[i - j - 1] == GSA2[SA[i]])
						occ2[SA[i]]++;
				}
			}
			minlcp = UINT8_MAX;
			for (uint64_t j = 0; (i + j <= begin) && ((GSA[i + j + 1] == GSA[i]) || (GSA[i + j + 1] == GSA2[SA[i]])); j++) {
				minlcp = std::min(minlcp, LCP[i + j + 1]);
				if (minlcp > LCP0[i]) {
					if (GSA[i + j + 1] == GSA[i])
						occ[SA[i]]++;
					if (GSA[i + j + 1] == GSA2[SA[i]])
						occ2[SA[i]]++;				
				}
			}
		}
	}
 
	//for (uint64_t i = 10000000; i <= 10000050; i++)
	//	fprintf(stderr, "i = %lu; GSA[i] = %u; LCP[i] = %u; LCP0[i] = %u; GSA2[i] = %u; occ[i] = %u; occ2[i] = %u\n", i, 
	//		GSA[i], LCP[i], LCP0[i], GSA2[SA[i]], occ[SA[i]], occ2[SA[i]]);

	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing generalized LCP array: %lu ms.\n", duration);
	}
}


void SuffixArray::computeMinUnique(uint64_t n, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	uint64_t *MU = REV;
	memset(MU, 0, sizeof(uint64_t) * (n + 2));
	#pragma omp parallel for
	for (uint64_t i = 0; i < n; i++) {
		uint64_t SA_i = SA[i];
		MU[SA_i + LCP0[i] + 1] = std::max(MU[SA_i + LCP0[i] + 1], SA_i + 1);
	}
	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing minimum unique substrings: %lu ms.\n", duration);
	}
}

uint64_t* SuffixArray::run(uint8_t* s, std::vector<uint64_t> &spos, 
				std::vector<uint32_t> &sid, uint64_t n, int mode, 
				uint8_t ulmax, uint8_t el, bool debug,
				bool debug_sa, bool print_avg) {
	switch (mode) {
		case 0:
			computeSuffixArray(s, n, debug, debug_sa);
			computeRevSuffixArray(n, debug);
			computeGnrSuffixArray(spos, sid, n, debug);
			computeLcpArray(s, n, debug, print_avg);
			return NULL;
		case 1:
			computeGnrLcpArray(n, debug);
			break;
		case 2:
			computeGnrLcpArray_d(n, debug);
			break;
		case 3:
			computeGnrLcpArray_ext(n, el, debug);
			break;
		case 4:
			computeGnrLcpArray_dext(n, el, ulmax, debug);
			break;
	}

	computeMinUnique(n, debug);
	return REV;
}
