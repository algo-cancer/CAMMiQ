
#include <cstring>
#include <chrono>
#include "divsufsort.h"
#include "gsa.hpp"

SuffixArray::SuffixArray(uint64_t n) {
	SA = new int64_t[n + 2];
	REV = new uint64_t[n + 2];
	GSA = new uint32_t[n + 2];
	LCP = new uint16_t[n + 2];
	GLCP = new uint16_t[n + 2];
}

void SuffixArray::computeSuffixArray(uint8_t* s, uint64_t n) {
	auto start = std::chrono::high_resolution_clock::now();
	if (divsufsort(s, SA, n) != 0) {
		fprintf(stderr, "Error occurred in computing suffix array.\n");
		abort();
	} else {

		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing suffix array: %lu ms.\n", duration);
	}
	if (debug)
		if (sufcheck(s, SA, n, true) != 0) {
			fprintf(stderr, "Sufcheck failed.\n");
			abort();
		}
}

void SuffixArray::computeRevSuffixArray(uint64_t n) {
	auto start = std::chrono::high_resolution_clock::now();
	memset(REV, 0, sizeof(uint64_t) * (n + 2));
	#pragma omp parallel for
	for (uint64_t i = 0; i < n; i++) {
		uint64_t SA_i = SA[i];
		if (REV[SA_i] == 0) {
			REV[SA_i] = i; //if (SA_i == 8788480807 || SA_i == 2506818793) fprintf(stderr, "%lu\n", i);
		}
		else {
			fprintf(stderr, "i: %lu; SA[i]: %lu; REV: %lu.\n", i, SA_i, REV[SA_i]);
			fprintf(stderr, "Error in suffix array.\n");
			abort();
		}
	}
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - start).count();
	fprintf(stderr, "Time for computing reverse suffix array: %lu ms.\n", duration);
}

void SuffixArray::computeGnrSuffixArray(std::vector<uint64_t> &spos, 
					std::vector<uint32_t> &sid, uint64_t n) {
	auto start = std::chrono::high_resolution_clock::now();
	//memset(GSA, 0, sizeof(uint32_t) * (n + 2));
	//#pragma omp parallel for
	for (uint64_t i = 0, j = 0; i < n; i++, j += (i >= spos[j]))
		GSA[REV[i]] = sid[j];
	GSA[n] = 0;
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - start).count();
	fprintf(stderr, "Time for computing generalized suffix array: %lu ms.\n", duration);
}

void SuffixArray::computeLcpArray(uint8_t* s, uint64_t n) {
	/*for (uint64_t i = 10080737490; i < 10080737500; i++) {
		uint64_t SA_i = SA[i];
		for (int j = 0; j < 80; j++)
			fprintf(stderr, "%c", (s + SA_i)[j] - 165);
		fprintf(stderr, "\t%u\n", GSA[i]);
	}
	for (uint64_t i = 9147092977; i < 9147092987; i++) {
		uint64_t SA_i = SA[i];
		for (int j = 0; j < 80; j++)
			fprintf(stderr, "%c", (s + SA_i)[j] - 165);
		fprintf(stderr, "\t%u\n", GSA[i]);
	}*/
		/*for (uint64_t i = 11268149250; i < 11268149290; i++) {
			uint64_t SA_i = SA[i];
			fprintf(stderr, "i: %lu\t", i);
			for (int j = 0; j < 50; j++)
				fprintf(stderr, "%c", (s + SA_i)[j] - 165);
			fprintf(stderr, "\t%u\n", GSA[i]);
		}*/

	auto start = std::chrono::high_resolution_clock::now();
	uint64_t totalLCP = 0;
	double avgLCP = 0.0;
	uint64_t len = 0;
	//memset(LCP, 0, sizeof(uint16_t) * n);
	//#pragma omp parallel firstprivate(len)
	for (uint64_t i = 0; i < n; i++) {
		//k: 0 - n - 1; SA: 0 - n - 1
		uint64_t k = REV[i];
		if (k == 0)
			continue;
		uint64_t j = SA[k - 1];
		for(; i + len < n && j + len < n && s[i + len] == s[j + len]; len++);
		totalLCP += len;
		LCP[k] = (len >= 0xFFFFull) ? UINT16_MAX : (uint16_t)(len & 0xFFFF);
		if (len > 0) len--;
	}
	LCP[n] = 0;
	/*for (uint64_t i = 10080737490; i < 10080737500; i++) {
		fprintf(stderr, "i: %lu; LCP[i]: %u\n", i, LCP[i]);
	}*/

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - start).count();
	fprintf(stderr, "Time for computing LCP array: %lu ms.\n", duration);
	avgLCP = 1.0 * totalLCP / n;
	fprintf(stderr, "The average LCP is %lf.\n", avgLCP);
}

void SuffixArray::computeGnrLcpArray(uint64_t n) {
	auto start = std::chrono::high_resolution_clock::now();
	uint16_t minlcp = UINT16_MAX;
	uint64_t nextd = 0;
	memset(GLCP, 0, sizeof(uint16_t) * (n + 2));
	for (uint64_t i = 0; i < n - 1; i += (nextd + 1), minlcp = UINT16_MAX, nextd = 0) {
		for (; GSA[i + nextd] == GSA[i + nextd + 1]; nextd++);
		for (int64_t j = nextd; j >= 0; j--) {
			minlcp = std::min(minlcp, LCP[i + j + 1]);
			GLCP[i + j] = minlcp;
		}
		//if (i % 10000000 == 0)
		//	fprintf(stderr, "%lu.\n", i);
	}
	GLCP[n - 1] = 0;
	minlcp = UINT16_MAX;
	nextd = 0;
	for (int64_t i = n - 1; i >= 1; i -= (nextd + 1), minlcp = UINT16_MAX, nextd = 0) {
		for (; GSA[i - nextd] == GSA[i - nextd - 1]; nextd++);
		for (int64_t j = nextd; j >= 0; j--) {
			minlcp = std::min(minlcp, LCP[i - j]);
			GLCP[i - j] = std::max(GLCP[i - j], minlcp);
		}
		//if (i % 10000000 == 0)
		//	fprintf(stderr, "%lu.\n", i);
	}

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - start).count();
	fprintf(stderr, "Time for computing generalized LCP array: %lu ms.\n", duration);
}

void SuffixArray::computeGnr2LcpArray(uint64_t n) {
	auto start = std::chrono::high_resolution_clock::now();
	uint16_t minlcp = UINT16_MAX;
	uint64_t nextd = 0, begin = 0, end = n - 1;
	for (; GSA[end] == GSA[end - 1]; end--);
	//memset(GLCP, 255, sizeof(uint16_t) * (n + 2));
	for (uint64_t i = begin; i < end; i += (nextd + 1), minlcp = UINT16_MAX, nextd = 0) {
		for (; GSA[i + nextd] == GSA[i + nextd + 1]; nextd++);
		for (int64_t j = nextd; j >= 0; j--) {
			minlcp = std::min(minlcp, LCP[i + j + 1]);
			GLCP[i + j] = minlcp;
		}
	}
	for (uint64_t i = end; i < n; i++)
		GLCP[i] = 0;
	minlcp = UINT16_MAX;
	nextd = 0;
	//for (begin = n - 1; GSA[begin] == GSA[begin - 1]; begin--);
	//for (end = 0; GSA[end] == GSA[end + 1]; end++);
	begin = n - 1;
	end = 1;
	for (int64_t i = begin; i > 1; i -= (nextd + 1), minlcp = UINT16_MAX, nextd = 0) {
		for (; GSA[i - nextd] == GSA[i - nextd - 1]; nextd++);
		for (int64_t j = nextd; j >= 0; j--) {
			minlcp = std::min(minlcp, LCP[i - j]);
			if (GLCP[i - j] < minlcp) {
				uint16_t min2lcp = UINT16_MAX;
				uint64_t i_ = i - nextd - 1;
				for (; GSA[i_] == GSA[i_ - 1]; i_--)	
					 min2lcp = std::min(min2lcp, LCP[i_]);
				min2lcp = std::min(min2lcp, LCP[i_]);
				min2lcp = std::min(min2lcp, minlcp);
				GLCP[i - j] = std::max(GLCP[i - j], min2lcp);
			}
			if (GLCP[i - j] > minlcp) {
				uint16_t min2lcp = UINT16_MAX;
				uint64_t i_ = i;
				for (; GSA[i_] == GSA[i_ + 1]; i_++)
					min2lcp = std::min(min2lcp, LCP[i_ + 1]);
				min2lcp = std::min(min2lcp, LCP[i_ + 1]);
				for (i_ = i_ + 1; GSA[i_] == GSA[i_ + 1]; i_++)
					min2lcp = std::min(min2lcp, LCP[i_ + 1]);
				min2lcp = std::min(min2lcp, LCP[i_ + 1]);
				GLCP[i - j] = std::max(minlcp, min2lcp);
			}
				
			//GLCP[i - j] = std::min(GLCP[i - j], minlcp);
		}
	}
	/*for (uint64_t i = 10080737490; i < 10080737500; i++) {
		fprintf(stderr, "i: %lu; GLCP[i]: %u\n", i, GLCP[i]);
	}*/
	

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - start).count();
	fprintf(stderr, "Time for computing generalized LCP array: %lu ms.\n", duration);
}

void SuffixArray::computeMinUnique(uint64_t n) {
	auto start = std::chrono::high_resolution_clock::now();
	uint64_t *MU = REV;
	memset(MU, 0, sizeof(uint64_t) * (n + 2));
	//#pragma omp parallel
	for (uint64_t i = 0; i < n; i++) {
		uint64_t SA_i = SA[i];
		MU[SA_i + GLCP[i] + 1] = std::max(MU[SA_i + GLCP[i] + 1], SA_i + 1);
	}

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - start).count();
	fprintf(stderr, "Time for computing minimum unique substrings: %lu ms.\n", duration);
}

uint64_t* SuffixArray::run(uint8_t* s, std::vector<uint64_t> &spos, 
				std::vector<uint32_t> &sid, uint64_t n) {
	fprintf(stderr, "Bgein computing suffix array.\n");
	computeSuffixArray(s, n);
	fprintf(stderr, "Begin computing reverse suffix array.\n");
	computeRevSuffixArray(n);
	fprintf(stderr, "Begin computing generalized suffix array.\n");
	computeGnrSuffixArray(spos, sid, n);
	fprintf(stderr, "Begin computing LCP array.\n");
	computeLcpArray(s, n);
	fprintf(stderr, "Begin computing generalized LCP array.\n");
	//computeGnrLcpArray(n);
	computeGnr2LcpArray(n);
	
	/*for (uint64_t i = 0; i < n; i++) {
		if (GLCP[i] <= 4) {
			bool special_character_ = false;
			for (int j = 0; j < 5; j++)
				if (((s + SA[i])[j] - 165 != 65) && 
			   	 ((s + SA[i])[j] - 165 != 67) &&
			    	((s + SA[i])[j] - 165 != 71) &&
			    	((s + SA[i])[j] - 165 != 84)) {
					special_character_ = true;
					break;
				}
			if (!special_character_) {
				fprintf(stderr, "i: %lu; SA[i]: %lu; GLCP[i]: %u\t", i, SA[i], GLCP[i]);
				for (int j = 0; j < 5; j++)
					fprintf(stderr, "%c", (s + SA[i])[j] - 165);
				fprintf(stderr, "\n");
			}
		}
	}*/
	fprintf(stderr, "Begin computing minimum unique substrings.\n");
	computeMinUnique(n);
	return REV;
}
