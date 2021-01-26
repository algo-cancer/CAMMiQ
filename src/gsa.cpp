#include <cstring>
#include <fstream>
#include <chrono>
#include <omp.h>

#include "divsufsort.h"
#include "gsa.hpp"

SuffixArray::SuffixArray(uint64_t n, int gsa_unit) {
	SA = new int64_t[n + 2];
	REV = new uint64_t[n + 10];
	gsa_unit_ = gsa_unit;
}

void SuffixArray::allocBuffer(uint64_t n) {
	buffer = new uint64_t[n];
	//memset(buffer, 0, sizeof(uint64_t) * (n));
}

void SuffixArray::computeSuffixArray(uint8_t *s, uint64_t n, bool debug, bool debug_sa) {
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
		if (REV[SA_i] == 0)
			REV[SA_i] = i;
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

void SuffixArray::fillGnrSuffixArray16(uint16_t *gsa, uint64_t *rev, std::vector<uint64_t> &spos, 
					std::vector<uint32_t> &sid, uint64_t l, uint64_t r) {
	uint64_t j = 0;
	#pragma omp parallel for firstprivate(j)
	for (uint64_t i = l; i < r; i++) {
		while (j < spos.size() && i >= spos[j])
			j++;
		gsa[rev[i]] = sid[j];
	}
}

void SuffixArray::fillGnrSuffixArray32(uint32_t *gsa, uint64_t *rev, std::vector<uint64_t> &spos, 
					std::vector<uint32_t> &sid, uint64_t l, uint64_t r) {
	uint64_t j = 0;
	#pragma omp parallel for firstprivate(j)
	for (uint64_t i = l; i < r; i++) {
		while (j < spos.size() && i >= spos[j])
			j++;
		gsa[rev[i]] = sid[j];
	}
}

void SuffixArray::computeGnrSuffixArray16(std::vector<uint64_t> &spos, 
					std::vector<uint32_t> &sid, uint64_t n, 
					bool debug, bool work_in_buffer) {
	auto start = std::chrono::high_resolution_clock::now();
	
	/* GSA will be stored in SA or buffer. */
	std::string gsa_fn = "gsa.bin", sa_fn = "sa0.bin";
	if (work_in_buffer) {
		GSA16 = (uint16_t*) buffer;
		fillGnrSuffixArray16(GSA16, REV, spos, sid, 0, n);
		GSA16[n] = 0;
		writeArray16(GSA16, n + 1, gsa_fn);
	} else {
		writeArray64((uint64_t *) SA, n / 4 + 1, sa_fn);
		GSA16 = (uint16_t*) SA;	
		fillGnrSuffixArray16(GSA16, REV, spos, sid, 0, n);
		GSA16[n] = 0;
		writeArray16(GSA16, n + 1, gsa_fn);
		readArray64((uint64_t *) SA, n / 4 + 1, sa_fn);
		if (remove("sa0.bin") == 0)
			fprintf(stderr, "Deleted temp file sa0.bin.\n");
		else
			fprintf(stderr, "Error in deleting file sa0.bin.\n");
	}
	
	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing generalized suffix array: %lu ms.\n", duration);
	}
}

void SuffixArray::computeGnrSuffixArray32(std::vector<uint64_t> &spos, 
					std::vector<uint32_t> &sid, uint64_t n, 
					bool debug, bool work_in_buffer) {
	auto start = std::chrono::high_resolution_clock::now();

	/* GSA will be stored in SA or buffer. */
	std::string gsa_fn = "gsa.bin", sa_fn = "sa0.bin";
	if (work_in_buffer) {
		GSA32 = (uint32_t*) buffer;
		fillGnrSuffixArray32(GSA32, REV, spos, sid, 0, n);
		GSA32[n] = 0;
		writeArray32(GSA32, n + 1, gsa_fn);
	} else {
		writeArray64((uint64_t *) SA, n / 2 + 1, sa_fn);
		GSA32 = (uint32_t*) SA;	
		fillGnrSuffixArray32(GSA32, REV, spos, sid, 0, n);
		GSA32[n] = 0;
		writeArray32(GSA32, n + 1, gsa_fn);
		readArray64((uint64_t *) SA, n / 2 + 1, sa_fn);
		if (remove("sa0.bin") == 0)
			fprintf(stderr, "Deleted temp file sa0.bin.\n");
		else
			fprintf(stderr, "Error in deleting file sa0.bin.\n");
	}
	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing generalized suffix array: %lu ms.\n", duration);
	}
}

void SuffixArray::computeLcpArray(uint8_t *s, uint64_t n, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	uint64_t len = 0;
	LCP = (uint16_t*) buffer;
	//memset(LCP, 0, sizeof(uint16_t) * n);
	#pragma omp parallel for firstprivate(len)
	for (uint64_t i = 0; i < n; i++) {
		//k: 0 - n - 1; SA: 0 - n - 1
		uint64_t k = REV[i];
		if (k == 0)
			continue;
		uint64_t j = SA[k - 1];
		for(; i + len < n && j + len < n && s[i + len] == s[j + len]; len++);
		LCP[k] = (len >= 0xFFFFull) ? UINT16_MAX : (uint16_t)(len & 0xFFFF);
		if (len > 0) len--;
	}
	LCP[n] = 0;
	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing LCP array: %lu ms.\n", duration);
	}
}

void SuffixArray::computeAvgLcp(uint8_t* s, uint64_t n, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	uint64_t totalLCP = 0;
	double avgLCP = 0.0;
	uint64_t len = 0;
	for (uint64_t i = 0; i < n; i++) {
		//k: 0 - n - 1; SA: 0 - n - 1
		uint64_t k = REV[i];
		if (k == 0)
			continue;
		uint64_t j = SA[k - 1];
		for(; i + len < n && j + len < n && s[i + len] == s[j + len]; len++);
		totalLCP += len;
		if (len > 0) len--;
	}
	avgLCP = 1.0 * totalLCP / n;
	fprintf(stderr, "The average LCP is %lf.\n", avgLCP);
	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing avg LCP: %lu ms.\n", duration);
	}
}

void SuffixArray::reloadLCP(uint64_t n, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	std::string lcp_fn = "lcp.bin";
	readArray16(LCP, n + 1, lcp_fn);
	if (remove("lcp.bin") == 0)
		fprintf(stderr, "Deleted temp file lcp.bin.\n");
	else
		fprintf(stderr, "Error in deleting file lcp.bin.\n");
	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for reloading LCP array: %lu ms.\n", duration);
	}
}

void SuffixArray::prepareGnrLCP(uint64_t n, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();

	lcp_prepared__ = 1;
	std::string gsa_fn = "gsa.bin";
	if (gsa_unit_ == 32) {
		readArray32((uint32_t *) REV, n + 1, gsa_fn);
		GSA32 = (uint32_t *) REV;
		GSA32_ = new uint32_t[n + 2];
	} else {
		readArray16((uint16_t *) REV, n + 1, gsa_fn);
		GSA16 = (uint16_t *) REV;
		GSA16_ = (uint16_t *) (REV + n / 4 + 2);
	} 
	if (remove("gsa.bin") == 0)
		fprintf(stderr, "Deleted temp file gsa.bin.\n");
	else
		fprintf(stderr, "Error in deleting file gsa.bin.\n");

	memcpy(REV + n / 2 + 4, LCP, (n + 1) * sizeof(uint16_t));
	LCP = (uint16_t*) (REV + n / 2 + 4);
	LCP0 = (uint16_t*) buffer;

	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for resetting memory for LCP0 computation: %lu ms.\n", 
			duration);
	}
}

void SuffixArray::computeGnrLcpArray16(uint64_t n, uint16_t el, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	uint16_t minlcp = UINT16_MAX;
	uint64_t nextd = 0, begin = 0, end = n - 1;
	for (; GSA16[end] == GSA16[end - 1]; end--);

	for (uint64_t i = begin; i < end; i += (nextd + 1), minlcp = UINT16_MAX, nextd = 0) {
		for (; GSA16[i + nextd] == GSA16[i + nextd + 1]; nextd++);
		for (int64_t j = nextd; j >= 0; j--) {
			minlcp = std::min(minlcp, LCP[i + j + 1]);
			LCP0[i + j] = std::max(el, minlcp);
		}
	}

	for (uint64_t i = end; i < n; i++)
		LCP0[i] = 0;
	minlcp = UINT16_MAX;
	nextd = 0;
	for (end = 0; GSA16[end] == GSA16[end + 1]; end++);
	begin = n - 1;

	for (uint64_t i = begin; i > end; i -= (nextd + 1), minlcp = UINT16_MAX, nextd = 0) {
		for (; GSA16[i - nextd] == GSA16[i - nextd - 1]; nextd++);
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

void SuffixArray::computeGnrLcpArray32(uint64_t n, uint16_t el, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	uint16_t minlcp = UINT16_MAX;
	uint64_t nextd = 0, begin = 0, end = n - 1;
	for (; GSA32[end] == GSA32[end - 1]; end--);

	for (uint64_t i = begin; i < end; i += (nextd + 1), minlcp = UINT16_MAX, nextd = 0) {
		for (; GSA32[i + nextd] == GSA32[i + nextd + 1]; nextd++);
		for (int64_t j = nextd; j >= 0; j--) {
			minlcp = std::min(minlcp, LCP[i + j + 1]);
			LCP0[i + j] = std::max(el, minlcp);
		}
	}
	
	for (uint64_t i = end; i < n; i++)
		LCP0[i] = 0;
	minlcp = UINT16_MAX;
	nextd = 0;
	for (end = 0; GSA32[end] == GSA32[end + 1]; end++);
	begin = n - 1;

	for (uint64_t i = begin; i > end; i -= (nextd + 1), minlcp = UINT16_MAX, nextd = 0) {
		for (; GSA32[i - nextd] == GSA32[i - nextd - 1]; nextd++);
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

void SuffixArray::computeGnrLcpArray16_d(uint64_t n, uint16_t el, uint16_t ulmax, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	uint16_t minlcp = UINT16_MAX;
	uint64_t nextd = 0, begin = 0, end = n - 1;
	for (; GSA16[end] == GSA16[end - 1]; end--);
	memset(GSA16_, 0, sizeof(uint16_t) * (n + 1));

	#pragma omp parallel private(minlcp, nextd)
	{
		int threadnum = omp_get_thread_num(), 
			numthreads = omp_get_num_threads();
		uint64_t l = n * threadnum / numthreads, 
			r = std::min(end, n * (threadnum + 1) / numthreads);
		
		if (l > begin)
			for (; GSA16[l] == GSA16[l - 1]; l--);

		minlcp = UINT16_MAX;
		nextd = 0;
		for (uint64_t i = l; i < r; i += (nextd + 1), minlcp = UINT16_MAX, nextd = 0) {
			for (; GSA16[i + nextd] == GSA16[i + nextd + 1]; nextd++);
			for (int64_t j = nextd; j >= 0; j--) {
				minlcp = std::min(minlcp, LCP[i + j + 1]);
				LCP0[i + j] = minlcp;
				GSA16_[SA[i + j]] = GSA16[i + nextd + 1];
			}
		}
	}

	for (uint64_t i = end; i < n; i++)
		LCP0[i] = 0;
	minlcp = UINT16_MAX;
	nextd = 0;
	begin = n - 1;
	for (end = 0; GSA16[end] == GSA16[end + 1]; end++);
	for (; GSA16[end] == GSA16[end + 1]; end++);

	#pragma omp parallel private(minlcp, nextd)
	{
		int threadnum = omp_get_thread_num(), 
			numthreads = omp_get_num_threads();
		uint64_t l = std::max(end, n * threadnum / numthreads), 
			r = n * (threadnum + 1) / numthreads;
		if (r == n)
			r = begin;
		else
			if (GSA16[r] == GSA16[r + 1])
				for (; GSA16[r] == GSA16[r - 1]; r--);

		minlcp = UINT16_MAX;
		nextd = 0;
		for (uint64_t i = r; i > l; i -= (nextd + 1), minlcp = UINT16_MAX, nextd = 0) {
			for (; GSA16[i - nextd] == GSA16[i - nextd - 1]; nextd++);
			for (int64_t j = nextd; j >= 0; j--) {
				minlcp = std::min(minlcp, LCP[i - j]);
				if (LCP0[i - j] < minlcp) {
					uint16_t min2lcp = UINT16_MAX;
					uint64_t i_ = i - nextd - 1;
					for (; GSA16[i_] == GSA16[i_ - 1]; i_--)	
						min2lcp = std::min(min2lcp, LCP[i_]);
					min2lcp = std::min(min2lcp, LCP[i_]);
					min2lcp = std::min(min2lcp, minlcp);
					LCP0[i - j] = std::max(LCP0[i - j], min2lcp);
					LCP0[i - j] = std::max(LCP0[i - j], el);
					GSA16_[SA[i - j]] = GSA16[i - nextd - 1];
					if (LCP0[i - j] >= minlcp)
						LCP0[i - j] = ulmax + 2;
				} else {
					if (LCP0[i - j] > minlcp) {
						uint16_t min2lcp = UINT16_MAX;
						uint64_t i_ = i;
						for (; GSA16[i_] == GSA16[i_ + 1] && i_ < n; i_++)
							min2lcp = std::min(min2lcp, LCP[i_ + 1]);
						min2lcp = std::min(min2lcp, LCP[i_ + 1]); 
						for (i_ = i_ + 1; GSA16[i_] == GSA16[i_ + 1] && i_ < n; i_++)
							min2lcp = std::min(min2lcp, LCP[i_ + 1]);
						min2lcp = std::min(min2lcp, LCP[i_ + 1]);
						uint16_t lcp0 = std::max(minlcp, min2lcp);
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
	}

	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing LCP0-D array: %lu ms.\n", duration);
	}
}

void SuffixArray::computeGnrLcpArray32_d(uint64_t n, uint16_t el, uint16_t ulmax, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	uint16_t minlcp = UINT16_MAX;
	uint64_t nextd = 0, begin = 0, end = n - 1;
	for (; GSA32[end] == GSA32[end - 1]; end--);
	memset(GSA32_, 0, sizeof(uint32_t) * (n + 2));
	
	#pragma omp parallel private(minlcp, nextd)
	{
		int threadnum = omp_get_thread_num(), 
			numthreads = omp_get_num_threads();
		uint64_t l = n * threadnum / numthreads, 
			r = std::min(end, n * (threadnum + 1) / numthreads);

		if (l > begin)
			for (; GSA32[l] == GSA32[l - 1]; l--);

		minlcp = UINT16_MAX;
		nextd = 0;
		for (uint64_t i = l; i < r; i += (nextd + 1), minlcp = UINT16_MAX, nextd = 0) {
			for (; GSA32[i + nextd] == GSA32[i + nextd + 1]; nextd++);
			for (int64_t j = nextd; j >= 0; j--) {
				minlcp = std::min(minlcp, LCP[i + j + 1]);
				LCP0[i + j] = minlcp;
				GSA32_[SA[i + j]] = GSA32[i + nextd + 1];
			}
		}
	}

	for (uint64_t i = end; i < n; i++)
		LCP0[i] = 0;
	minlcp = UINT16_MAX;
	nextd = 0;
	begin = n - 1;
	for (end = 0; GSA32[end] == GSA32[end + 1]; end++);
	for (; GSA32[end] == GSA32[end + 1]; end++);

	#pragma omp parallel private(minlcp, nextd)
	{
		int threadnum = omp_get_thread_num(), 
			numthreads = omp_get_num_threads();
		uint64_t l = std::max(end, n * threadnum / numthreads), 
			r = n * (threadnum + 1) / numthreads;
		if (r == n)
			r = begin;
		else
			if (GSA32[r] == GSA32[r + 1])
				for (; GSA32[r] == GSA32[r - 1]; r--);

		minlcp = UINT16_MAX;
		nextd = 0;
		for (uint64_t i = r; i > l; i -= (nextd + 1), minlcp = UINT16_MAX, nextd = 0) {
			for (; GSA32[i - nextd] == GSA32[i - nextd - 1]; nextd++);
			for (int64_t j = nextd; j >= 0; j--) {
				minlcp = std::min(minlcp, LCP[i - j]);
				if (LCP0[i - j] < minlcp) {
					uint16_t min2lcp = UINT16_MAX;
					uint64_t i_ = i - nextd - 1;
					for (; GSA32[i_] == GSA32[i_ - 1]; i_--)	
						min2lcp = std::min(min2lcp, LCP[i_]);
					min2lcp = std::min(min2lcp, LCP[i_]);
					min2lcp = std::min(min2lcp, minlcp);
					LCP0[i - j] = std::max(LCP0[i - j], min2lcp);
					LCP0[i - j] = std::max(LCP0[i - j], el);
					GSA32_[SA[i - j]] = GSA32[i - nextd - 1];
					if (LCP0[i - j] >= minlcp)
						LCP0[i - j] = ulmax + 2;
				} else {
					if (LCP0[i - j] > minlcp) {
						uint16_t min2lcp = UINT16_MAX;
						uint64_t i_ = i;
						for (; GSA32[i_] == GSA32[i_ + 1] && i_ < n; i_++)
							min2lcp = std::min(min2lcp, LCP[i_ + 1]);
						min2lcp = std::min(min2lcp, LCP[i_ + 1]); 
						for (i_ = i_ + 1; GSA32[i_] == GSA32[i_ + 1] && i_ < n; i_++)
							min2lcp = std::min(min2lcp, LCP[i_ + 1]);
						min2lcp = std::min(min2lcp, LCP[i_ + 1]);
						uint16_t lcp0 = std::max(minlcp, min2lcp);
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
	}

	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing LCP0-D array: %lu ms.\n", duration);
	}
}

void SuffixArray::computeMinUnique(uint64_t n, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	uint16_t *MU = LCP;
	memset(MU, UINT8_MAX, sizeof(uint16_t) * (n + 1));

	#pragma omp parallel for
	for (uint64_t i = 0; i < n; i++) {
		uint64_t SA_i = SA[i];
		MU[SA_i + LCP0[i] + 1] = std::min(MU[SA_i + LCP0[i] + 1], LCP0[i]);
	}

	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing minimum unique substrings: %lu ms.\n", duration);
	}
}

void SuffixArray::computeMinUnique(uint64_t n, uint16_t ulmax, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	uint16_t *MU = LCP;
	memset(MU, UINT8_MAX, sizeof(uint16_t) * (n + 1));

	#pragma omp parallel for
	for (uint64_t i = 0; i < n; i++) {
		//if (SA[i] + LCP0[i] + 1 > n)
		//	fprintf(stderr, "%lu, %ld, %lu, %u\n", i, SA[i], SA[i] + LCP0[i] + 1, LCP0[i]);
		uint64_t SA_i = SA[i];
		if (LCP0[i] < ulmax)
			MU[SA_i + LCP0[i] + 1] = std::min(MU[SA_i + LCP0[i] + 1], LCP0[i]);
	}

	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing minimum unique substrings: %lu ms.\n", duration);
	}
}

void SuffixArray::computeOCC16(uint64_t n, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	occ = (uint8_t*) (REV + 3 * n / 4 + 6);
	memset(occ, 1, sizeof(uint8_t) * (n + 1));

	#pragma omp parallel for
	for (uint64_t i = 0; i < n - 1; i++) {
		uint64_t SA_i = SA[i];
		uint16_t minlcp = LCP[i + 1];
		uint64_t j = 0;
		while ((i + j <= n - 1) && (GSA16[i + j + 1] == GSA16[i]) && (minlcp > LCP0[i])) {
			occ[SA_i]++;
			j++;
			minlcp = std::min(minlcp, LCP[i + j + 1]);
		}
	}

	#pragma omp parallel for
	for (uint64_t i = n - 1; i > 0; i--) {
		uint64_t SA_i = SA[i];
		uint16_t minlcp = LCP[i];
		int64_t j = 0;
		while ((i - j > 0) && (GSA16[i - j - 1] == GSA16[i]) && (minlcp > LCP0[i])) {
			occ[SA_i]++;
			j++;
			minlcp = std::min(minlcp, LCP[i - j]);
		}
	}

	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing OCC array: %lu ms.\n", duration);
	} 
}

void SuffixArray::computeOCC32(uint64_t n, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	occ = (uint8_t*) (REV + 3 * n / 4 + 6);
	memset(occ, 1, sizeof(uint8_t) * (n + 1));
	
	#pragma omp parallel for
	for (uint64_t i = 0; i < n - 1; i++) {
		uint64_t SA_i = SA[i];
		uint16_t minlcp = LCP[i + 1];
		uint64_t j = 0;
		while ((i + j <= n - 1) && (GSA32[i + j + 1] == GSA32[i]) && (minlcp > LCP0[i])) {
			occ[SA_i]++;
			j++;
			minlcp = std::min(minlcp, LCP[i + j + 1]);
		}
	}

	#pragma omp parallel for
	for (uint64_t i = n - 1; i > 0; i--) {
		uint64_t SA_i = SA[i];
		uint16_t minlcp = LCP[i];
		int64_t j = 0;
		while ((i - j > 0) && (GSA32[i - j - 1] == GSA32[i]) && (minlcp > LCP0[i])) {
			occ[SA_i]++;
			j++;
			minlcp = std::min(minlcp, LCP[i - j]);
		}
	}

	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing OCC array: %lu ms.\n", duration);
	} 
}

void SuffixArray::computeOCC16_d(uint64_t n, uint16_t ulmax, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	occ = (uint8_t*) (REV + 3 * n / 4 + 6);
	occ2 = occ + n + 2;
	memset(occ, 0, sizeof(uint8_t) * (n + 1));
	memset(occ2, 0, sizeof(uint8_t) * (n + 1));

	uint16_t minlcp = UINT16_MAX;
	uint64_t begin = n - 1, end = 0;
	for (; GSA16[end] == GSA16[end + 1]; end++);
	for (; GSA16[end] == GSA16[end + 1]; end++);
	
	#pragma omp parallel for private(minlcp)
	for (uint64_t i = begin; i > end; i--) {
		if (LCP0[i] <= ulmax) {
			uint64_t SA_i = SA[i];
			occ[SA_i] = 1;
			minlcp = UINT16_MAX;
			for (uint64_t j = 0; (i - j > end) && ((GSA16[i - j - 1] == GSA16[i]) 
					|| (GSA16[i - j - 1] == GSA16_[SA_i])); j++) {
				minlcp = std::min(minlcp, LCP[i - j]);
				if (minlcp > LCP0[i]) {
					if (GSA16[i - j - 1] == GSA16[i])
						occ[SA_i]++;
					if (GSA16[i - j - 1] == GSA16_[SA_i])
						occ2[SA_i]++;
				}
			}
			minlcp = UINT16_MAX;
			for (uint64_t j = 0; (i + j <= begin) && ((GSA16[i + j + 1] == GSA16[i]) 
					|| (GSA16[i + j + 1] == GSA16_[SA_i])); j++) {
				minlcp = std::min(minlcp, LCP[i + j + 1]);
				if (minlcp > LCP0[i]) {
					if (GSA16[i + j + 1] == GSA16[i])
						occ[SA_i]++;
					if (GSA16[i + j + 1] == GSA16_[SA_i])
						occ2[SA_i]++;				
				}
			}
		}
	}

	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing OCC array: %lu ms.\n", duration);
	}
}

void SuffixArray::computeOCC32_d(uint64_t n, uint16_t ulmax, bool debug) {
	auto start = std::chrono::high_resolution_clock::now();
	occ = (uint8_t*) (REV + 3 * n / 4 + 6);
	occ2 = occ + n + 2;
	memset(occ, 0, sizeof(uint8_t) * (n + 1));
	memset(occ2, 0, sizeof(uint8_t) * (n + 1));

	uint16_t minlcp = UINT16_MAX;
	uint64_t begin = n - 1, end = 0;
	for (end = 0; GSA32[end] == GSA32[end + 1]; end++);
	for (; GSA32[end] == GSA32[end + 1]; end++);
	
	#pragma omp parallel for private(minlcp)
	for (uint64_t i = begin; i > end; i--) {
		if (LCP0[i] <= ulmax) {
			uint64_t SA_i = SA[i];
			occ[SA_i] = 1;
			minlcp = UINT16_MAX;
			for (uint64_t j = 0; (i - j > end) && ((GSA32[i - j - 1] == GSA32[i]) 
					|| (GSA32[i - j - 1] == GSA32_[SA_i])); j++) {
				minlcp = std::min(minlcp, LCP[i - j]);
				if (minlcp > LCP0[i]) {
					if (GSA32[i - j - 1] == GSA32[i])
						occ[SA_i]++;
					if (GSA32[i - j - 1] == GSA32_[SA_i])
						occ2[SA_i]++;
				}
			}
			minlcp = UINT16_MAX;
			for (uint64_t j = 0; (i + j <= begin) && ((GSA32[i + j + 1] == GSA32[i]) 
					|| (GSA32[i + j + 1] == GSA32_[SA_i])); j++) {
				minlcp = std::min(minlcp, LCP[i + j + 1]);
				if (minlcp > LCP0[i]) {
					if (GSA32[i + j + 1] == GSA32[i])
						occ[SA_i]++;
					if (GSA32[i + j + 1] == GSA32_[SA_i])
						occ2[SA_i]++;				
				}
			}
		}
	}

	if (debug) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
				(std::chrono::high_resolution_clock::now() - start).count();
		fprintf(stderr, "Time for computing OCC array: %lu ms.\n", duration);
	}
}

void SuffixArray::readArray64(uint64_t *arr, uint64_t n, std::string &fn) {
	std::ifstream istream;
	try { 
		istream.open(fn.c_str(), std::ios::in | std::ios::binary);
		istream.read((char*) arr, n * sizeof(uint64_t));
		istream.close();
	} catch (...) {
		fprintf(stderr, "Error in reading array %s.\n", fn.c_str());
		abort();
	}
}

void SuffixArray::readArray32(uint32_t *arr, uint64_t n, std::string &fn) {
	std::ifstream istream;
	try { 
		istream.open(fn.c_str(), std::ios::in | std::ios::binary);
		istream.read((char*) arr, n * sizeof(uint32_t));
		istream.close();
	} catch (...) {
		fprintf(stderr, "Error in reading array %s.\n", fn.c_str());
		abort();
	}
}

void SuffixArray::readArray16(uint16_t *arr, uint64_t n, std::string &fn) {
	std::ifstream istream;
	try { 
		istream.open(fn.c_str(), std::ios::in | std::ios::binary);
		istream.read((char*) arr, n * sizeof(uint16_t));
		istream.close();
	} catch (...) {
		fprintf(stderr, "Error in reading array %s.\n", fn.c_str());
		abort();
	}
}

void SuffixArray::writeArray64(uint64_t *arr, uint64_t n, std::string &fn) {
	std::ofstream ostream;
	try { 
		ostream.open(fn.c_str(), std::ios::out | std::ios::binary);
		ostream.write((char*) arr, n * sizeof(uint64_t));
		ostream.close();
	} catch (...) {
		fprintf(stderr, "Error in writing array %s.\n", fn.c_str());
		abort();
	}
}

void SuffixArray::writeArray32(uint32_t *arr, uint64_t n, std::string &fn) {
	std::ofstream ostream;
	try { 
		ostream.open(fn.c_str(), std::ios::out | std::ios::binary);
		ostream.write((char*) arr, n * sizeof(uint32_t));
		ostream.close();
	} catch (...) {
		fprintf(stderr, "Error in writing array %s.\n", fn.c_str());
		abort();
	}
}

void SuffixArray::writeArray16(uint16_t *arr, uint64_t n, std::string &fn) {
	std::ofstream ostream;
	try { 
		ostream.open(fn.c_str(), std::ios::out | std::ios::binary);
		ostream.write((char*) arr, n * sizeof(uint16_t));
		ostream.close();
	} catch (...) {
		fprintf(stderr, "Error in writing array %s.\n", fn.c_str());
		abort();
	}
}

uint16_t* SuffixArray::run(uint8_t* s, std::vector<uint64_t> &spos, 
				std::vector<uint32_t> &sid, uint64_t n, int mode, 
				uint16_t ulmax, uint16_t el, bool debug, 
				bool debug_sa, bool write_lcp) {
	switch (mode) {
		case 0:
			computeSuffixArray(s, n, debug, debug_sa);
			computeRevSuffixArray(n, debug);
			allocBuffer(n / 4 + 2);
			if (gsa_unit_ == 32)
				computeGnrSuffixArray32(spos, sid, n, debug, 0);
			else
				computeGnrSuffixArray16(spos, sid, n, debug, 1);
			computeLcpArray(s, n, debug);
			return NULL;
		case 1:
			prepareGnrLCP(n, debug);
			if (gsa_unit_ == 32) {
				computeGnrLcpArray32(n, el, debug);
				computeOCC32(n, debug);
			} else {
				computeGnrLcpArray16(n, el, debug);
				computeOCC16(n, debug);
			}
			if (write_lcp) {
				std::string lcp_fn = "lcp.bin";
				writeArray16(LCP, n + 1, lcp_fn);
			}
			computeMinUnique(n, debug);
			break;
		case 2:
			if (lcp_prepared__ == 0)
				prepareGnrLCP(n, debug);
			else
				reloadLCP(n, debug);
			if (gsa_unit_ == 32) {
				computeGnrLcpArray32_d(n, el, ulmax, debug);
				computeOCC32_d(n, ulmax, debug);
			} else {
				computeGnrLcpArray16_d(n, el, ulmax, debug);
				computeOCC16_d(n, ulmax, debug);
			}
			computeMinUnique(n, ulmax, debug);
			break;
	}
	
	return LCP;
}

