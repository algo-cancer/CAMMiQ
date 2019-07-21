#include <dirent.h>
#include <pthread.h>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <cstdlib>

#include "util.hpp"
#include "build.hpp"

int FastaReader::mode_;

FastaReader::FastaReader(int L_, int K, int t) {
	L = L_;
	minL = K;
	num_threads = t;
	if (!filenames.empty()) filenames.clear();
	if (!contig_pos.empty()) contig_pos.clear();
	if (!ref_pos.empty()) ref_pos.clear();
	if (!refID.empty()) refID.clear();
	if (!seqs.empty()) seqs.clear();
	uint8_t *seq0 = new uint8_t[block_size_];
	seqs.push_back(seq0);
}

FastaReader::FastaReader(int L_, int K, int t, std::vector<std::string>& fns_) {
	L = L_;
	minL = K;
	num_threads = t;
	if (!filenames.empty()) filenames.clear();
	for (auto fn : fns_)
		filenames[fn] = 0;
	if (!contig_pos.empty()) contig_pos.clear();
	if (!ref_pos.empty()) ref_pos.clear();
	if (!refID.empty()) refID.clear();
	if (!seqs.empty()) seqs.clear();
	uint8_t *seq0 = new uint8_t[block_size_];
	seqs.push_back(seq0);
}

void FastaReader::allocSeq() {
	uint8_t *seq = new uint8_t[block_size_];
	seqs.push_back(seq);
	cur++;
	pos = 0;
}

void FastaReader::prepFasta(std::string &INDIR) {
	if (!filenames.empty()) {
		filenames.clear();
		fprintf(stderr, "Flushed existing file names.\n");
	}
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir(INDIR.c_str())) != NULL) {
		while ((ent = readdir(dir)) != NULL) {
			std::string filename = ent->d_name;
			if (filename.length() >= 7) {
				std::string fext = filename.substr(filename.length() - 6);
				if (strcmp(fext.c_str(), ".fasta") == 0)
					filenames[INDIR + filename] = 0;
				else {
					fext = fext.substr(fext.length() - 4);
					if ((strcmp(fext.c_str(), ".fna") == 0) ||
						(strcmp(fext.c_str(), ".ffn") == 0))
						filenames[INDIR + filename] = 0;
				}
			}
		}
		closedir(dir);
	} else {
		/* could not open directory */
		fprintf(stderr, "Input directory not exists.\n");
		abort();
	}
}

void FastaReader::readAllFasta() {
	for (auto it : filenames) {
		std::string fn = it.first;
		readFasta(fn);
		refID.push_back(it.second);
	}
	assert(refID.size() == filenames.size());
	fprintf(stderr, "\n****************************\n");
	fprintf(stderr, "Total num bases: %lu\n", N_);
	fprintf(stderr, "Total num species: %u\n", M_);
	fprintf(stderr, "Total num contigs: %u\n", C_);
	fprintf(stderr, "****************************\n");
}

void FastaReader::readFnMap(std::string &INDIR, std::string &INFILE) {
	std::string line, fn, sp;
	std::ifstream inputFile;

	/* Read in fasta file. */ 
	inputFile.open(INFILE);
	while (std::getline(inputFile, line)) {
		std::istringstream lines(line);
		lines >> fn;
		lines >> sp;
		filenames[INDIR + fn] = stoi(sp);
	}
	inputFile.close();
	fprintf(stderr, "Prepared species information of every fasta file.\n");
}

void FastaReader::readFasta(std::string &INFILE) { 
 	std::string bases;
	std::ifstream inputFile;

	/* Read in fasta file. */ 
	inputFile.open(INFILE);
	while (std::getline(inputFile, bases)) {
		if (bases[0] == '>') {
			if (N_ > 0 && seqs[cur][pos - 1] >= base_offset) {
				if (!insertContig()) {
					inputFile.close();
					abort();
				}
			}
		} else {
			if (!insertBase(bases)) {
				inputFile.close();
				abort();
			}
		}
	}
	if (!insertContig()) {
		inputFile.close();
		abort();
	}
	ref_pos.push_back(cur * block_size_ + pos);
	M_++;
	if (M_ >= maxM) {
		fprintf(stderr, "Number of reference genomes exceeds limit.\n");
		inputFile.close();
		abort();
	}
	inputFile.close();
	
	//fprintf(stderr,	"Reference sequence %s loaded.\r", INFILE.c_str()); 
}

bool FastaReader::insertBase(std::string& bases) {
	if (N_ + bases.length() >= maxN) {
		fprintf(stderr, "Total number of symbols exceeds limit.\n");
		return false;
	}
	if (pos + bases.length() > block_size_) {
		size_t i = 0, j = block_size_ - pos;
		for (; i < j; i++)
			seqs[cur][pos++] = bases[i] + base_offset;
		allocSeq();
		for (; i < bases.length(); i++)
			seqs[cur][pos++] = bases[i] + base_offset;
	} else {
		for (size_t i = 0; i < bases.length(); i++)
			seqs[cur][pos++] = bases[i] + base_offset;
	}
	N_ += bases.length();
	return true;
}

bool FastaReader::insertContig() {
	if (pos + 4 > block_size_) {
		size_t i = 3, j = pos + 3 - block_size_;
		for (; i > j; i--)
			seqs[cur][pos++] = (C_ >> (i * 7)) & 0x7F;
		allocSeq();
		for (int i_ = i; i_ >= 0; i_--, i--)
			seqs[cur][pos++] = (C_ >> (i * 7)) & 0x7F;
	} else {
		for (int i = 3; i >= 0; i--)
			seqs[cur][pos++] = (C_ >> (i * 7)) & 0x7F;
	}

	contig_pos.push_back(cur * block_size_ + pos);
	//fprintf(stderr, "cl: %lu\n",  contig_pos[contig_pos.size() - 1] - contig_pos[contig_pos.size() - 2]);

	C_++;
	if (C_ >= maxC) {
		fprintf(stderr, "Number of contigs exceeds limit.\n");
		return false;
	}
	return true;
}

void FastaReader::allocSuffixArray() {
	/* Allocate the memory space of input string. */
	uint64_t N = cur * block_size_ + pos;
	seqs_ = new uint8_t[N + 2];
	for (size_t i = 0; i < cur; i++)
		memcpy(seqs_ + i * block_size_, seqs[i], block_size_);
	memcpy(seqs_ + cur * block_size_, seqs[cur], pos);
	fprintf(stderr, "Input string get prepared.\n");
	/*size_t j = 0, cc = 0;
	for (size_t i = 0; i < ref_pos.size(); i++, cc = 0) {
		
		while (contig_pos[j] < ref_pos[i]) {
			cc++;
			j++;
		}

		if (i == 0)
			fprintf(stderr, "%u\t%lu\n", refID[i], ref_pos[i] - (cc + 1) * L - (cc + 1) * 4);
		else
			fprintf(stderr, "%u\t%lu\n", refID[i], ref_pos[i] - (cc + 1) * L - ref_pos[i - 1] - (cc + 1) * 4);
	}*/

	//for (size_t i = 0; i < 10; i++)
	//	fprintf(stderr, "%lu\n", ref_pos[i]);
	/*std::string test = "GGACCCCTAATTA";
	uint8_t *test_ = new uint8_t[test.length()];
	for (size_t i = 0; i < test.length(); i++)
		test_[i] = test[i] + 165;
	for (uint64_t i = 0; i < N; i++)
		if (memcmp(seqs_ + i, test_, test.length()) == 0) {
			fprintf(stderr, "Found at pos %lu\t", i);
			for (int j = 0; j < (int)(test.length() + 3); j++)
				fprintf(stderr, "%c", (seqs_ + i)[j - 1] - 165);
			fprintf(stderr, "\n");
		}*/

	/* Compute suffix arrays. */	
	seqs_[N] = 0;
	seqs_[N + 1] = 0;
	sa = new SuffixArray(N);
	mu_index = sa->run(seqs_, ref_pos, refID, N);
}

void* FastaReader::computeIndexstats() {
	int tid = 0;
	pthread_mutex_lock(&thread_lock);
	tid = current_tid++;
	pthread_mutex_unlock(&thread_lock);
	fprintf(stderr, "Begin computing the statistics by thread %d.\n", tid + 1);

	uint64_t i = 1, j = 0, nexti = 1;
	size_t nref = ref_pos.size() / num_threads;
	if (tid > 0)
		i = ref_pos[tid * nref - 1];
	nexti = (tid == num_threads - 1) ? ref_pos[ref_pos.size() - 1] : ref_pos[(tid + 1) * nref - 1];
	size_t ci = 0, ri = tid * nref;
	for (; i >= contig_pos[ci]; ci++);

	for (; i < nexti; i++, ci += (i >= contig_pos[ci]), ri += (i >= ref_pos[ri])) {
		j = mu_index[i];

		// Index not assigned or at a separation region between two contigs.
		if ((j == 0) || (i >= contig_pos[ci] - 4)) {
			i++;
			continue;
		}

		// Substring spans two different contigs.
		if (ci > 0 && j - 1 < contig_pos[ci - 1]) {
			i++;
			continue;
		}

		// Substring have special characters "RYKMSWBDHVN-".
		bool special_character_ = false;
		for (uint64_t spos = j - 1; spos < i; spos++)
			if ((seqs_[spos] - base_offset != 65) && 
			    (seqs_[spos] - base_offset != 67) &&
			    (seqs_[spos] - base_offset != 71) &&
			    (seqs_[spos] - base_offset != 84)) {
				special_character_ = true;
				break;
			}
		if (special_character_) {
			i++;
			continue;
		}

		// Substring length exceeds limit.
		if (i - j >= 0xFFFFul) {
			i++;
			continue;
		}
		
		// Update the statistics.
		uint64_t length = i - j + 1;
		minusl[ri] = std::min(minusl[ri], length);
		maxusl[ri] = std::max(maxusl[ri], length);
		sumusl[ri] += length;
		uscount[ri]++;
	
		// Print only positions.
		/*if (!show_substring_) {
			fprintf(pFile, "i: %lu, ref: %lu, Length: %lu\n", mu_index[i] - 1, rpos, length);
		} else { // Print the actual substrings.
			fprintf(pFile, "i: %lu, ref: %lu, Length: %lu, ", mu_index[i] - 1, rpos, length);
			for (uint64_t start = mu_index[i] - 1; start < i; start++)
				fprintf(pFile, "%c", seqs_[start] - base_offset);
			fprintf(pFile, "\n");
		}*/
	}
	return 0;
}

void FastaReader::insert32(uint64_t i, uint32_t length, uint32_t rid) {
	uint32_t bucket = hasht->computeHashVal(seqs_ + i);
	int lock_num = ((int) bucket) & (num_locks - 1);
	pthread_spin_lock(hasht_access[lock_num]);
	hasht->insert32(bucket, seqs_ + i, length, rid);
	pthread_spin_unlock(hasht_access[lock_num]);
}

void FastaReader::insert32_d(uint64_t i, uint32_t length, uint32_t rid) {
	uint32_t bucket = hasht->computeHashVal(seqs_ + i);
	int lock_num = ((int) bucket) & (num_locks - 1);
	pthread_spin_lock(hasht_access[lock_num]);
	/*for (int j = 0; j < length; j++)
		fprintf(stderr, "%c", (seqs_ + i)[j] - 165);
	fprintf(stderr, "\t%u\n", rid);*/
	hasht->insert32_d(bucket, seqs_ + i, length, rid);
	pthread_spin_unlock(hasht_access[lock_num]);
}

void FastaReader::insert64(uint64_t i, uint32_t length, uint32_t rid) {
	uint64_t bucket = hasht->computeHashVal64(seqs_ + i);
	int lock_num = ((int) bucket) & (num_locks - 1);
	pthread_spin_lock(hasht_access[lock_num]);
	hasht->insert64(bucket, seqs_ + i, length, rid);
	pthread_spin_unlock(hasht_access[lock_num]);
}

void* FastaReader::analIndexCLARK() {
	int tid = 0;
	pthread_mutex_lock(&thread_lock);
	tid = current_tid++;
	pthread_mutex_unlock(&thread_lock);
	fprintf(stderr, "Begin counting unique L-mers by thread %d, pid %lu.\n", tid + 1, pthread_self());

	uint64_t i = 1, j = 0, nexti = 1, start = 0;
	size_t nref = ref_pos.size() / num_threads;
	if (tid > 0)
		i = ref_pos[tid * nref - 1];
	nexti = (tid == num_threads - 1) ? ref_pos[ref_pos.size() - 1] : ref_pos[(tid + 1) * nref - 1];
	size_t ci = 0, ri = tid * nref, lastr = ri;
	for (; i >= contig_pos[ci]; ci++);

	while (i < nexti) {
		// Index not assigned.
		j = mu_index[i];
		if (j == 0) {
			i++;
			continue;
		}

		// Reached a separation region between two contigs.
		if (i >= contig_pos[ci] - 4) {
			if ((start + L + 2 >= contig_pos[ci]) && exist_unique_[ci]) {
				if (ri == lastr)
					uLmcount_clark[ri] -= (start + L + 3 - contig_pos[ci]);
				else
					uLmcount_clark[lastr] -= (start + L + 3 - contig_pos[ci]);
			}
			start = std::max(contig_pos[ci], i - L);
			ci++;
			if (i >= ref_pos[ri] - 4)
				ri++;
			if (start + L + 2 >= contig_pos[ci])
				exist_unique_[ci] = false;
		}

		/*j = mu_index[i];

		// Index not assigned or at a separation region between two contigs.
		if ((j == 0) || (i >= contig_pos[ci] - 4)) {
			i++;
			continue;
		}*/

		// Substring spans two different contigs.
		if (ci > 0 && j - 1 < contig_pos[ci - 1]) {
			i++;
			continue;
		}

		// Substring have special characters "RYKMSWBDHVN-".
		bool special_character_ = false;
		for (uint64_t spos = j - 1; spos < i; spos++)
			if ((seqs_[spos] - base_offset != 65) && 
			    (seqs_[spos] - base_offset != 67) &&
			    (seqs_[spos] - base_offset != 71) &&
			    (seqs_[spos] - base_offset != 84)) {
				special_character_ = true;
				break;
			}
		if (special_character_) {
			i++;
			continue;
		}

		// Substring length exceeds L, is not considered within the index.
		uint64_t length = i - j + 1;
		if (length > L) {
			i++;
			continue;
		}

		// Process next reference species.
		/*if (j > ref_pos[ri])
			ri++;

		// Process next contig.
		if (j > contig_pos[ci]) {
			if ((start + L + 2 >= contig_pos[ci]) && exist_unique_[ci]) {
				if (ri == lastr)
					uLmcount_clark[ri] -= (start + L + 3 - contig_pos[ci]);
				else
					uLmcount_clark[lastr] -= (start + L + 3 - contig_pos[ci]);
			}
			start = std::max(contig_pos[ci], i - L);
			ci++;
			if (start + L + 2 >= contig_pos[ci]) {
				exist_unique_[ci] = false;
			}
		}*/

		// Aggregate unique L-mers.
		if (length <= K_) {
			if (i <= start + L) {
				uLmcount_clark[ri] += (j - start);
				start = j;
			} else {
				uLmcount_clark[ri] += (j + L - i);
				start = j;
			}
		}
		
		i++;
		lastr = ri;		
	}
	return 0;
}	

void* FastaReader::computeIndexno() {
	int tid = 0;
	pthread_mutex_lock(&thread_lock);
	tid = current_tid++;
	pthread_mutex_unlock(&thread_lock);
	fprintf(stderr, "Begin computing non-overlapping substrings by thread %d.\n", tid + 1);

	uint64_t i = 1, j = 0, nexti = 1, start = 0;
	size_t nref = ref_pos.size() / num_threads;
	if (tid > 0)
		i = ref_pos[tid * nref - 1];
	nexti = (tid == num_threads - 1) ? ref_pos[ref_pos.size() - 1] : ref_pos[(tid + 1) * nref - 1];
	size_t ci = 0, ri = tid * nref, lastr = ri;
	for (; i >= contig_pos[ci]; ci++);
	uint64_t lasti = 0;

	//fprintf(stderr, "NEXTI: %lu.\n", nexti);
	while (i < nexti) {
		// Index not assigned.
		j = mu_index[i];
		if (j == 0) {
			i++;
			continue;
		}

		// Reached a separation region between two contigs.
		if (i >= contig_pos[ci] - 4) {
			if ((start + L + 2 >= contig_pos[ci]) && exist_unique_[ci]) {
				if (ri == lastr)
					uLmcount[ri] -= (start + L + 3 - contig_pos[ci]);
				else
					uLmcount[lastr] -= (start + L + 3 - contig_pos[ci]);
			}
			start = std::max(contig_pos[ci], i - L);
			ci++;
			if (i >= ref_pos[ri] - 4)
				ri++;
			if (start + L + 2 >= contig_pos[ci]) {
				exist_unique_[ci] = false;
			}
		}

		// Substring spans two different contigs.
		if (ci > 0 && j - 1 < contig_pos[ci - 1]) {
			i++;
			continue;
		}

		// Substring have special characters "RYKMSWBDHVN-".
		bool special_character_ = false;
		for (uint64_t spos = j - 1; spos < i; spos++)
			if ((seqs_[spos] - base_offset != 65) && 
			    (seqs_[spos] - base_offset != 67) &&
			    (seqs_[spos] - base_offset != 71) &&
			    (seqs_[spos] - base_offset != 84)) {
				special_character_ = true;
				break;
			}
		if (special_character_) {
			i++;
			continue;
		}

		// Substring length exceeds L, is not considered within the index.
		uint64_t length = i - j + 1;
		if (length > L) {
			i++;
			continue;
		}

		// Process next reference species.
		/*if (j > ref_pos[ri])
			ri++;

		// Process next contig.
		if (j > contig_pos[ci]) {
			if ((start + L + 2 >= contig_pos[ci]) && exist_unique_[ci]) {
				if (ri == lastr)
					uLmcount[ri] -= (start + L + 3 - contig_pos[ci]);
				else
					uLmcount[lastr] -= (start + L + 3 - contig_pos[ci]);
			}
			start = std::max(contig_pos[ci], i - L);
			ci++;
			if (start + L + 2 >= contig_pos[ci]) {
				exist_unique_[ci] = false;
			}
		}*/

		// Index: (max) non-overlapping unique substrings on each contig.
		if (lasti == 0) {
			insert32(j - 1, length, refID[ri]);
			lasti = i;
			//std::pair<uint64_t, uint64_t> us_interval(j - 1, i);
			//nous[ri].push_back(us_interval);
		} else {
			if (j > lasti) {
				insert32(j - 1, length, refID[ri]);
				lasti = i;
			}
		}
		//fprintf(stderr, "%lu %lu %d\r", j, i, L);
		// Aggregate unique L-mers.
		if (i <= start + L) {
			uLmcount[ri] += (j - start);
			start = j;
		} else {
			uLmcount[ri] += (j + L - i);
			start = j;
		}
			
		i++;
		lastr = ri;		
	}

	// 
	/*for (ri = tid * nref; ri < (tid == num_threads - 1) ? ref_pos.size() : (tid + 1) * nref; ri++) {
		std::sort(nous[ri].begin(), nous[ri].end(), pair_compare);
		uint64_t k = 0;
		insert32(nous[ri][k].first, nous[ri][k].second - nous[ri][k].first, refID[ri]);
		for (size_t j = 1; j < nous[ri].size(); j++) {
			if (nous[ri][j].first >= nous[ri][k].second) {
				insert32(nous[ri][j].first, nous[ri][j].second - nous[ri][j].first, refID[ri]);
				k = j;
			}
		}
	}*/

	return 0;
}	

void* FastaReader::computeIndexmin() {
	int tid = 0;
	pthread_mutex_lock(&thread_lock);
	tid = current_tid++;
	pthread_mutex_unlock(&thread_lock);
	fprintf(stderr, "Begin computing covering substrings by thread %d, pid %lu.\n", tid + 1, pthread_self());

	uint64_t i = 1, j = 0, nexti = 1, start = 0;
	size_t nref = ref_pos.size() / num_threads;
	if (tid > 0)
		i = ref_pos[tid * nref - 1];
	nexti = (tid == num_threads - 1) ? ref_pos[ref_pos.size() - 1] : ref_pos[(tid + 1) * nref - 1];
	size_t ci = 0, ri = tid * nref, lastr = ri;
	for (; i >= contig_pos[ci]; ci++);
	uint64_t start_ = 0, lastj = 0, lastl = 0; // For computing covering substrings.

	while (i < nexti) {
		// Index not assigned.
		j = mu_index[i];
		if (j == 0) {
			i++;
			continue;
		}

		// Reached a separation region between two contigs.
		if (i >= contig_pos[ci] - 4) {
			if ((start + L + 2 >= contig_pos[ci]) && exist_unique_[ci]) {
				if (ri == lastr)
					uLmcount[ri] -= (start + L + 3 - contig_pos[ci]);
				else
					uLmcount[lastr] -= (start + L + 3 - contig_pos[ci]);
			}
			start = std::max(contig_pos[ci], i - L);
			ci++;
			if (i >= ref_pos[ri] - 4)
				ri++;
			if (start + L + 2 >= contig_pos[ci]) 
				exist_unique_[ci] = false;
		}

		// Substring spans two different contigs.
		if (ci > 0 && j - 1 < contig_pos[ci - 1]) {
			i++;
			continue;
		}

		// Substring have special characters "RYKMSWBDHVN-".
		bool special_character_ = false;
		for (uint64_t spos = j - 1; spos < i; spos++)
			if ((seqs_[spos] - base_offset != 65) && 
			    (seqs_[spos] - base_offset != 67) &&
			    (seqs_[spos] - base_offset != 71) &&
			    (seqs_[spos] - base_offset != 84)) {
				special_character_ = true;
				break;
			}
		if (special_character_) {
			i++;
			continue;
		}

		// Substring length exceeds L, is not considered within the index.
		uint64_t length = i - j + 1;
		if (length > L) {
			i++;
			continue;
		}

		// Process next reference species.
		/*if (j > ref_pos[ri])
			ri++;

		// Process next contig.
		if (j > contig_pos[ci]) {
			if ((start + L + 2 >= contig_pos[ci]) && exist_unique_[ci]) {
				if (ri == lastr)
					uLmcount[ri] -= (start + L + 3 - contig_pos[ci]);
				else
					uLmcount[lastr] -= (start + L + 3 - contig_pos[ci]);
			}
			start = std::max(contig_pos[ci], i - L);
			ci++;
			if (start + L + 2 >= contig_pos[ci]) {
				exist_unique_[ci] = false;
			}
		}*/

		// Index: minimum # unique substrings that cover every unique L-mer.
		if ((i > start_ + L) && (lastl > 0)) {
			insert32(lastj - 1, lastl, refID[lastr]);
			start_ = lastj;
		}

		// Aggregate unique L-mers.
		if (i <= start + L) {
			uLmcount[ri] += (j - start);
			start = j;
		} else {
			uLmcount[ri] += (j + L - i);
			start = j;
		}
		
		i++;
		lastr = ri;
		lastl = length;
		lastj = j;		
	}

	return 0;
}

//This has to be modified since some of the doubly unique substrings are not sampled!!!
//Modified on 07/21/2019
//Testing.
void* FastaReader::computeIndexmin_d() {
	int tid = 0;
	pthread_mutex_lock(&thread_lock);
	tid = current_tid++;
	pthread_mutex_unlock(&thread_lock);
	fprintf(stderr, "Begin computing covering substrings by thread %d, pid %lu.\n", tid + 1, pthread_self());

	uint64_t i = 1, j = 0, nexti = 1, start = 0;
	size_t nref = ref_pos.size() / num_threads;
	if (tid > 0)
		i = ref_pos[tid * nref - 1];
	nexti = (tid == num_threads - 1) ? ref_pos[ref_pos.size() - 1] : ref_pos[(tid + 1) * nref - 1];
	size_t ci = 0, ri = tid * nref, lastr = ri;
	for (; i >= contig_pos[ci]; ci++);
	uint64_t start_ = 0, lastj = 0, lastl = 0; // For computing covering substrings.

	while (i < nexti) {
		// Index not assigned.
		j = mu_index[i];
		if (j == 0) {
			i++;
			continue;
		}

		// Reached a separation region between two contigs.
		if (i >= contig_pos[ci] - 4) {
			if ((start + L + 2 >= contig_pos[ci]) && exist_unique_[ci]) {
				if (ri == lastr)
					uLmcount[ri] -= (start + L + 3 - contig_pos[ci]);
				else
					uLmcount[lastr] -= (start + L + 3 - contig_pos[ci]);
			}
			start = std::max(contig_pos[ci], i - L);
			ci++;
			if (i >= ref_pos[ri] - 4)
				ri++;
			if (start + L + 2 >= contig_pos[ci]) 
				exist_unique_[ci] = false;
		}

		// Substring spans two different contigs.
		if (ci > 0 && j - 1 < contig_pos[ci - 1]) {
			i++;
			continue;
		}

		// Substring have special characters "RYKMSWBDHVN-".
		bool special_character_ = false;
		for (uint64_t spos = j - 1; spos < i; spos++)
			if ((seqs_[spos] - base_offset != 65) && 
			    (seqs_[spos] - base_offset != 67) &&
			    (seqs_[spos] - base_offset != 71) &&
			    (seqs_[spos] - base_offset != 84)) {
				special_character_ = true;
				break;
			}
		if (special_character_) {
			i++;
			continue;
		}

		// Substring length exceeds L, is not considered within the index.
		uint64_t length = i - j + 1;
		if (length > L) {
			i++;
			continue;
		}

		// Process next reference species.
		/*if (j > ref_pos[ri])
			ri++;

		// Process next contig.
		if (j > contig_pos[ci]) {
			if ((start + L + 2 >= contig_pos[ci]) && exist_unique_[ci]) {
				if (ri == lastr)
					uLmcount[ri] -= (start + L + 3 - contig_pos[ci]);
				else
					uLmcount[lastr] -= (start + L + 3 - contig_pos[ci]);
			}
			start = std::max(contig_pos[ci], i - L);
			ci++;
			if (start + L + 2 >= contig_pos[ci]) {
				exist_unique_[ci] = false;
			}
		}*/

		// Index: minimum # unique substrings that cover every unique L-mer.
		if ((i > start_ + L) && (lastl > 0)) {
			insert32_d(lastj - 1, lastl, refID[lastr]);
			start_ = lastj;
		}

		// Aggregate unique L-mers.
		if (i <= start + L) {
			uLmcount[ri] += (j - start);
			start = j;
		} else {
			uLmcount[ri] += (j + L - i);
			start = j;
		}
		
		i++;
		lastr = ri;
		lastl = length;
		lastj = j;		
	}

	return 0;
}

void FastaReader::computeIndex(int mode) {
	mode_ = mode;
	switch (mode) {
		// 
		case 0:
			for (uint32_t i = 0; i < M_; i++) {
				minusl.push_back(0xFFFFul);
				maxusl.push_back(0);
				sumusl.push_back(0);
				uscount.push_back(0);
			}
			break;
		//
		case 1:
		case 2:
		case 4:
			for (uint32_t i = 0; i < M_; i++)
				uLmcount.push_back(0);
			for (uint32_t i = 0; i < C_; i++)
				exist_unique_.push_back(true);
			for (int i = 0; i < num_locks; i++) {
				pthread_spinlock_t *lock_i = new pthread_spinlock_t();
				hasht_access.push_back(lock_i);
				pthread_spin_init(lock_i, PTHREAD_PROCESS_SHARED);
			}
			hasht = new Hash(HASH_LEN_);

			break;
		default:
			for (uint32_t i = 0; i < M_; i++)
				uLmcount_clark.push_back(0);
			for (uint32_t i = 0; i < C_; i++)
				exist_unique_.push_back(true);
			break;
	}


	pthread_t threads[num_threads];
	pthread_mutex_init(&thread_lock, NULL);
	current_tid = 0;
	for (int i = 0; i < num_threads; i++)
		pthread_create(&threads[i], NULL, computeIndex_t, this);
	for (int i = 0; i < num_threads; i++)
		pthread_join(threads[i], NULL);
	

	if (mode == 0) {
		//fprintf(stderr, "Completed the output of unique substrings.\n");
		FILE *pFile = fopen("./unique_stats.out", "w");
		fprintf(pFile, "REFID\tMIN\tMAX\tAVG\tCNT\n");
		for (uint32_t i = 0; i < M_; i++)
			fprintf(pFile, "%u\t%lu\t%lu\t%.2f\t%u\n", refID[i],
				minusl[i], maxusl[i], sumusl[i] * 1.0 / uscount[i], uscount[i]);
		fclose (pFile);
	}
	if (mode == 1 || mode == 2) {
		std::string index_file = "index.bin";
		hasht->encodeIdx32(index_file);
		FILE *pFile = fopen("./unique_lmer_count.out", "w");
		fprintf(pFile, "REFID\tCNT\n");
		for (uint32_t i = 0; i < M_; i++)
			fprintf(pFile, "%u\t%u\n", refID[i], uLmcount[i]);
		fclose (pFile);
	}
	if (mode == 3) {
		FILE *pFile = fopen("./unique_klmer_count.out", "w");
		fprintf(pFile, "REFID\tCNT\n");
		for (uint32_t i = 0; i < M_; i++)
			fprintf(pFile, "%u\t%u\n", refID[i], uLmcount_clark[i]);
		fclose (pFile);
	}
	if (mode == 4) {
		hasht->encodeIdx32_d();
		FILE *pFile = fopen("./unique_lmer_count.out", "w");
		fprintf(pFile, "REFID\tCNT\n");
		for (uint32_t i = 0; i < M_; i++)
			fprintf(pFile, "%u\t%u\n", refID[i], uLmcount[i]);
		fclose (pFile);
	}
}

void FastaReader::setNumLocks(int nl) {
	if (nl > 0xFFFF) {
		num_locks = 1 << 16;
		return;
	}
	if (nl <= 0xFF) {
		num_locks = 1 << 8;
		return;
	}
	num_locks = 1 << (31 - __builtin_clz(nl));
}

void FastaReader::setHashLength(uint32_t hl) {
	HASH_LEN_ = hl;
}

