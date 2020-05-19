#include <dirent.h>
#include <pthread.h>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <algorithm>

#include "util.hpp"
#include "build.hpp"

int FastaReader::mode_;

FastaReader::FastaReader(int L_, int K, int t) {
	L = L_;
	minuL = K;
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
	minuL = K;
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
	fprintf(stderr, "Total num genomes: %u\n", M_);
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
		int sp_ = stoi(sp);
		filenames[INDIR + fn] = sp_;
		if (sp_ > max_rid)
			max_rid = sp_;
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
				} else {
					if (!insertRC()) {
						inputFile.close();
						abort();
					}
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
	} else {
		if (!insertRC()) {
			inputFile.close();
			abort();
		}
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
	size_t l = bases.length();
	size_t i = 0, j = 0;
	if (N_ + l >= maxN) {
		fprintf(stderr, "Total number of symbols exceeds limit.\n");
		return false;
	}
	N_ += l;
	while (pos + l > block_size_) {
		j = block_size_ - pos;
		for (size_t j_ = 0; j_ < j; j_++, i++)
			seqs[cur][pos++] = bases[i] + base_offset;
		allocSeq();
		l -= j;
	}
	for (size_t i_ = 0; i_ < l; i_++, i++)
		seqs[cur][pos++] = bases[i] + base_offset;
	rc_contig += bases;
	return true;
}

bool FastaReader::insertBaseRC(std::string& rc_bases) {
	size_t l = rc_bases.length();
	size_t i = 0, j = 0;
	if (N_ + l >= maxN) {
		fprintf(stderr, "Total number of symbols exceeds limit.\n");
		return false;
	}
	N_ += l;
	while (pos + l > block_size_) {
		j = block_size_ - pos;
		for (size_t j_ = 0; j_ < j; j_++, i++)
			seqs[cur][pos++] = rc_bases[i] + base_offset;
		allocSeq();
		l -= j;
	} 
	for (size_t i_ = 0; i_ < l; i_++, i++)
		seqs[cur][pos++] = rc_bases[i] + base_offset;
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

	C_++;
	if (C_ >= maxC) {
		fprintf(stderr, "Number of contigs exceeds limit.\n");
		return false;
	}
	return true;
}

char complement(char c) {   
	switch(c) {   
		case 'A':
			return 'T';
		case 'T':
			return 'A';
		case 'G':
			return 'C';
		case 'C':
			return 'G';
		default:
			return 'N';
	}
	return 'N';
}

bool FastaReader::insertRC() {
	std::reverse(rc_contig.begin(), rc_contig.end());
	std::transform(std::begin(rc_contig), std::end(rc_contig), std::begin(rc_contig), complement);
	if (!insertBaseRC(rc_contig))
		return false;
	if (!insertContig())
		return false;
	rc_contig.clear();
	return true;
}

void FastaReader::allocSuffixArray(int doubly_unique) {
	uint64_t N = cur * block_size_ + pos;
	switch (doubly_unique) {
		case 0:
			/* Allocate the memory space of input string. */
			seqs_ = new uint8_t[N + 2];
			for (size_t i = 0; i < cur; i++)
				memcpy(seqs_ + i * block_size_, seqs[i], block_size_);
			memcpy(seqs_ + cur * block_size_, seqs[cur], pos);
			fprintf(stderr, "Input string get prepared.\n");

			/* Compute suffix arrays. */	
			seqs_[N] = 0;
			seqs_[N + 1] = 0;
			if (max_rid > 0xFFFE)
				sa = new SuffixArray(N, 32);
			else
				sa = new SuffixArray(N, 16);
			sa->run(seqs_, ref_pos, refID, N, 0, L, 0, 1, 0, 0);
			return;
		case 1:
			mu_index = sa->run(seqs_, ref_pos, refID, N, 1, L, minuL - 1, 1, 0, 0);
			occ = sa->occ;
			return;
		case 2:
			mu_index = sa->run(seqs_, ref_pos, refID, N, 2, L, minuL - 1, 1, 0, 0);
			occ = sa->occ;
			occ2 = sa->occ2;
			if (max_rid > 0xFFFE)
				GSA2 = sa->GSA32_;
			else
				GSA2_ = sa->GSA16_;
			return;
		case 3:
			mu_index = sa->run(seqs_, ref_pos, refID, N, 1, L, minuL - 1, 1, 0, 1);
			occ = sa->occ;
			return;
	}	
}

void FastaReader::insert32(uint64_t i, uint32_t length, uint32_t rid, uint8_t occ) {
	uint32_t bucket = hasht->computeHashVal(seqs_ + i);
	pthread_spin_lock(hasht_access);
	hasht->insert32(bucket, seqs_ + i, length, rid, occ);
	pthread_spin_unlock(hasht_access);
}

void FastaReader::insert32_d(uint64_t i, uint32_t length, uint32_t rid, uint32_t rid2, uint8_t occ, uint8_t occ2) {
	uint32_t bucket = hasht->computeHashVal(seqs_ + i);
	pthread_spin_lock(hasht_access);
	hasht->insert32_d(bucket, seqs_ + i, length, rid, occ, rid2, occ2);
	pthread_spin_unlock(hasht_access);
}

void FastaReader::insert64(uint64_t i, uint32_t length, uint32_t rid, uint8_t occ) {
	uint64_t bucket = hasht->computeHashVal64(seqs_ + i);
	pthread_spin_lock(hasht_access);
	hasht->insert64(bucket, seqs_ + i, length, rid, occ);
	pthread_spin_unlock(hasht_access);
}

void FastaReader::insert64_d(uint64_t i, uint32_t length, uint32_t rid, uint32_t rid2, uint8_t occ, uint8_t occ2) {
	uint64_t bucket = hasht->computeHashVal64(seqs_ + i);
	pthread_spin_lock(hasht_access);
	hasht->insert64_d(bucket, seqs_ + i, length, rid, occ, rid2, occ2);
	pthread_spin_unlock(hasht_access);
}	

void* FastaReader::computeIndexmin() {
	int tid = 0;
	pthread_mutex_lock(&thread_lock);
	tid = current_tid++;
	pthread_mutex_unlock(&thread_lock);
	//fprintf(stderr, "L: %d; maxuL: %d\n", L, maxuL);
	//fprintf(stderr, "Begin computing unique substrings by thread %d.\n", tid + 1);

	uint64_t i = 1, j = 0, nexti = 1, start = 0;
	size_t nref = ref_pos.size() / num_threads;
	if (tid > 0)
		i = ref_pos[tid * nref - 1];
	nexti = (tid == num_threads - 1) ? ref_pos[ref_pos.size() - 1] : ref_pos[(tid + 1) * nref - 1];
	size_t ci = 0, ri = tid * nref, lastr = ri;
	for (; i >= contig_pos[ci]; ci++);
	uint64_t start_ = 0, lastj = 0, lastl = 0; // For computing a minimum sample of substrings.

	while (i < nexti) {
		// Index not assigned.
		if (mu_index[i] == UINT16_MAX) {
			i++;
			continue;
		}
		j = i - mu_index[i];

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
			if (ci >= contig_pos.size())
				break;
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

		// Substring length exceeds maxuL, is not considered within the index.
		uint64_t length = i - j + 1;
		if (length > maxuL) {
			i++;
			continue;
		}

		// Index: minimum # unique substrings that cover every unique L-mer.
		if ((i > start_ + L) && (lastl > 0)) {
			if (HASH_LEN_ < 16)
				insert32(lastj - 1, lastl, refID[lastr], occ[lastj - 1]);
			else
				insert64(lastj - 1, lastl, refID[lastr], occ[lastj - 1]);
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

void* FastaReader::computeIndexmin_d() {
	int tid = 0;
	pthread_mutex_lock(&thread_lock);
	tid = current_tid++;
	pthread_mutex_unlock(&thread_lock);
	//fprintf(stderr, "L: %d; maxuL: %d\n", L, maxuL);
	//fprintf(stderr, "Begin computing doubly-unique substrings by thread %d.\n", tid + 1);

	uint64_t i = 1, j = 0, nexti = 1, start = 0;
	size_t nref = ref_pos.size() / num_threads;
	if (tid > 0)
		i = ref_pos[tid * nref - 1];
	nexti = (tid == num_threads - 1) ? ref_pos[ref_pos.size() - 1] : ref_pos[(tid + 1) * nref - 1];
	size_t ci = 0, ri = tid * nref, lastr = ri;
	for (; i >= contig_pos[ci]; ci++);
	uint64_t start_ = 0, lastj = 0, lastl = 0; // For computing a minimum sample of substrings.

	while (i < nexti) {
		// Index not assigned.
		if (mu_index[i] == UINT16_MAX) {
			i++;
			continue;
		}
		j = i - mu_index[i];

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
			if (ci >= contig_pos.size())
				break;
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

		// Substring length exceeds maxuL, is not considered within the index.
		uint64_t length = i - j + 1;
		if (length > maxuL) {
			i++;
			continue;
		}

		// Index: minimum # unique substrings that cover every unique L-mer.
		if ((i > start_ + L) && (lastl > 0)) {
			if (HASH_LEN_ < 16)
				insert32_d(lastj - 1, lastl, refID[lastr], GSA2[lastj - 1], occ[lastj - 1], occ2[lastj - 1]);
			else
				insert64_d(lastj - 1, lastl, refID[lastr], GSA2[lastj - 1], occ[lastj - 1], occ2[lastj - 1]);
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

void* FastaReader::computeIndexmin_d_() {
	int tid = 0;
	pthread_mutex_lock(&thread_lock);
	tid = current_tid++;
	pthread_mutex_unlock(&thread_lock);
	//fprintf(stderr, "L: %d; maxuL: %d\n", L, maxuL);
	//fprintf(stderr, "Begin computing doubly-unique substrings by thread %d.\n", tid + 1);

	uint64_t i = 1, j = 0, nexti = 1, start = 0;
	size_t nref = ref_pos.size() / num_threads;
	if (tid > 0)
		i = ref_pos[tid * nref - 1];
	nexti = (tid == num_threads - 1) ? ref_pos[ref_pos.size() - 1] : ref_pos[(tid + 1) * nref - 1];
	size_t ci = 0, ri = tid * nref, lastr = ri;
	for (; i >= contig_pos[ci]; ci++);
	uint64_t start_ = 0, lastj = 0, lastl = 0; // For computing a minimum sample of substrings.

	while (i < nexti) {
		// Index not assigned.
		if (mu_index[i] == UINT16_MAX) {
			i++;
			continue;
		}
		j = i - mu_index[i];

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
			if (ci >= contig_pos.size())
				break;
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

		// Substring length exceeds maxuL, is not considered within the index.
		uint64_t length = i - j + 1;
		if (length > maxuL) {
			i++;
			continue;
		}

		// Index: minimum # unique substrings that cover every unique L-mer.
		if ((i > start_ + L) && (lastl > 0)) {
			if (HASH_LEN_ < 16)
				insert32_d(lastj - 1, lastl, refID[lastr], GSA2_[lastj - 1], occ[lastj - 1], occ2[lastj - 1]);
			else
				insert64_d(lastj - 1, lastl, refID[lastr], GSA2_[lastj - 1], occ[lastj - 1], occ2[lastj - 1]);
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
	if ((max_rid > 0xFFFE) || (mode == 1))
		mode_ = mode;
	else
		mode_ = mode + 2;
	switch (mode) {
		case 1:
		case 2:
			for (uint32_t i = 0; i < M_; i++)
				uLmcount.push_back(0);
			for (uint32_t i = 0; i < C_; i++)
				exist_unique_.push_back(true);
			hasht = new Hash(HASH_LEN_);
			hasht_access = new pthread_spinlock_t();
			pthread_spin_init(hasht_access, PTHREAD_PROCESS_SHARED);

			break;
		case 3:
			for (uint32_t i = 0; i < M_; i++)
				uLmcount[i] = 0;
			for (uint32_t i = 0; i < C_; i++)
				exist_unique_[i] = true;
			hasht->clear();
			pthread_spin_unlock(hasht_access);
		default:
			break;
	}

	auto start = std::chrono::high_resolution_clock::now();
	pthread_t threads[num_threads];
	current_tid = 0;
	for (int i = 0; i < num_threads; i++)
		pthread_create(&threads[i], NULL, computeIndex_t, this);
	for (int i = 0; i < num_threads; i++)
		pthread_join(threads[i], NULL);
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
					(std::chrono::high_resolution_clock::now() - start).count();
	fprintf(stderr, "Time for organizing index: %lu ms.\n", duration);

	if (mode == 1) {
		std::string index_file = "index_u.bin1";
		if (HASH_LEN_ < 16)
			hasht->encodeIdx32(index_file, 0);
		else
			hasht->encodeIdx64(index_file, 0);
		FILE *lcFile = fopen("./unique_lmer_count_u.out", "w");
		fprintf(lcFile, "REFID\tCNT\n");
		for (uint32_t i = 0; i < M_; i++)
			fprintf(lcFile, "%u\t%u\n", refID[i], uLmcount[i]);
		fclose (lcFile);
		FILE *glFile = fopen("./genome_lengths.out", "w");
		fprintf(glFile, "REFID\tLENGTH\n");
		uint32_t gl = 0, j = 0;
		for (uint32_t i = 0; i < C_; i++) {
			gl += (contig_pos[i] - ((i == 0) ? 0 : contig_pos[i - 1]) - 4);
			if (contig_pos[i] >= ref_pos[j]) {
				fprintf(glFile, "%u\t%u\n", refID[j], gl / 2);
				j++;
				gl = 0;
			}	
		}
		fclose (glFile);
	}
	if (mode == 2) {
		std::string index_file = "index_d.bin2";
		if (HASH_LEN_ < 16)
			hasht->encodeIdx32_d(index_file, 0);
		else
			hasht->encodeIdx64_d(index_file, 0);
		FILE *lcFile = fopen("./unique_lmer_count_d.out", "w");
		fprintf(lcFile, "REFID\tCNT\n");
		for (uint32_t i = 0; i < M_; i++)
			fprintf(lcFile, "%u\t%u\n", refID[i], uLmcount[i]);
		fclose (lcFile);
		FILE *glFile = fopen("./genome_lengths.out", "w");
		fprintf(glFile, "REFID\tLENGTH\n");
		uint32_t gl = 0, j = 0;
		for (uint32_t i = 0; i < C_; i++) {
			gl += (contig_pos[i] - ((i == 0) ? 0 : contig_pos[i - 1]) - 4);
			if (contig_pos[i] >= ref_pos[j]) {
				fprintf(glFile, "%u\t%u\n", refID[j], gl / 2);
				j++;
				gl = 0;
			}
		}
		fclose (glFile);
	}
	if (mode == 3) {
		std::string index_file = "index_d.bin2";
		if (HASH_LEN_ < 16)
			hasht->encodeIdx32_d(index_file, 0);
		else
			hasht->encodeIdx64_d(index_file, 0);
		FILE *pFile = fopen("./unique_lmer_count_d.out", "w");
		fprintf(pFile, "REFID\tCNT\n");
		for (uint32_t i = 0; i < M_; i++)
			fprintf(pFile, "%u\t%u\n", refID[i], uLmcount[i]);
		fclose (pFile);
	}
}

void FastaReader::setHashLength(uint32_t hl) {
	HASH_LEN_ = hl;
}

