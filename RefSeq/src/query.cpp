#include <dirent.h>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <unordered_map>
#include <set>
#include <algorithm>

#include "query.hpp"
#include "binaryio.hpp"
#include "hashtrie.hpp"

FqReader::FqReader(int d, int ho, uint32_t hl, std::string &idx_fn, std::string &map_fn) {
	doubly_unique = d;
	hash_option = ho;
	//fprintf(stderr, "%d\n", d);
	if (d == 0) {
		ht_u = new Hash(hl);
		hash_len_u = hl;
		IDXFILEU = idx_fn;
	} else {
		ht_d = new Hash(hl);
		hash_len_d = hl;
		IDXFILED = idx_fn;
	}
	MAPFILE = map_fn;
}

FqReader::FqReader(int ho, uint32_t hl_u, std::string &idx_fn_u, 
		uint32_t hl_d, std::string &idx_fn_d, std::string &map_fn) {
	doubly_unique = 2;
	hash_option = ho;
	ht_u = new Hash(hl_u);
	ht_d = new Hash(hl_d);
	hash_len_u = hl_u;
	hash_len_d = hl_d;
	IDXFILEU = idx_fn_u;
	IDXFILED = idx_fn_d;
	MAPFILE = map_fn;
}

void FqReader::resetReads() { 
 	if (!reads.empty()) {
		for (auto read : reads)
			if (read != NULL)
				delete []read;
	}
}

FqReader::~FqReader() {
	resetReads();
	if (ht_u != NULL)
		delete ht_u;
	if (ht_d != NULL)
		delete ht_d;
}

void FqReader::loadIdx() {
	switch (doubly_unique) {
		case 0:
			if (hash_option == 32)
				ht_u->loadIdx32(IDXFILEU);
			else
				ht_u->loadIdx64(IDXFILEU);
			break;
		case 1:
			if (hash_option == 32)
				ht_d->loadIdx32_d(IDXFILED);
			else
				ht_d->loadIdx64_d(IDXFILED);
			break;
		case 2:
			if (hash_option == 32) {
				ht_u->loadIdx32(IDXFILEU);
				ht_d->loadIdx32_d(IDXFILED);
			} else {
				ht_u->loadIdx64(IDXFILEU);
				ht_d->loadIdx64_d(IDXFILED);
			}
			break;
	}
		
	fprintf(stderr, "Loaded index file.\n");
	//ht->testbucket(6421955);
}

void FqReader::loadSmap() {
	std::string line, fn, sp;
	std::ifstream inputFile;

	/* Read in fasta file. */ 
	inputFile.open(MAPFILE);
	while (std::getline(inputFile, line)) {
		std::istringstream lines(line);
		lines >> fn;
		lines >> sp;
		read_cnts[stoi(sp)] = 0;
	}
	inputFile.close();
	fprintf(stderr, "Loaded species map file.\n");
}

void FqReader::getFqList(std::string &INDIR) {
	if (!qfilenames.empty()) {
		qfilenames.clear();
		fprintf(stderr, "Flushed existing file names.\n");
	}
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir(INDIR.c_str())) != NULL) {
		while ((ent = readdir(dir)) != NULL) {
			std::string filename = ent->d_name;
			if (filename.length() >= 7) {
				std::string fext = filename.substr(filename.find_last_of(".") + 1);
				if (fext == "fq" || fext == "fastq")
					qfilenames.push_back(INDIR + filename);
			}
		}
		closedir(dir);
	} else {
		/* could not open directory */
		fprintf(stderr, "Input directory not exists.\n");
		abort();
	}
}

void FqReader::queryFastq(std::vector<std::string> &qfilenames_) {
	switch (doubly_unique) {
		case 0:
			if (hash_option == 32)
				for (auto fq_file : qfilenames_) {
					readFastq(fq_file);
					query32();
					resetReads();
				}
			else
				for (auto fq_file : qfilenames_) {
					readFastq(fq_file);
					query64();
					resetReads();
				}
			break;
		case 1:
			if (hash_option == 32)
				for (auto fq_file : qfilenames_) {
					readFastq(fq_file);
					query32_d();
					resetReads();
				}
			else
				for (auto fq_file : qfilenames_) {
					readFastq(fq_file);
					query64_d();
					resetReads();
				}
			break;
		case 2:
			if (hash_option == 32)
				for (auto fq_file : qfilenames_) {
					readFastq(fq_file);
					query32_ud();
					resetReads();
				}
			else
				for (auto fq_file : qfilenames_) {
					readFastq(fq_file);
					//query64_ud();
					resetReads();
				}
			break;
	}

	fprintf(stderr, "-Number of unlabeled reads: %lu.\n", nundet);
	fprintf(stderr, "-Number of reads with conflict labels: %lu.\n", nconf);
}

void FqReader::queryAllFastq() {
	switch (doubly_unique) {
		case 0:
			if (hash_option == 32)
				for (auto fq_file : qfilenames) {
					readFastq(fq_file);
					query32();
					resetReads();
				}
			else
				for (auto fq_file : qfilenames) {
					readFastq(fq_file);
					query64();
					resetReads();
				}
			break;
		case 1:
			if (hash_option == 32)
				for (auto fq_file : qfilenames) {
					readFastq(fq_file);
					query32_d();
					resetReads();
				}
			else
				for (auto fq_file : qfilenames) {
					readFastq(fq_file);
					query64_d();
					resetReads();
				}
			break;
		case 2:
			if (hash_option == 32)
				for (auto fq_file : qfilenames) {
					readFastq(fq_file);
					query32_ud();
					resetReads();
				}
			else
				for (auto fq_file : qfilenames) {
					readFastq(fq_file);
					//query64_ud();
					resetReads();
				}
			break;
	}
}

void FqReader::readFastq(std::string &INFILE) { 
 	std::string bases;
	std::ifstream inputFile;
	size_t rl;

	/* Read in fasta file. */ 
	inputFile.open(INFILE);
	while (std::getline(inputFile, bases)) {
		std::getline(inputFile, bases);
		rl = bases.length();
		//max_rl = std::max(max_rl, rl);
		uint8_t *read = new uint8_t[rl];
		strncpy((char*)read, bases.c_str(), rl);
		//fprintf(stderr, "%lu, %lu.\n", sizeof(read), rl);
		reads.push_back(read);
		std::getline(inputFile, bases);
		std::getline(inputFile, bases);
	}
	
	inputFile.close();
	fprintf(stderr, "%lu\n", reads.size());
}

void FqReader::getRC(uint8_t *dest, uint8_t *srcs, size_t length) {
	for(size_t i = 0; i < length; i++) 
		dest[i] = rcIdx[srcs[length - i - 1]];
}

/* Query L-mers covered by unique substrings. */
/* Return: read counts from each species (stored in READ_CNTS). */
void FqReader::query32() {
	/* Use Rabin-Karp algorithm. */
	uint32_t hv = 0, rid = 0, tmp_rid = 0, hs = 2 * hash_len_u -2;
	size_t rl = 0;
	uint8_t *rc_read = new uint8_t[max_rl];
	bool conflict_flag = 0;
	size_t nrd = 0;
	for (auto read : reads) {
		rid = 0;
		conflict_flag = 0;
		rl = 100;
		//fprintf(stderr, "Processed %lu reads.%lu\n", nrd, rl);

		/* First hash. */
		hv = 0;
		for (size_t i = 0; i < hash_len_u; i++)
			hv = ((hv << 2) | symbolIdx[read[i]]);

		//for (int j = 0; j < 100; j++)
		//	fprintf(stderr, "%c", (char) read[j]);
		//fprintf(stderr, "\n");

		/* Forward strand. */		
		for (size_t i = 0; i <= rl - hash_len_u; i++) {
			/*if (nrd == 111985) {
				for (int j = 0; j < 100; j++)
					fprintf(stderr, "%c", (char) read[j]);
				fprintf(stderr, "\n%lu, %u\n", i, hv);
				tmp_rid = ht->find32_test(hv, read + i, rl - hash_len_ - i);
			}*/
			
			//fprintf(stderr, "%u; %lu; %lu\n", hv, hash_len_, i);
			tmp_rid = ht_u->find32(hv, read + i + hash_len_u, rl - hash_len_u - i);

			if (tmp_rid != 0) {
				if (rid == 0)
					rid = tmp_rid;
				else if (tmp_rid != rid) {
					fprintf(stderr, "%d\n", symbolIdx[read[76]]);
					fprintf(stderr, "\n%lu, %u, %u\n", nrd, rid, tmp_rid);
					abort();
					conflict_flag = 1;
					break;
				}
			}
			/* Next hash. */
			if (i < rl - hash_len_u) {
			hv = hv - (symbolIdx[read[i]] << hs);
			hv = ((hv << 2) | symbolIdx[read[i + hash_len_u]]);}
		}
		if (conflict_flag) {
			nconf++;
			continue;
		}
		//fprintf(stderr, "Processed %lu reads.\n", nrd);

		/* Reverse complement. */
		getRC(rc_read, read, rl);
		hv = 0;
		for (size_t i = 0; i < hash_len_u; i++)
			hv = ((hv << 2) | symbolIdx[rc_read[i]]);
		for (size_t i = 0; i <= rl - hash_len_u; i++) {
			/*if (nrd == 352412) {
				for (int j = 0; j < 100; j++)
					fprintf(stderr, "%c", (char) rc_read[j]);
				fprintf(stderr, "\n%lu, %u\n", i, hv);
				uint32_t hv_ = 0;
				for (size_t j = 0; j < hash_len_; j++)
					hv_ = ((hv_ << 2) | symbolIdx[(rc_read + i)[j]]);
				fprintf(stderr, "%u\n", hv_);
			}*/
			tmp_rid = ht_u->find32(hv, rc_read + i + hash_len_u, rl - hash_len_u - i);
			if (tmp_rid != 0) {
				if (rid == 0)
					rid = tmp_rid;
				else if (tmp_rid != rid) {
					conflict_flag = 1;
					break;
				}
			}
			//Next hash. 
			if (i < rl - hash_len_u) {
			hv = hv - (symbolIdx[rc_read[i]] << hs);
			hv = ((hv << 2) | symbolIdx[rc_read[i + hash_len_u]]);}
		}
		
		/* Record results. */
		if (conflict_flag == 1)
			nconf++;
		else {
			if (rid == 0)
				nundet++;
			else
				read_cnts[rid]++;
		}
		if (nrd++ % 100000 == 0)
			fprintf(stderr, "Processed %lu reads.\n", nrd);
	}
	fprintf(stderr, "Number of unlabeled reads: %lu.\n", nundet);
	fprintf(stderr, "Number of reads with conflict labels: %lu.\n", nconf);

	delete[] rc_read;
}

void FqReader::query64() {
}

/* Query L-mers covered by doubly-unique substrings. */
void FqReader::query32_d() {
	/* Use Rabin-Karp algorithm. */
	uint32_t hv = 0, hs = 2 * hash_len_d - 2;
	std::pair<uint32_t, uint32_t> tmp_rid;
	std::set<std::pair<uint32_t, uint32_t>> rid_pairs;
	size_t rl = 0;
	uint8_t *rc_read = new uint8_t[max_rl];
	size_t nrd = 0;
	for (auto read : reads) {
		rid_pairs.clear();
		rl = 100;
		//conflict_flag = 0;

		/* First hash. */
		hv = 0;
		for (size_t i = 0; i < hash_len_d; i++)
			hv = ((hv << 2) | symbolIdx[read[i]]);

		/* Forward strand. */	
		for (size_t i = 0; i <= rl - hash_len_d; i++) {
			tmp_rid = ht_d->find32_d(hv, read + i + hash_len_d, rl - hash_len_d - i);
			/*if (tmp_rid.first != 0)	{
				fprintf(stderr, "%lu; %u\n", i, tmp_rid.first);
				rids.insert(tmp_rid.first);
			}
			if (tmp_rid.second != 0) {
				fprintf(stderr, "%lu; %u\n", i, tmp_rid.second);
				rids.insert(tmp_rid.second);
			}*/
			if (tmp_rid.first != 0 || tmp_rid.second != 0) {
				rid_pairs.insert(tmp_rid);
			}
			/*if (rids.size() > 2) {
				conflict_flag = 1;
				fprintf(stderr, "%u\t%lu\n", hash_len_, nrd);
				for (int j = 0; j < 100; j++)
					fprintf(stderr, "%c", (char) read[j]);
				fprintf(stderr, "\n");

				abort();
				//break;
			}*/
			/* Next hash. */
			if (i < rl - hash_len_d) {
				hv = hv - (symbolIdx[read[i]] << hs);
				hv = ((hv << 2) | symbolIdx[read[i + hash_len_d]]);
			}
		}

		/* Reverse complement. */
		getRC(rc_read, read, rl);
		hv = 0;
		for (size_t i = 0; i < hash_len_d; i++)
			hv = ((hv << 2) | symbolIdx[rc_read[i]]);
		//GACAACGATTACGGCTTCCAATATACTCCGATTCTTGAAACAAATGCCAAAGGTGTATGGATTGAAAAAGAACAAACAGATTTACAGGAATCACCAGTAG
		for (size_t i = 0; i <= rl - hash_len_d; i++) {
			tmp_rid = ht_d->find32_d(hv, rc_read + i + hash_len_d, rl - hash_len_d - i);
			/*if (tmp_rid.first != 0) {
				fprintf(stderr, "%lu; %u, %u\n", i, tmp_rid.first, hv);
				rids.insert(tmp_rid.first);
			}
			if (tmp_rid.second != 0) {
				fprintf(stderr, "%lu; %u, %u\n", i, tmp_rid.second, hv);
				rid_pairs.insert(tmp_rid.second);
			}
			if (rids.size() > 2) {
				fprintf(stderr, "%u\t%lu\n", hash_len_, nrd);
				for (int j = 0; j < 100; j++)
					fprintf(stderr, "%c", (char) rc_read[j]);
				fprintf(stderr, "-RC\n");

				conflict_flag = 1;
				abort();
				//break;
			}*/
			if (nrd == 3 && tmp_rid.first != 0) {
				for (int j = 0; j < 100; j++)
					fprintf(stderr, "%c", (char) rc_read[j]);
				fprintf(stderr, "-RC\n");
				fprintf(stderr, "%u, %u, %u", tmp_rid.first, temp_rid.second, hash_len_);
			}
			if (tmp_rid.first != 0 || tmp_rid.second != 0) {
				rid_pairs.insert(tmp_rid);
			}
			/* Next hash. */
			if (i < rl - hash_len_d) {
				hv = hv - (symbolIdx[rc_read[i]] << hs);
				hv = ((hv << 2) | symbolIdx[rc_read[i + hash_len_d]]);
			}
		}

		/* Record results. */
		switch (rid_pairs.size()) {
			case 0:
				nundet++;
				break;
			case 1:
				for (auto rid_pair : rid_pairs) {
					if (rid_pair.first != 0)
						read_cnts[rid_pair.first]++;
					if (rid_pair.second != 0)
						read_cnts[rid_pair.second]++;
				}
				break;
			default:
				std::set<uint32_t> intersection;
				int i = 0;
				for (auto rid_pair : rid_pairs) {
					if (nrd < 20)
					fprintf(stderr, "(%u, %u)\t", rid_pair.first, rid_pair.second);
					if (i == 0) {
						/*if (rid_pair.first != 0)
							intersection.insert(rid_pair.first);
						if (rid_pair.second != 0)
							intersection.insert(rid_pair.second);*/
						intersection.insert(rid_pair.first);
						intersection.insert(rid_pair.second);
						intersection.erase(0);
					} else {
						for (auto rid : intersection) {
							if (rid_pair.first != rid && rid_pair.second != rid)
								intersection.erase(rid);
						}
						/*intersection.erase(
							std::remove_if(
								intersection.begin(), 
								intersection.end(),
								[rid_pair] (uint32_t rid) { 
									//return (rid_pair.first != rid && 
									//	rid_pair.second != rid); 
									return true;
								}
							), 
							intersection.end()
						);*/
					}
				}
				if (nrd < 20)
				fprintf(stderr, "%u\t%lu\n", nrd, intersection.size());
				switch (intersection.size()) {
					case 0: 
						nundet++;
						break;
					case 1:
						for (auto rid : intersection)
							read_cnts[rid]++;
						break;
					default:
						nconf++;
						break;
				}
				break;
		}
		if (nrd++ % 100000 == 0) {
			fprintf(stderr, "Processed %lu reads.\n", nrd);
			//break;
		}
	}

	for (auto rc : read_cnts)
	{
		fprintf(stderr, "%u\t%u\n", rc.first, rc.second);
	}
	fprintf(stderr, "Number of unlabeled reads: %lu.\n", nundet);
	fprintf(stderr, "Number of reads with conflict labels: %lu.\n", nconf);
	//delete []rc_read;
	fprintf(stderr, "BB\n");
}

void FqReader::query64_d() {
}


void FqReader::query32_ud() {
	uint32_t hv = 0, hs = 0, rid = 0;
	std::pair<uint32_t, uint32_t> tmp_rid;
	std::set<uint32_t> rids;
	std::set<std::pair<uint32_t, uint32_t>> rid_pairs;
	size_t rl = 0;
	uint8_t *rc_read = new uint8_t[max_rl];
	size_t nrd = 0;
	for (auto read : reads) {
		rid_pairs.clear();
		rl = 100;

		/* Forward strand. */
		/* Query unique substrings. */
		hs = 2 * hash_len_u - 2;
		hv = 0;
		for (size_t i = 0; i < hash_len_u; i++)
			hv = ((hv << 2) | symbolIdx[read[i]]);

		for (size_t i = 0; i <= rl - hash_len_u; i++) {
			rid = ht_u->find32_u(hv, read + i + hash_len_u, rl - hash_len_u - i);
			
			if (rid != 0)
				rids.insert(rid);
			
			/* Next hash. */
			if (i < rl - hash_len_u) {
				hv = hv - (symbolIdx[read[i]] << hs);
				hv = ((hv << 2) | symbolIdx[read[i + hash_len_u]]);
			}
		}

		if (rids.size() > 1) {
			nconf++;
			continue;
		}

		/* Query unique substrings. */
		hs = 2 * hash_len_u - 2;
		hv = 0;
		for (size_t i = 0; i < hash_len_u; i++)
			hv = ((hv << 2) | symbolIdx[read[i]]);

		for (size_t i = 0; i <= rl - hash_len_u; i++) {
			rid = ht_u->find32_u(hv, read + i + hash_len_u, rl - hash_len_u - i);
			
			if (rid != 0)
				rids.insert(rid);
			
			/* Next hash. */
			if (i < rl - hash_len_u) {
				hv = hv - (symbolIdx[read[i]] << hs);
				hv = ((hv << 2) | symbolIdx[read[i + hash_len_u]]);
			}
		}

		if (rids.size() > 1) {
			nconf++;
			continue;
		}

		/* Reverse complement. */
		getRC(rc_read, read, rl);
		hv = 0;
		for (size_t i = 0; i < hash_len_d; i++)
			hv = ((hv << 2) | symbolIdx[rc_read[i]]);
		//GACAACGATTACGGCTTCCAATATACTCCGATTCTTGAAACAAATGCCAAAGGTGTATGGATTGAAAAAGAACAAACAGATTTACAGGAATCACCAGTAG
		for (size_t i = 0; i <= rl - hash_len_d; i++) {
			tmp_rid = ht_d->find32_d(hv, rc_read + i + hash_len_d, rl - hash_len_d - i);
			/*if (tmp_rid.first != 0) {
				fprintf(stderr, "%lu; %u, %u\n", i, tmp_rid.first, hv);
				rids.insert(tmp_rid.first);
			}
			if (tmp_rid.second != 0) {
				fprintf(stderr, "%lu; %u, %u\n", i, tmp_rid.second, hv);
				rid_pairs.insert(tmp_rid.second);
			}
			if (rids.size() > 2) {
				fprintf(stderr, "%u\t%lu\n", hash_len_, nrd);
				for (int j = 0; j < 100; j++)
					fprintf(stderr, "%c", (char) rc_read[j]);
				fprintf(stderr, "-RC\n");

				conflict_flag = 1;
				abort();
				//break;
			}*/
			if (nrd == 3 && tmp_rid.first != 0) {
				for (int j = 0; j < 100; j++)
					fprintf(stderr, "%c", (char) rc_read[j]);
				fprintf(stderr, "-RC\n");
				fprintf(stderr, "%u, %u, %u", tmp_rid.first, temp_rid.second, hash_len_);
			}
			if (tmp_rid.first != 0 || tmp_rid.second != 0) {
				rid_pairs.insert(tmp_rid);
			}
			/* Next hash. */
			if (i < rl - hash_len_d) {
				hv = hv - (symbolIdx[rc_read[i]] << hs);
				hv = ((hv << 2) | symbolIdx[rc_read[i + hash_len_d]]);
			}
		}

		/* Record results. */
		switch (rid_pairs.size()) {
			case 0:
				nundet++;
				break;
			case 1:
				for (auto rid_pair : rid_pairs) {
					if (rid_pair.first != 0)
						read_cnts[rid_pair.first]++;
					if (rid_pair.second != 0)
						read_cnts[rid_pair.second]++;
				}
				break;
			default:
				std::set<uint32_t> intersection;
				int i = 0;
				for (auto rid_pair : rid_pairs) {
					if (nrd < 20)
					fprintf(stderr, "(%u, %u)\t", rid_pair.first, rid_pair.second);
					if (i == 0) {
						/*if (rid_pair.first != 0)
							intersection.insert(rid_pair.first);
						if (rid_pair.second != 0)
							intersection.insert(rid_pair.second);*/
						intersection.insert(rid_pair.first);
						intersection.insert(rid_pair.second);
						intersection.erase(0);
					} else {
						for (auto rid : intersection) {
							if (rid_pair.first != rid && rid_pair.second != rid)
								intersection.erase(rid);
						}
						/*intersection.erase(
							std::remove_if(
								intersection.begin(), 
								intersection.end(),
								[rid_pair] (uint32_t rid) { 
									//return (rid_pair.first != rid && 
									//	rid_pair.second != rid); 
									return true;
								}
							), 
							intersection.end()
						);*/
					}
				}
				if (nrd < 20)
				fprintf(stderr, "%u\t%lu\n", nrd, intersection.size());
				switch (intersection.size()) {
					case 0: 
						nundet++;
						break;
					case 1:
						for (auto rid : intersection)
							read_cnts[rid]++;
						break;
					default:
						nconf++;
						break;
				}
				break;
		}
		if (nrd++ % 100000 == 0) {
			fprintf(stderr, "Processed %lu reads.\n", nrd);
			//break;
		}
	}

	for (auto rc : read_cnts)
	{
		fprintf(stderr, "%u\t%u\n", rc.first, rc.second);
	}
	fprintf(stderr, "Number of unlabeled reads: %lu.\n", nundet);
	fprintf(stderr, "Number of reads with conflict labels: %lu.\n", nconf);
	//delete []rc_read;
	fprintf(stderr, "BB\n");
}

/*
// The following is from HARC
void readDnaFile(std::bitset<2*readlen> *read) 
{ 
	#pragma omp parallel 
	{ 
	int tid = omp_get_thread_num(); 
	uint32_t i, stop;	 
	//doing initial setup and first read 
	i = uint64_t(tid)*numreads/omp_get_num_threads();//spread out first read equally 
	stop = uint64_t(tid+1)*numreads/omp_get_num_threads(); 
	if(tid == omp_get_num_threads()-1) 
		stop = numreads; 
	std::ifstream f(infile, std::ifstream::in); 
	f.seekg(uint64_t(i)*(readlen+1), f.beg); 
	std::string s; 
	while(i < stop) 
	{ 
		std::getline(f,s); 
		read[i] = stringtobitset(s); 
		i++; 
	} 
	f.close(); 
	} 
	return; 
}*/

int FqReader::symbolIdx[256] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, 1, -1, -1, \
-1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, -1, -1, -1, -1, -1, \
-1, -1, -1, -1, -1, -1, -1, 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, \
-1, -1, -1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, \
-1, -1, -1, -1, -1, -1};

int FqReader::rcIdx[128] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 84, -1, 71, -1, -1, \
-1, 67, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 65, -1, -1, -1, -1, -1, \
-1, -1, -1, -1, -1, -1, -1, 84, -1, 71, -1, -1, -1, 67, -1, -1, -1, -1, -1, -1, \
-1, -1, -1, -1, -1, -1, 65, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

