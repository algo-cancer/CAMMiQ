#include <dirent.h>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <chrono>
#include <ilcplex/ilocplex.h>

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

FqReader::FqReader(int ho, uint32_t hl_u, std::string &idx_fn_u, 
		uint32_t hl_d, std::string &idx_fn_d, std::string &map_fn,
		std::string &output_fn, float erate) {
	doubly_unique = 2;
	hash_option = ho;
	ht_u = new Hash(hl_u);
	ht_d = new Hash(hl_d);
	hash_len_u = hl_u;
	hash_len_d = hl_d;
	IDXFILEU = idx_fn_u;
	IDXFILED = idx_fn_d;
	MAPFILE = map_fn;
	OUTPUTFILE = output_fn;
	erate_ = erate;
}

FqReader::FqReader(int ho, uint32_t hl, std::string &idx_fn_u, 
		std::string &idx_fn_d, std::string &map_fn,
		std::string &output_fn, float erate) {
	doubly_unique = 2;
	hash_option = ho;
	ht_u = new Hash(hl);
	hash_len_u = hl;
	IDXFILEU = idx_fn_u;
	IDXFILED = idx_fn_d;
	MAPFILE = map_fn;
	OUTPUTFILE = output_fn;
	erate_ = erate;
}

void FqReader::resetReads() { 
 	if (!reads.empty()) {
		for (auto read : reads)
			if (read != NULL)
				delete []read;
	}
	reads.clear();
}

FqReader::~FqReader() {
	//resetReads();
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

void FqReader::loadIdx_p() {
	fprintf(stderr, "%s\n", IDXFILEU.c_str());
	//std::string line, sp, cnt;
	//if (ho == )
	//ht_u->loadIdx64_b(IDXFILEU, IDXFILED);
	//std::ifstream inputFile;
 
	/*inputFile.open("unique_substring_count_u_0228.out");

	while (std::getline(inputFile, line)) {
		std::istringstream lines(line);
		lines >> sp;
		lines >> cnt;
		uint32_t sp_ = ((uint32_t) stoi(sp));
		uint32_t cnt_ = ((uint32_t) stoi(cnt));
		species_tus[sp_] = cnt_;
		for (uint32_t i = 0; i < cnt_; i++)
			ht_u->map_sp[sp_].push_back(NULL);
	}
	inputFile.close();

	inputFile.open("unique_substring_count_d_0228.out");

	while (std::getline(inputFile, line)) {
		std::istringstream lines(line);
		lines >> sp;
		lines >> cnt;
		uint32_t sp_ = ((uint32_t) stoi(sp));
		uint32_t cnt_ = ((uint32_t) stoi(cnt));
		species_tds[sp_] = cnt_;
		for (uint32_t i = 0; i < cnt_; i++)
			ht_d->map_sp[sp_].push_back(NULL);
	}
	inputFile.close();*/

	//ht_u->loadIdx64_test(IDXFILEU);
	//ht_d->loadIdx64_test(IDXFILED);
	ht_u->loadIdx64_p(IDXFILEU);
	ht_d->loadIdx64_p(IDXFILED);
		
	fprintf(stderr, "Loaded index file.\n");

	//std::string uf_main = "index_all_u_main.bin1", uf_extra = "index_all_u_extra.bin1";
	//std::string df_main = "index_all_d_main.bin2", df_extra = "index_all_d_extra.bin2";
	//ht_u->encodeIdx64_temp(uf_main, uf_extra, 0);
	//ht_d->encodeIdx64_temp(df_main, df_extra, 1);
	//ht->testbucket(6421955);
}

void FqReader::loadSmap() {
	std::string line, fn, sp;
	std::ifstream inputFile;

	/* Read in fasta file. */ 
	inputFile.open(MAPFILE);
	//fprintf(stderr, "M: %s\n", MAPFILE.c_str());

	size_t i = 0;
	while (std::getline(inputFile, line)) {
		//fprintf(stderr, "%s\n", line.c_str());
		std::istringstream lines(line);
		lines >> fn;
		lines >> sp;
		uint32_t sp_ = ((uint32_t) stoi(sp));
		read_cnts_u[sp_] = 0;
		read_cnts_d[sp_] = 0;
		species_order[sp_] = i++;
		//ht_u->map_sp[sp_] = NULL;
		//ht_u->map_cnt[sp_] = std::make_pair(0, 0);
	}
	inputFile.close();

	fprintf(stderr, "Loaded species map file.\n");

	//for (auto rc : read_cnts)
	//	fprintf(stderr, "%u : %lu\n", rc.first, rc.second);
}

void FqReader::loadGenomeLength() {
	std::string line, sp, gl, nul;
	std::ifstream inputGLFile, inputULFile, inputDLFile;

	inputGLFile.open("../Human_Gut/db_test/genome_lengths_gut.out");
	while (std::getline(inputGLFile, line)) {
		std::istringstream lines(line);
		lines >> sp;
		uint32_t sp_ = ((uint32_t) stoi(sp));
		lines >> gl;
		genome_lengths[sp_] = ((uint32_t) stoi(gl));
	}
	inputGLFile.close();

	inputULFile.open("unique_lmer_count_u_humangut.out");
	while (std::getline(inputULFile, line)) {
		std::istringstream lines(line);
		lines >> sp;
		uint32_t sp_ = ((uint32_t) stoi(sp));
		lines >> nul;
		species_nus[sp_] = ((uint32_t) stoi(nul));
	}
	inputULFile.close();

	inputDLFile.open("unique_lmer_count_d_humangut.out");
	while (std::getline(inputDLFile, line)) {
		std::istringstream lines(line);
		lines >> sp;
		uint32_t sp_ = ((uint32_t) stoi(sp));
		lines >> nul;
		species_nds[sp_] = ((uint32_t) stoi(nul));
	}
	inputDLFile.close();

	fprintf(stderr, "Loaded genome length file.\n");
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
					//query32_m();
					//resetReads();
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

void FqReader::queryFastq_p(std::vector<std::string> &qfilenames_) {		
	for (auto fq_file : qfilenames_) {
		readFastq(fq_file);
		query64_p();
		//reformatIdx32_p();
		//resetReads();
		//runILP_a(100, 50, 30, 0.04);
		loadGenomeLength();
		runILP_p(125, 100, 10000, erate_, 100.0, 0.0001, 0.01);
	}

	//fprintf(stderr, "-Number of unlabeled reads: %lu.\n", nundet);
	//fprintf(stderr, "-Number of reads with conflict labels: %lu.\n", nconf);
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
					//query32_m();
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
 	std::string bases, sp;
	std::ifstream inputFile;
	size_t rl;

	/* Read in fasta file. */ 
	inputFile.open(INFILE);
	while (std::getline(inputFile, bases)) {
		//std::istringstream line(bases);
		//line >> sp;
		//true_labels.push_back((uint32_t) stoi(sp.substr(1)));
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
/*std::set<uint32_t> &queryUnique() {
}*/

void FqReader::query32() {
	fprintf(stderr, "Query 32. 0112\n");
	uint32_t trid = 556258;
	int fp_cnt = 0;
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
			/*if (nrd == 104) {
				for (int j = 0; j < 100; j++)
					fprintf(stderr, "%c", (char) read[j]);
				fprintf(stderr, "\n%lu, %u\n", i, hv);
				tmp_rid = ht_u->find32_test(hv, read + i + hash_len_u, rl - hash_len_u - i);
			}*/
			
			//fprintf(stderr, "%u; %lu; %lu\n", hv, hash_len_, i);
			tmp_rid = ht_u->find32(hv, read + i + hash_len_u, rl - hash_len_u - i);
			
			//GACAGCGGATCTGATGGTGCCGCCCCGTCTTCAAATCGATCTCCAACAGATAATAATTATCCGAATGGGCGATCAATTTATAACGCAGGATGGCTTTCTT
			if (tmp_rid != 0) {
				if (rid == 0)
					rid = tmp_rid;
				else if (tmp_rid != rid) {
					//fprintf(stderr, "%d\n", symbolIdx[read[76]]);
					//fprintf(stderr, "\n%lu, %u, %u\n", nrd, rid, tmp_rid);
					//abort();
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
			else {
				read_cnts_u[rid]++;
				if (rid != trid)
					fp_cnt++;
			}
		}
		if (nrd++ % 100000 == 0)
			fprintf(stderr, "Processed %lu reads.\n", nrd);
	}

	for (auto rc : read_cnts_u)
	{
		fprintf(stderr, "%u\t%lu\n", rc.first, rc.second);
	}

	fprintf(stderr, "Number of unlabeled reads: %lu.\n", nundet);
	fprintf(stderr, "Number of reads with conflict labels: %lu.\n", nconf);
	fprintf(stderr, "Number of reads tp labels: %lu.\n", read_cnts_u[trid]);
	fprintf(stderr, "Number of reads fp labels: %u.\n", fp_cnt);

	delete[] rc_read;
}

void FqReader::query64() {
}

/* Query L-mers covered by doubly-unique substrings. */
void FqReader::query32_d() {
	fprintf(stderr, "Query 32 du. 0112\n");
	uint32_t trid = 556258;
	int fp_cnt = 0;

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
			/*if (nrd == 3 && tmp_rid.first != 0) {
				for (int j = 0; j < 100; j++)
					fprintf(stderr, "%c", (char) rc_read[j]);
				fprintf(stderr, "-RC\n");
				fprintf(stderr, "%u, %u, %u", tmp_rid.first, temp_rid.second, hash_len_);
			}*/
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
						read_cnts_d[rid_pair.first]++;
					if (rid_pair.second != 0)
						read_cnts_d[rid_pair.second]++;
					if (rid_pair.first != 0 && rid_pair.first != trid && rid_pair.second != 0 && rid_pair.second != trid)
						fp_cnt += 1;
				}
				break;
			default:
				std::set<uint32_t> intersection;
				int i = 0;
				for (auto rid_pair : rid_pairs) {
					//if (nrd < 20)
					//fprintf(stderr, "(%u, %u)\t", rid_pair.first, rid_pair.second);
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
				//if (nrd < 20)
				//fprintf(stderr, "%u\t%lu\n", nrd, intersection.size());
				switch (intersection.size()) {
					case 0: 
						nundet++;
						break;
					case 1:
						for (auto rid : intersection) {
							read_cnts_d[rid]++;
							if (rid != trid)
								fp_cnt++;
						}
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

	for (auto rc : read_cnts_d)
	{
		fprintf(stderr, "%u\t%lu\n", rc.first, rc.second);
	}
	fprintf(stderr, "Number of unlabeled reads: %lu.\n", nundet);
	fprintf(stderr, "Number of reads with conflict labels: %lu.\n", nconf);
	fprintf(stderr, "Number of reads tp labels: %lu.\n", read_cnts_d[trid]);
	fprintf(stderr, "Number of reads fp labels: %u.\n", fp_cnt);

	//delete []rc_read;
	fprintf(stderr, "BB\n");
}

void FqReader::query64_d() {
}

void FqReader::query32_p() {
	uint32_t hv = 0, hs = 0;
	//uint32_t nerror = 0;
	pleafNode *pln;
	std::set<uint32_t> intersection;
	std::set<pleafNode*> pnodes;
	std::set<uint32_t> rids;
	std::set<std::pair<uint32_t, uint32_t>> rid_pairs;
	//std::unordered_set<pleafNode*> results;
	size_t rl = 0;
	uint8_t *rc_read = new uint8_t[max_rl];
	size_t nrd = 0;

	for (auto read : reads) {
		rids.clear();
		rid_pairs.clear();
		pnodes.clear();
		intersection.clear();
		rl = 100;

		/* Forward strand. */
		/* Query unique substrings. */
		hs = 2 * hash_len_u - 2;
		hv = 0;
		for (size_t i = 0; i < hash_len_u; i++)
			hv = ((hv << 2) | symbolIdx[read[i]]);
		for (size_t i = 0; i < rl - hash_len_u; i++) {
			pln = ht_u->find32_p(hv, read + i + hash_len_u, rl - hash_len_u - i);
			if (pln != NULL)
				pnodes.insert(pln);
			hv = hv - (symbolIdx[read[i]] << hs); /* Next hash. */
			hv = ((hv << 2) | symbolIdx[read[i + hash_len_u]]);
		}
		pln = ht_u->find32_p(hv, read, 0);
		if (pln != NULL)
			pnodes.insert(pln);

		/* Query doubly unique substrings. */
		/*hs = 2 * hash_len_d - 2;
		hv = 0;
		for (size_t i = 0; i < hash_len_d; i++)
			hv = ((hv << 2) | symbolIdx[read[i]]);
		for (size_t i = 0; i < rl - hash_len_d; i++) {
			pln = ht_d->find32_p(hv, read + i + hash_len_d, rl - hash_len_d - i);
			if (pln != NULL)
				pnodes.insert(pln);
			hv = hv - (symbolIdx[read[i]] << hs); */
			/* Next hash. */
			/*hv = ((hv << 2) | symbolIdx[read[i + hash_len_d]]);
		}
		pln = ht_d->find32_p(hv, read, 0);
		if (pln != NULL)
			pnodes.insert(pln);*/

		/* Reverse complement. */
		getRC(rc_read, read, rl);

		/* Query unique substrings. */
		hs = 2 * hash_len_u - 2;
		hv = 0;
		for (size_t i = 0; i < hash_len_u; i++)
			hv = ((hv << 2) | symbolIdx[rc_read[i]]);

		for (size_t i = 0; i < rl - hash_len_u; i++) {
			pln = ht_u->find32_p(hv, rc_read + i + hash_len_u, rl - hash_len_u - i);
			if (pln != NULL)
				pnodes.insert(pln);
			hv = hv - (symbolIdx[rc_read[i]] << hs); /* Next hash. */
			hv = ((hv << 2) | symbolIdx[rc_read[i + hash_len_u]]);
		}
		pln = ht_u->find32_p(hv, rc_read, 0);
		if (pln != NULL)
			pnodes.insert(pln);
				
		/* Query doubly unique substrings. */
		/*hs = 2 * hash_len_d - 2;
		hv = 0;
		for (size_t i = 0; i < hash_len_d; i++)
			hv = ((hv << 2) | symbolIdx[rc_read[i]]);
		for (size_t i = 0; i < rl - hash_len_d; i++) {
			pln = ht_d->find32_p(hv, rc_read + i + hash_len_d, rl - hash_len_d - i);
			if (pln != NULL)
				pnodes.insert(pln);
			hv = hv - (symbolIdx[rc_read[i]] << hs); */
			/* Next hash. */
			/*hv = ((hv << 2) | symbolIdx[rc_read[i + hash_len_d]]);
		}
		pln = ht_d->find32_p(hv, rc_read, 0);
		if (pln != NULL)
			pnodes.insert(pln);*/
		
		/* Record results. */
		int i = 0;
		for (auto pn : pnodes)
			if (pn->refID2 == 0)
				rids.insert(pn->refID1);
			else
				if (pn->refID1 < pn->refID2)
					rid_pairs.insert(std::make_pair(pn->refID1, pn->refID2));
				else
					rid_pairs.insert(std::make_pair(pn->refID2, pn->refID1));
		/*bool fff = 0;
		for (auto pn : pnodes)
			if (pn->refID.size() == 1) {
				if (pn->refID[0] == 1226322 || pn->refID[0] == 1226323)
					fff = 1;
			} else {
				if (pn->refID[0] == 1226322 || pn->refID[1] == 1226322 || pn->refID[0] == 1226323 || pn->refID[1] == 1226323)
					fff = 1;
			}
		if (fff) {
			for (auto pn : pnodes) {
				if (pn->refID.size() == 1)
					fprintf(stderr, "rid = %u\t", pn->refID[0]);
				else
					fprintf(stderr, "rid1 = %u; rid2 = %u\t", pn->refID[0], pn->refID[1]);
			}
			fprintf(stderr, "\n");
		}*/

		switch (rid_pairs.size()) {
			case 0:
				if (rids.empty())
					nundet++;
				else {
					if (rids.size() == 1) {
						for (auto rid : rids)
							read_cnts_u[rid]++;
						for (auto pn : pnodes)
							pn->rcount += 1.0;
					} else
						nconf++;
				}
				break;
			case 1:
				if (rids.empty()) {
					for (auto rid_pair : rid_pairs) {
						read_cnts_d[rid_pair.first]++;
						read_cnts_d[rid_pair.second]++;
						//if (rid_pair.first != 658661 && rid_pair.second != 658661)
						//	nerror++;
					}
					for (auto pn : pnodes)
						pn->rcount += 1.0;
				} else {
					if (rids.size() > 1)
						nconf++;
					else {
						for (auto rid : rids)
							for (auto rid_pair : rid_pairs)
								if (rid_pair.first != rid && rid_pair.second != rid) {
									nconf++;
								} else {
									read_cnts_u[rid]++;
									read_cnts_d[rid]++;
									for (auto pn : pnodes)
										pn->rcount += 1.0;
								}
					}
				}
				break;
			default:
				if (!rids.empty()) {
					if (rids.size() > 1) {
						nconf++;
						break;
					}
					for (auto rid : rids) {
						bool conf = 0;
						for (auto rid_pair : rid_pairs)
							if (rid_pair.first != rid && rid_pair.second != rid) {
								conf = 1;
								break;
							}
						if (conf)
							nconf++;
						else {
							read_cnts_u[rid]++;
							read_cnts_d[rid]++;
							for (auto pn : pnodes)
								pn->rcount += 1.0;
						}
					}
				} else {
					i = 0;
					for (auto rid_pair : rid_pairs) {
						if (i == 0) {	
							intersection.insert(rid_pair.first);
							intersection.insert(rid_pair.second);
						} else {
							std::vector<uint32_t> to_be_removed;
							for (auto rid : intersection) {
								if (rid_pair.first != rid && rid_pair.second != rid)
									to_be_removed.push_back(rid);
							}
							for (auto rid : to_be_removed)
								intersection.erase(rid);
						}
						i++;
					}
					switch (intersection.size()) {
						case 0: 
							nconf++;
							break;
						case 1:
							for (auto rid : intersection)
								read_cnts_d[rid]++;
							for (auto pn : pnodes)
								pn->rcount += 1.0;
							break;
						default:
							nconf++;
							break;
					}
				}
				break;
		}
		if (nrd++ % 100000 == 0) {
			fprintf(stderr, "Processed %lu reads.\n", nrd);
		}
	}

	fprintf(stderr, "Number of unlabeled reads: %lu.\n", nundet);
	fprintf(stderr, "Number of reads with conflict labels: %lu.\n", nconf);
	//fprintf(stderr, "Number of reads with incorrect labels: %lu.\n", nerror);
	//delete []rc_read;
	for (auto rc : read_cnts_u)
		fprintf(stderr, "%u\t%lu\t%lu\n", rc.first, rc.second, read_cnts_d[rc.first]);

	//fprintf(stderr, "Result Size: %lu", results.size());
	/*for (auto pn : results) {
		us_lengths.push_back(pn->depth);
		us_refids.push_back(pn->refID);
		us_ucounts.push_back(pn->ucount);
		us_rcounts.push_back(pn->rcount);
	}*/
	fprintf(stderr, "Completed query.\n");
}

void FqReader::query64_p() {
	auto start = std::chrono::high_resolution_clock::now();

	uint64_t hv = 0;
	uint32_t hs = 0;
	pleafNode *pln;
	std::set<uint32_t> intersection;
	std::set<pleafNode*> pnodes;
	std::set<uint32_t> rids;
	std::set<std::pair<uint32_t, uint32_t>> rid_pairs;
	
	size_t rl = 0;
	uint8_t *rc_read = new uint8_t[max_rl];
	size_t nrd = 0;
	/*For measuring accuracy*/ size_t misclass = 0;

	for (auto read : reads) {
		rids.clear();
		rid_pairs.clear();
		pnodes.clear();
		intersection.clear();
		rl = 125;	
		
		//fprintf(stderr, "Query forward strand.\n");
		/* Forward strand. */
		/* Query unique substrings. */
		hs = 2 * hash_len_u - 2;
		hv = 0;
		for (size_t i = 0; i < hash_len_u; i++)
			hv = ((hv << 2) | symbolIdx[read[i]]);
		for (size_t i = 0; i < rl - hash_len_u; i++) {
			pln = ht_u->find64_p(hv, read + i + hash_len_u, rl - hash_len_u - i);
			if (pln != NULL)
				pnodes.insert(pln);
			pln = ht_d->find64_p(hv, read + i + hash_len_u, rl - hash_len_u - i);
			if (pln != NULL)
				pnodes.insert(pln);
			//if (nrd == 0)
			//	fprintf(stderr, "%lu\n", hv);
			hv = hv - ((uint64_t) symbolIdx[read[i]] << hs); /* Next hash. */
			hv = ((hv << 2) | symbolIdx[read[i + hash_len_u]]);
		}
		pln = ht_u->find64_p(hv, read, 0);
		if (pln != NULL)
			pnodes.insert(pln);
		pln = ht_d->find64_p(hv, read, 0);
		if (pln != NULL)
			pnodes.insert(pln);

		//fprintf(stderr, "Query rc strand.\n");
		/* Reverse complement. */
		getRC(rc_read, read, rl);

		/* Query unique substrings. */
		hs = 2 * hash_len_u - 2;
		hv = 0;
		for (size_t i = 0; i < hash_len_u; i++)
			hv = ((hv << 2) | symbolIdx[rc_read[i]]);

		for (size_t i = 0; i < rl - hash_len_u; i++) {
			pln = ht_u->find64_p(hv, rc_read + i + hash_len_u, rl - hash_len_u - i);
			if (pln != NULL) {
				//fprintf(stderr, "-ID: %u\n", pln->refID1);
				pnodes.insert(pln);
			}
			pln = ht_d->find64_p(hv, rc_read + i + hash_len_u, rl - hash_len_u - i);
			if (pln != NULL) {
				//fprintf(stderr, "-ID: %u\n", pln->refID1);
				pnodes.insert(pln);
			}
			hv = hv - ((uint64_t) symbolIdx[rc_read[i]] << hs); /* Next hash. */
			hv = ((hv << 2) | symbolIdx[rc_read[i + hash_len_u]]);
		}
		pln = ht_u->find64_p(hv, rc_read, 0);
		if (pln != NULL) {
			//fprintf(stderr, "-ID: %u\n", pln->refID1);
			pnodes.insert(pln);
		}
		pln = ht_d->find64_p(hv, rc_read, 0);
		if (pln != NULL) {
			//fprintf(stderr, "-ID: %u\n", pln->refID1);
			pnodes.insert(pln);
		}
				
		//fprintf(stderr, "%lu\n", pnodes.size());
		
		/* Record results. */
		int i = 0;
		for (auto pn : pnodes) {
			//fprintf(stderr, "HERE; %lu\n", nrd);
			//fprintf(stderr, "ID: %u\n", pn->refID1);
			if (pn->refID2 == 0)
				rids.insert(pn->refID1);
			else {
				//fprintf()
				if (pn->refID1 < pn->refID2)
					rid_pairs.insert(std::make_pair(pn->refID1, pn->refID2));
				else
					rid_pairs.insert(std::make_pair(pn->refID2, pn->refID1));
			}
		}
		
		switch (rid_pairs.size()) {
			case 0:
				if (rids.empty())
					nundet++;
				else {
					if (rids.size() == 1) {
						for (auto rid : rids) {
							read_cnts_u[rid]++;
							//if (rid != true_labels[nrd])
							//	misclass++;
						}
						for (auto pn : pnodes)
							pn->rcount += 1;
					} else
						nconf++;
				}
				break;
			case 1:
				if (rids.empty()) {
					for (auto rid_pair : rid_pairs) {
						read_cnts_d[rid_pair.first]++;
						read_cnts_d[rid_pair.second]++;
						//if (rid_pair.first != 658661 && rid_pair.second != 658661)
						//	nerror++;
						//if (rid_pair.first != true_labels[nrd] && rid_pair.second != true_labels[nrd])
						//	misclass++;
					}
					for (auto pn : pnodes)
						pn->rcount += 1;
				} else {
					if (rids.size() > 1)
						nconf++;
					else {
						for (auto rid : rids)
							for (auto rid_pair : rid_pairs)
								if (rid_pair.first != rid && rid_pair.second != rid) {
									nconf++;
								} else {
									read_cnts_u[rid]++;
									read_cnts_d[rid]++;
									for (auto pn : pnodes)
										pn->rcount += 1;
									//if (rid != true_labels[nrd])
									//	misclass++;
								}
					}
				}
				break;
			default:
				if (!rids.empty()) {
					if (rids.size() > 1) {
						nconf++;
						break;
					}
					for (auto rid : rids) {
						bool conf = 0;
						for (auto rid_pair : rid_pairs)
							if (rid_pair.first != rid && rid_pair.second != rid) {
								conf = 1;
								break;
							}
						if (conf)
							nconf++;
						else {
							read_cnts_u[rid]++;
							read_cnts_d[rid]++;
							for (auto pn : pnodes)
								pn->rcount += 1;
							//if (rid != true_labels[nrd])
							//	misclass++;
						}
					}
				} else {
					i = 0;
					for (auto rid_pair : rid_pairs) {
						if (i == 0) {	
							intersection.insert(rid_pair.first);
							intersection.insert(rid_pair.second);
						} else {
							std::vector<uint32_t> to_be_removed;
							for (auto rid : intersection) {
								if (rid_pair.first != rid && rid_pair.second != rid)
									to_be_removed.push_back(rid);
							}
							for (auto rid : to_be_removed)
								intersection.erase(rid);
						}
						i++;
					}
					switch (intersection.size()) {
						case 0: 
							nconf++;
							break;
						case 1:
							for (auto rid : intersection) {
								read_cnts_d[rid]++;
								//if (rid != true_labels[nrd])
								//	misclass++;
							}
							for (auto pn : pnodes)
								pn->rcount += 1;
							break;
						default:
							nconf++;
							break;
					}
				}
				break;
		}
		if (nrd++ % 100000 == 0) {
			fprintf(stderr, "Processed %lu reads.\n", nrd);
		}
	}

	fprintf(stderr, "Number of unlabeled reads: %lu.\n", nundet);
	fprintf(stderr, "Number of reads with conflict labels: %lu.\n", nconf);
	fprintf(stderr, "Number of misclassified reads: %lu.\n", misclass);
	//fprintf(stderr, "Number of reads with incorrect labels: %lu.\n", nerror);
	//delete []rc_read;
	//for (auto rc : read_cnts_u)
	//	fprintf(stderr, "%u\t%lu\t%lu\n", rc.first, rc.second, read_cnts_d[rc.first]);

	fprintf(stderr, "Completed query.\n");
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - start).count();
	fprintf(stderr, "Time for query: %lu ms.\n", duration);
}

void FqReader::runILP_p(int rl, int read_cnt_thres, uint32_t unique_thres, 
			double erate, double max_cov, double resolution, double epsilon) {
	auto start = std::chrono::high_resolution_clock::now();
	
	// Initialize cplex environment. 
	IloEnv env;
	IloModel model(env);
	IloExpr objective(env);
	
	// Binary variable EXIST[l]: Existence of species l.
	size_t n_species = read_cnts_u.size();
	exist = new int[n_species];
	memset(exist, 1, sizeof(int) * n_species);
	IloBoolVarArray EXIST(env, n_species);
	fprintf(stderr, "ILP environment in CPLEX initialized.\n");

	// Constraint 0: If the number of reads being assigned to a particular 
	//	species l is less than READ_CNT_THRES, then EXIST[l] is set to 0
	for (auto rc : read_cnts_u) {
		uint32_t rid = rc.first;
		double d1 = rc.second * 1.0, d2 = rc.second * 1.0;
		d1 -= read_cnt_thres;
		d2 -= (1.0 * species_nus[rid]) * resolution;
		//fprintf(stderr, "%u; %.4f; %.4f; %u; %u\n", rid, d1, d2, species_nus[rid], unique_thres);
		if (species_nus[rid] >= unique_thres) {
			if (d1 < 0.0 || d2 < 0.0) {
				exist[species_order[rid]] = 0;
				model.add(EXIST[species_order[rid]] == 0);
			} 
		} else {
			if (d2 < 0.0) {
				exist[species_order[rid]] = 0;
				model.add(EXIST[species_order[rid]] == 0);
			}
		}
	}
	for (auto rc : read_cnts_d) {
		uint32_t rid = rc.first;
		float d1 = rc.second * 1.0, d2 = rc.second * 1.0;
		d1 -= read_cnt_thres;
		d2 -= (1.0 * species_nds[rid]) * resolution;
		//fprintf(stderr, "%u; %.4f; %.4f\n", rid, d1, d2);
		if (species_nus[rid] >= unique_thres) {
			if (d1 < 0.0 || d2 < 0.0) {
				exist[species_order[rid]] = 0;
				model.add(EXIST[species_order[rid]] == 0);
			} 
		} else {
			if (d2 < 0.0) {
				exist[species_order[rid]] = 0;
				model.add(EXIST[species_order[rid]] == 0);
			} 
		}
	}
	size_t n_species_exist = 0;
	for (auto rc : read_cnts_u) {
		uint32_t rid = rc.first;
		if (exist[species_order[rid]]) {
			n_species_exist++;
			fprintf(stderr, "Species %u exists. \n", rid);
			//if (species_nus[rid] >= unique_thres)
			//	model.add(EXIST[species_order[rid]] == 1);
		}
	}
	fprintf(stderr, "\tConstructed constraint 0.\n");
	for (auto rc : read_cnts_u)
		if (exist[species_order[rc.first]])
			fprintf(stderr, "%u\t%lu\t%lu\n", rc.first, rc.second, read_cnts_d[rc.first]);


	size_t num_us_valid = 0, abs_index = 0;
	for (auto so : species_order) {
		uint32_t rid = so.first;
		if (exist[so.second] > 0) {
			//num_us_valid += ht_u->map_cnt[rid].first;
			//num_us_valid += ht_u->map_cnt[rid].second;
			num_us_valid += ht_u->map_sp[rid].size();
			num_us_valid += ht_d->map_sp[rid].size();
		}
	}
	fprintf(stderr, "Filtered out invalid species by read counts.\n");

	// Continuous variable C[l]: Coverage of species l.
	IloNumVarArray COV(env, n_species, 0.0, max_cov, ILOFLOAT);
	//IloNumVarArray ABS(env, num_us_valid, 0, 99999, ILOFLOAT);
	fprintf(stderr, "\tConstructed all other variables.\n");

	//for (size_t i = 0; i < num_us_valid; i++)
	//	objective += ABS[i];
	//model.add(IloMinimize(env, objective));
	fprintf(stderr, "\tConstructed the objective function.\n");
	
	// Minimize the total number of species. 
	for (auto so : species_order) {
		uint32_t rid = so.first;
		double u_factor = 1000.0 / ht_u->map_sp[rid].size();
		if (exist[so.second] > 0) {
			for (auto pn : ht_u->map_sp[rid]) {
				assert((rid == pn->refID1) || (rid == pn->refID2));
				size_t si = species_order[rid];
				double wcov = (pn->ucount1 * (rl - pn->depth) * 1.0 / rl);
				wcov = wcov * pow(1 - erate, pn->depth);
				//model.add(ABS[abs_index] >= (wcov * COV[si] - pn->rcount) * u_factor);
				//model.add(ABS[abs_index] >= (pn->rcount - wcov * COV[si]) * u_factor);
				//if (rid == 46354 || rid == 2009331)
				//	fprintf(stderr, "%u\t%u\n", rid, pn->rcount);	
				objective += (wcov * COV[si] - pn->rcount) * (wcov * COV[si] - pn->rcount) * u_factor;
				//objective += IloAbs(wcov * COV[si] - pn->rcount) * u_factor;
				abs_index++;
			} 
		}
	}
	
	/* Statistical test. */
	/*for (auto so : species_order) {
		uint32_t rid = so.first;
		std::vector<uint32_t> rcnts;
		std::vector<uint8_t> ulens;
		if (exist[so.second] > 0 && read_cnts_d[rid] >= read_cnt_thres) {
			uint32_t sum_rcnt = 0;
			uint32_t sum_ulen = 0;
			for (size_t pi = 0; pi < species_tds[rid]; pi++) {
				rcnts.push_back(ht_d->map_sp[rid][pi]->rcount);
				ulens.push_back(rl - ht_d->map_sp[rid][pi]->depth);
				sum_rcnt += ht_d->map_sp[rid][pi]->rcount;
				sum_ulen += (rl - ht_d->map_sp[rid][pi]->depth);
			}*/
			/*std::vector<uint32_t> dists;
			uint32_t pi = 0, last_pi = 0;
			for (auto cnt : rcnts) {
				if (cnt > 0) {
					dists.push_back(pi - last_pi);
					for (uint32_t j = 0; j <= cnt - 1; j++)
						dists.push_back(0);
					last_pi = pi;
				}
				pi++;
			}
			double t = 0.0, mean = 0.0, sigma = 0.0;
			for (auto dis : dists)
				t += 1.0 * dis;
			mean = t / dists.size();
			for (auto dis : dists)
				sigma += (dis - mean) * (dis - mean);
			sigma = sigma / dists.size();
			sigma = std::sqrt(sigma);
			t = (mean - species_tus[rid] * 1.0 / sum_rcnt * 1.0);
			t *= (std::sqrt(dists.size()) / sigma);
			fprintf(stderr, "%u\t%u\t%u\t%.4f\t%.4f\t%.4f\t%.4f\n", rid, species_tus[rid], sum_rcnt, mean, sigma, species_tus[rid] * 1.0 / sum_rcnt, t);*/
			/*double chi_sq_val = 0.0;
			double exp = sum_rcnt * 1.0 / species_tus[rid];
			for (auto cnt : rcnts) {
				if (fabs(exp - cnt) > 1.0)
					chi_sq_val += (exp - cnt) * (exp - cnt) / exp;
				else
					chi_sq_val += (exp - cnt) * (exp - cnt) * (exp - cnt) * (exp - cnt) / exp;
			}
			fprintf(stderr, "%u\t%u\t%u\t%.4f\t%.4f\n", rid, species_tus[rid], sum_rcnt, exp, chi_sq_val);*/
			//double ks_stat = 0.0;
			/*double ocnt = 0.0, ecnt = 0.0, max_gap = 0.0;
			for (size_t ci = 0; ci < rcnts.size(); ci++) {
				ocnt += rcnts[ci];
				double curr = ocnt / sum_rcnt;
				ecnt += ulens[ci];
				double ecurr = ecnt /  sum_ulen;
				if (fabs(curr - ecurr) > max_gap)
					max_gap = fabs(curr - ecurr);
			}
			fprintf(stderr, "%u\t%u\t%u\t%.4f\n", rid, species_tds[rid], sum_rcnt, max_gap);
		}
	}*/

	for (auto so : species_order) {
		uint32_t rid = so.first;
		double d_factor = 1000.0 / ht_d->map_sp[rid].size();
		if (exist[so.second] > 0) {
			for (auto pn : ht_d->map_sp[rid]) {
				assert((rid == pn->refID1) || (rid == pn->refID2));
				size_t si1 = species_order[pn->refID1], si2 = species_order[pn->refID2];
				double wcov1 = (pn->ucount1 * (rl - pn->depth) * 1.0 / rl);
				double wcov2 = (pn->ucount2 * (rl - pn->depth) * 1.0 / rl);
				wcov1 = wcov1 * pow(1 - erate, pn->depth);
				wcov2 = wcov2 * pow(1 - erate, pn->depth);
				//model.add(ABS[abs_index] >= (wcov1 * COV[si1] + wcov2 * COV[si2] - pn->rcount) * d_factor);
				//model.add(ABS[abs_index] >= (pn->rcount - wcov1 * COV[si1] - wcov2 * COV[si2]) * d_factor);
				objective += (wcov1 * COV[si1] + wcov2 * COV[si2] - pn->rcount) * (wcov1 * COV[si1] + wcov2 * COV[si2] - pn->rcount) * d_factor;
				//objective += IloAbs(wcov1 * COV[si1] + wcov2 * COV[si2] - pn->rcount) * d_factor;	
				abs_index++;
			} 
		}
	}
	model.add(IloMinimize(env, objective));
	fprintf(stderr, "\tConstructed constraint 1-2.\n");

	// Constraint 3: COV should match EXIST
	for (size_t i = 0; i < n_species; i++) {
		model.add(COV[i] <= max_cov * EXIST[i]);
		model.add(COV[i] >= 0.01 * EXIST[i]);
		//model.add(COV[i] == 10.0 * EXIST[i]);
		//model.add(COV[i] == 10.0 * EXIST[i]);
	}
	/*for (auto rc : read_cnts_u) {
		uint32_t rid = rc.first;
		double min_cov = ((2.0 - epsilon) * rc.second * rl) / species_nus[rid];
		if (min_cov > 0.0)
			fprintf(stderr, "U-%u; %u; %u; %.4f\n", rid, rc.second, species_nus[rid], min_cov);
		model.add(COV[species_order[rid]] >= min_cov * EXIST[species_order[rid]]);
	}
	for (auto rc : read_cnts_d) {
		uint32_t rid = rc.first;
		double min_cov = ((2.0 - epsilon) * rc.second * rl) / species_nds[rid];
		//if (min_cov > 0.0)
		//	fprintf(stderr, "D-%u; %u; %u; %.4f\n", rid, rc.second, species_nus[rid], min_cov);
		if (species_nus[rid] >= unique_thres) 
			model.add(COV[species_order[rid]] >= min_cov * EXIST[species_order[rid]]);
	}
	*/
	IloExprArray EXP1(env, n_species_exist);
	size_t exp_i = 0;
	for (; exp_i < n_species_exist; exp_i++)
		EXP1[exp_i] = IloExpr(env);
	//fprintf(stderr, "n_species_exist: %lu\n", n_species_exist);
	exp_i = 0;
	for (auto so : species_order) {
		uint32_t rid = so.first;
		if (exist[so.second] > 0) {
			for (auto pn : ht_u->map_sp[rid]) {
				size_t si = species_order[rid];
				double wcov = (pn->ucount1 * (rl - pn->depth) * 1.0 / rl);
				wcov = wcov * pow(1 - erate, pn->depth);
				EXP1[exp_i] += wcov * COV[si];
			} 
			//fprintf(stderr, "exp_i: %lu\n", exp_i);
			if (species_nus[rid] >= unique_thres)
				model.add(EXP1[exp_i] * (1.0 + epsilon) - read_cnts_u[rid] >= 0);
			else
				model.add(EXP1[exp_i] >= 0.0);
			exp_i++;	
		}
	}
	IloExprArray EXP2(env, n_species_exist);
	for (exp_i = 0; exp_i < n_species_exist; exp_i++)
		EXP2[exp_i] = IloExpr(env);

	exp_i = 0;
	for (auto so : species_order) {
		uint32_t rid = so.first;
		if (exist[so.second] > 0) {
			for (auto pn : ht_d->map_sp[rid]) {
				size_t si1 = species_order[pn->refID1], si2 = species_order[pn->refID2];
				double wcov1 = (pn->ucount1 * (rl - pn->depth) * 1.0 / rl);
				double wcov2 = (pn->ucount2 * (rl - pn->depth) * 1.0 / rl);
				wcov1 = wcov1 * pow(1 - erate, pn->depth);
				wcov2 = wcov2 * pow(1 - erate, pn->depth);
				EXP2[exp_i] += wcov1 * COV[si1] + wcov2 * COV[si2];
			}
			//fprintf(stderr, "exp_i: %lu\n", exp_i);
			if (species_nus[rid] >= unique_thres)
				model.add(EXP2[exp_i] * (1.0 + epsilon) - read_cnts_d[rid] >= 0);
			else
				model.add(EXP2[exp_i] >= 0.0);
			exp_i++;
		}
	}
	fprintf(stderr, "\tConstructed constraint 3.\n");

	// Constraint 4: refining the search space
	IloExpr TOTAL(env);
	for (auto gl_map : genome_lengths)
		TOTAL += (COV[species_order[gl_map.first]] * (1.0 * gl_map.second) / rl);
	//model.add(TOTAL >= 0.9 * reads.size());
	model.add(TOTAL <= (1.0 + epsilon) * reads.size());
	fprintf(stderr, "\tConstructed constraint 4.\n");

	IloExpr e1(env);
	for (size_t i = 0; i < n_species; i++) {
		e1 += EXIST[i];
	}
	//model.add(e1 <= 2);
	//fprintf(stderr, "\tConstructed constraint 6.\n");
	int k = n_species_exist;
	FILE * fout = fopen(OUTPUTFILE.c_str(), "w");

	// Solving the ILP.
	try {
		//while (k >= 20) {
			model.add(e1 <= k * 1.0);
			IloCplex cplex(model);
			cplex.setParam(IloCplex::IntParam::Threads, 32);
			cplex.setParam(IloCplex::NumParam::TiLim, 10800);
			cplex.setParam(IloCplex::NumParam::SolnPoolAGap, 0.0);
			
			//cplex.solve();
		
			if (cplex.populate()) {
				int nsolns = cplex.getSolnPoolNsolns();
				for (int s = 0; s < nsolns; s++) {
					fprintf(fout, "\nk = %d; Objective = %.4f.\n", k, cplex.getObjValue(s));
					fprintf(fout, "\nSolution %d:\n", s);
					for (auto sp : species_order) {
						bool ei = cplex.getValue(EXIST[sp.second], s);
						float ci = cplex.getValue(COV[sp.second], s);
						if (ei != 0)
							fprintf(fout, "%u\t%d\t%.4f\n", sp.first, ei, ci);
					}
				}
			}
			//k--;
		//}
		fclose(fout);
	} catch (IloException &ex) {}

	if (exist != NULL)
		delete [] exist;

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - start).count();
	fprintf(stderr, "Time for ilp: %lu ms.\n", duration);
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

