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

Genome::Genome(uint32_t taxid, std::string &gn) {
	read_cnts_u = 0;
	read_cnts_d = 0;
	glength = 0;
	nus = 0;
	nds = 0;
	taxID = taxid;
	name = gn;
}

FqReader::FqReader(int d, int ho, uint32_t hl, std::string &idx_fn, std::string &map_fn) {
	doubly_unique = d;
	hash_option = ho;
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
	if (ht_u != NULL)
		delete ht_u;
	if (ht_d != NULL)
		delete ht_d;
}

void FqReader::loadIdx_p() {
	fprintf(stderr, "%s\n", IDXFILEU.c_str());
	
	ht_u->loadIdx64_p(IDXFILEU);
	ht_d->loadIdx64_p(IDXFILED);
		
	fprintf(stderr, "Loaded index file.\n");
}

void FqReader::loadSmap() {
	genomes.push_back(NULL);
	std::string line, gname, id, taxid;
	std::ifstream inputFile;

	/* Read in fasta file. */ 
	inputFile.open(MAPFILE);
	//fprintf(stderr, "%s\n", MAPFILE.c_str());

	while (std::getline(inputFile, line)) {
		std::istringstream lines(line);
		lines >> gname;
		lines >> id;
		lines >> taxid;
		lines >> gname;
		Genome *new_genome = new Genome((uint32_t) stoi(taxid), gname);
		genomes.push_back(new_genome);
	}
	inputFile.close();
	//for (size_t i = 0; i < 4122; i++)
	//	fprintf(stderr, "%lu, %u\t%lu\t%lu\n", i, genomes[i + 1]->taxID, genomes[i + 1]->read_cnts_u, genomes[i + 1]->read_cnts_d);
	fprintf(stderr, "Loaded genome map file.\n");
}

void FqReader::loadGenomeLength() {
	std::string line, id, gl, nu;
	std::ifstream inputGLFile, inputULFile, inputDLFile;

	inputGLFile.open("genome_lengths.out");
	while (std::getline(inputGLFile, line)) {
		std::istringstream lines(line);
		lines >> id;
		lines >> gl;
		//genomes[(uint32_t) stoi(id)]->glength = ((uint32_t) stoi(gl));
		genomes[stoi(id)]->glength = ((uint32_t) stoi(gl));
	}
	inputGLFile.close();

	inputULFile.open("unique_lmer_count_u.out");
	while (std::getline(inputULFile, line)) {
		std::istringstream lines(line);
		lines >> id;
		lines >> nu;
		//genomes[(uint32_t) stoi(id)]->nus = ((uint32_t) stoi(nu));
		genomes[stoi(id)]->nus = ((uint32_t) stoi(nu));
	}
	inputULFile.close();

	inputDLFile.open("unique_lmer_count_d.out");
	while (std::getline(inputDLFile, line)) {
		std::istringstream lines(line);
		lines >> id;
		lines >> nu;
		//genomes[(uint32_t) stoi(id)]->nds = ((uint32_t) stoi(nu));
		genomes[stoi(id)]->nds = ((uint32_t) stoi(nu));
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

void FqReader::queryFastq_p(std::vector<std::string> &qfilenames_) {		
	for (auto fq_file : qfilenames_) {
		readFastq(fq_file);
		//for (size_t i = 0; i < 4122; i++)
		//	fprintf(stderr, "%lu, %u\t%lu\t%lu\n", i, genomes[i + 1]->taxID, genomes[i + 1]->read_cnts_u, genomes[i + 1]->read_cnts_d);
		query64_p(100);
		loadGenomeLength();
		runILP_p(100, 100, 10000, erate_, 100.0, 0.0001, 0.01);
	}
}

void FqReader::readFastq(std::string &INFILE) { 
 	std::string bases, sp;
	std::ifstream inputFile;
	size_t rl;

	/* Read in fasta file. */ 
	inputFile.open(INFILE);
	while (std::getline(inputFile, bases)) {
		std::getline(inputFile, bases);
		rl = bases.length();
		uint8_t *read = new uint8_t[rl];
		strncpy((char*)read, bases.c_str(), rl);
		reads.push_back(read);
		std::getline(inputFile, bases);
		std::getline(inputFile, bases);
	}
	
	inputFile.close();
	fprintf(stderr, "%lu\n", reads.size());
}

void FqReader::getRC(uint8_t *dest, uint8_t *srcs, size_t length) {
	//fprintf(stderr, "LENGTH: %lu\n", length);
	for(size_t i = 0; i < length; i++) {
		//fprintf(stderr, "%lu, %d, %u\n", i, rcIdx[srcs[length - i - 1]], srcs[length - i - 1]); 
		dest[i] = rcIdx[srcs[length - i - 1]];
		/*if (genomes[3]->read_cnts_u > 0) {
			fprintf(stderr, "Begin with error 3.\n");
			abort();
		}*/
	}
	/*if (genomes[3]->read_cnts_u > 0) {
		fprintf(stderr, "Begin with error 3.\n");
		abort();
	}*/
}

void FqReader::query64_p(size_t rl) {
	auto start = std::chrono::high_resolution_clock::now();
	uint64_t hv = 0;
	uint32_t hs = 0;
	pleafNode *pln;
	std::set<uint32_t> intersection;
	std::set<pleafNode*> pnodes;
	std::set<uint32_t> rids;
	std::set<std::pair<uint32_t, uint32_t>> rid_pairs;
	uint8_t *rc_read = new uint8_t[max_rl];
	size_t nrd = 0;

	//if (genomes[3]->read_cnts_u > 0)
	//	fprintf(stderr, "Begin with error.\n"); 

	for (auto read : reads) {
		rids.clear();
		rid_pairs.clear();
		pnodes.clear();
		intersection.clear();

		/*if (genomes[3]->read_cnts_u > 0) {
			fprintf(stderr, "Begin with error 1.\n");
			abort();
		}*/

		/* Forward strand. */
		/* Query unique and doubly-unique substrings. */
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
			hv = hv - ((uint64_t) symbolIdx[read[i]] << hs); /* Next hash. */
			hv = ((hv << 2) | symbolIdx[read[i + hash_len_u]]);
		}
		pln = ht_u->find64_p(hv, read, 0);
		if (pln != NULL)
			pnodes.insert(pln);
		pln = ht_d->find64_p(hv, read, 0);
		if (pln != NULL)
			pnodes.insert(pln);

		/*if (genomes[3]->read_cnts_u > 0) {
			fprintf(stderr, "Begin with error 2.\n");
			abort();
		}*/

		/* Reverse complement. */
		getRC(rc_read, read, rl);

		/*if (genomes[3]->read_cnts_u > 0) {
			fprintf(stderr, "Begin with error 3.\n");
			abort();
		}*/

		/* Query unique and doubly-unique substrings. */
		hs = 2 * hash_len_u - 2;
		hv = 0;
		for (size_t i = 0; i < hash_len_u; i++)
			hv = ((hv << 2) | symbolIdx[rc_read[i]]);

		for (size_t i = 0; i < rl - hash_len_u; i++) {
			pln = ht_u->find64_p(hv, rc_read + i + hash_len_u, rl - hash_len_u - i);
			if (pln != NULL)
				pnodes.insert(pln);
			pln = ht_d->find64_p(hv, rc_read + i + hash_len_u, rl - hash_len_u - i);
			if (pln != NULL)
				pnodes.insert(pln);
			hv = hv - ((uint64_t) symbolIdx[rc_read[i]] << hs); /* Next hash. */
			hv = ((hv << 2) | symbolIdx[rc_read[i + hash_len_u]]);
		}
		pln = ht_u->find64_p(hv, rc_read, 0);
		if (pln != NULL)
			pnodes.insert(pln);
		pln = ht_d->find64_p(hv, rc_read, 0);
		if (pln != NULL)
			pnodes.insert(pln);

		/*if (genomes[3]->read_cnts_u > 0) {
			fprintf(stderr, "Begin with error 4.\n");
			abort();
		}*/

		/* Record results. */
		int i = 0;
		for (auto pn : pnodes) {
			if (pn->refID2 == 0)
				rids.insert(pn->refID1);
			else {
				if (pn->refID1 < pn->refID2)
					rid_pairs.insert(std::make_pair(pn->refID1, pn->refID2));
				else
					rid_pairs.insert(std::make_pair(pn->refID2, pn->refID1));
			}
		}
		
		/*if (genomes[3]->read_cnts_u > 0) {
			fprintf(stderr, "Mid with error.\n");
			abort();
		}*/

		//fprintf(stderr, "Read %lu; %lu\n", nrd, genomes.size());
		switch (rid_pairs.size()) {
			case 0:
				if (rids.empty())
					nundet++;
				else {
					if (rids.size() == 1) {
						for (auto rid : rids) {
							genomes[rid]->read_cnts_u++;
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
						genomes[rid_pair.first]->read_cnts_d++;
						genomes[rid_pair.second]->read_cnts_d++;
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
									genomes[rid]->read_cnts_u++;
									genomes[rid]->read_cnts_d++;
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
							genomes[rid]->read_cnts_u++;
							genomes[rid]->read_cnts_d++;
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
								genomes[rid]->read_cnts_d++;
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
		//if (genomes[3]->read_cnts_u > 0 && nrd < 1000)
		//	fprintf(stderr, "%lu\n", nrd);
		if (nrd++ % 100000 == 0) {
			//fprintf(stderr, "%u\t%lu\t%lu\n", genomes[3]->taxID, genomes[3]->read_cnts_u, genomes[3]->read_cnts_d);
			fprintf(stderr, "\rProcessed %lu reads.\n", nrd);
		}
	}

	fprintf(stderr, "Number of unlabeled reads: %lu.\n", nundet);
	fprintf(stderr, "Number of reads with conflict labels: %lu.\n", nconf);
	//delete []rc_read;
	//for (size_t i = 0; i < 4122; i++)
	//	fprintf(stderr, "%lu, %u\t%lu\t%lu\n", i, genomes[i + 1]->taxID, genomes[i + 1]->read_cnts_u, genomes[i + 1]->read_cnts_d);

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
	size_t n_species = genomes.size() - 1;
	fprintf(stderr, "%lu\n", n_species);
	exist = new int[n_species];
	memset(exist, 1, sizeof(int) * n_species);
	IloBoolVarArray EXIST(env, n_species);
	fprintf(stderr, "ILP environment in CPLEX initialized.\n");

	// Constraint 0: If the number of reads being assigned to a particular 
	//	genome g is less than READ_CNT_THRES, then EXIST[g] is set to 0
	for (size_t i = 1; i <= n_species; i++) {
		double d1 = genomes[i]->read_cnts_u * 1.0, d2 = genomes[i]->read_cnts_u * 1.0;
		d1 -= read_cnt_thres;
		d2 -= (1.0 * genomes[i]->nus) * resolution;
		//fprintf(stderr, "%u; %.4f; %.4f; %u; %u\n", rid, d1, d2, species_nus[rid], unique_thres);
		if (genomes[i]->nus >= unique_thres) {
			if (d1 < 0.0 || d2 < 0.0) {
				exist[i - 1] = 0;
				model.add(EXIST[i - 1] == 0);
			} 
		} else {
			if (d2 < 0.0) {
				exist[i - 1] = 0;
				model.add(EXIST[i - 1] == 0);
			}
		}
	}
	for (size_t i = 1; i <= n_species; i++) {
		float d1 = genomes[i]->read_cnts_d * 1.0, d2 = genomes[i]->read_cnts_d * 1.0;
		d1 -= read_cnt_thres;
		d2 -= (1.0 * genomes[i]->nds) * resolution;
		//fprintf(stderr, "%u; %.4f; %.4f\n", rid, d1, d2);
		if (genomes[i]->nus >= unique_thres) { // we should try nds?
			if (d1 < 0.0 || d2 < 0.0) {
				exist[i - 1] = 0;
				model.add(EXIST[i - 1] == 0);
			} 
		} else {
			if (d2 < 0.0) {
				exist[i - 1] = 0;
				model.add(EXIST[i - 1] == 0);
			} 
		}
	}

	size_t n_species_exist = 0;
	for (size_t i = 0; i < n_species; i++) {
		if (exist[i]) {
			n_species_exist++;
			fprintf(stderr, "Species %u may exist. \n", genomes[i + 1]->taxID);
		}
	}
	fprintf(stderr, "\tConstructed constraint 0.\n");
	//for (auto rc : read_cnts_u)
	//	if (exist[species_order[rc.first]])
	//		fprintf(stderr, "%u\t%lu\t%lu\n", rc.first, rc.second, read_cnts_d[rc.first]);
	for (size_t i = 0; i < n_species; i++) {
		if (exist[i]) {
			fprintf(stderr, "%lu\t%u\t%lu\t%lu\n\n", i, genomes[i + 1]->taxID, genomes[i + 1]->read_cnts_u, genomes[i + 1]->read_cnts_d);
		}
	}

	size_t num_us_valid = 0, abs_index = 0;
	for (size_t i = 0; i < n_species; i++) {
		if (exist[i] > 0) {
			//num_us_valid += ht_u->map_cnt[rid].first;
			//num_us_valid += ht_u->map_cnt[rid].second;
			num_us_valid += ht_u->map_sp[i + 1].size();
			num_us_valid += ht_d->map_sp[i + 1].size();
		}
	}
	fprintf(stderr, "Filtered out invalid species by read counts.\n");

	// Continuous variable C[l]: Coverage of species l.
	IloNumVarArray COV(env, n_species, 0.0, max_cov, ILOFLOAT);
	fprintf(stderr, "\tConstructed all other variables.\n");

	fprintf(stderr, "\tConstructed the objective function.\n");
	// Minimize the total number of species. 
	for (size_t i = 0; i < n_species; i++) {
		double u_factor = 1000.0 / ht_u->map_sp[i + 1].size();
		if (exist[i] > 0) {
			for (auto pn : ht_u->map_sp[i + 1]) {
				assert((i + 1 == pn->refID1) || (i + 1 == pn->refID2));
				double wcov = (pn->ucount1 * (rl - pn->depth) * 1.0 / rl);
				wcov = wcov * pow(1 - erate, pn->depth);
				objective += (wcov * COV[i] - pn->rcount) * (wcov * COV[i] - pn->rcount) * u_factor;
				abs_index++;
			} 
		}
	}
	for (size_t i = 0; i < n_species; i++) {
		double d_factor = 1000.0 / ht_d->map_sp[i + 1].size();
		if (exist[i] > 0) {
			for (auto pn : ht_d->map_sp[i + 1]) {
				assert((i + 1 == pn->refID1) || (i + 1 == pn->refID2));
				uint32_t si1 = pn->refID1 - 1, si2 = pn->refID2 - 1;
				double wcov1 = (pn->ucount1 * (rl - pn->depth) * 1.0 / rl);
				double wcov2 = (pn->ucount2 * (rl - pn->depth) * 1.0 / rl);
				wcov1 = wcov1 * pow(1 - erate, pn->depth);
				wcov2 = wcov2 * pow(1 - erate, pn->depth);
				objective += (wcov1 * COV[si1] + wcov2 * COV[si2] - pn->rcount) * 
						(wcov1 * COV[si1] + wcov2 * COV[si2] - pn->rcount) * d_factor;
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
	}
	
	IloExprArray EXP1(env, n_species_exist);
	size_t exp_i = 0;
	for (; exp_i < n_species_exist; exp_i++)
		EXP1[exp_i] = IloExpr(env);
	//fprintf(stderr, "n_species_exist: %lu\n", n_species_exist);
	exp_i = 0;
	for (size_t i = 0; i < n_species; i++) {
		if (exist[i] > 0) {
			for (auto pn : ht_u->map_sp[i + 1]) {
				double wcov = (pn->ucount1 * (rl - pn->depth) * 1.0 / rl);
				wcov = wcov * pow(1 - erate, pn->depth);
				EXP1[exp_i] += wcov * COV[i];
			} 
			if (genomes[i + 1]->nus >= unique_thres)
				model.add(EXP1[exp_i] * (1.0 + epsilon) - genomes[i + 1]->read_cnts_u >= 0);
			else
				model.add(EXP1[exp_i] >= 0.0);
			exp_i++;	
		}
	}
	IloExprArray EXP2(env, n_species_exist);
	for (exp_i = 0; exp_i < n_species_exist; exp_i++)
		EXP2[exp_i] = IloExpr(env);

	exp_i = 0;
	for (size_t i = 0; i < n_species; i++) {
		if (exist[i] > 0) {
			for (auto pn : ht_d->map_sp[i + 1]) {
				uint32_t si1 = pn->refID1 - 1, si2 = pn->refID2 - 1;
				double wcov1 = (pn->ucount1 * (rl - pn->depth) * 1.0 / rl);
				double wcov2 = (pn->ucount2 * (rl - pn->depth) * 1.0 / rl);
				wcov1 = wcov1 * pow(1 - erate, pn->depth);
				wcov2 = wcov2 * pow(1 - erate, pn->depth);
				EXP2[exp_i] += wcov1 * COV[si1] + wcov2 * COV[si2];
			}
			if (genomes[i + 1]->nus >= unique_thres)
				model.add(EXP2[exp_i] * (1.0 + epsilon) - genomes[i + 1]->read_cnts_d >= 0);
			else
				model.add(EXP2[exp_i] >= 0.0);
			exp_i++;
		}
	}
	fprintf(stderr, "\tConstructed constraint 3.\n");

	// Constraint 4: refining the search space
	//fprintf(stderr, "EXPR\n");
	IloExpr TOTAL(env);
	//fprintf(stderr, "EXPR\n");
	//fprintf(stderr, "%lu\n", n_species);
	for (size_t i = 0; i < n_species; i++) {
		//fprintf(stderr, "%lu, %u\n", i, genomes[i]->glength);
		TOTAL += (COV[i] * (1.0 * genomes[i + 1]->glength) / rl);
	}
	//model.add(TOTAL >= 0.9 * reads.size());
	model.add(TOTAL <= (1.0 + epsilon) * reads.size());
	fprintf(stderr, "\tConstructed constraint 4.\n");

	IloExpr e1(env);
	for (size_t i = 0; i < n_species; i++) {
		e1 += EXIST[i];
	}
	//model.add(e1 <= 2);
	//fprintf(stderr, "\tConstructed constraint 6.\n");
	model.add(e1 <= n_species_exist * 1.0);
	FILE * fout = fopen(OUTPUTFILE.c_str(), "w");

	// Solving the ILP.
	try {	
		IloCplex cplex(model);
		//cplex.setParam(IloCplex::IntParam::Threads, 32);
		cplex.setParam(IloCplex::NumParam::TiLim, 10800);
		cplex.setParam(IloCplex::NumParam::SolnPoolAGap, 0.0);

		if (cplex.populate()) {
			int nsolns = cplex.getSolnPoolNsolns();
			for (int s = 0; s < nsolns; s++) {
				//fprintf(fout, "\nk = %d; Objective = %.4f.\n", k, cplex.getObjValue(s));
				fprintf(fout, "\nSolution %d:\n", s);
				for (size_t i = 0; i < n_species; i++) {
					bool ei = cplex.getValue(EXIST[i], s);
					float ci = cplex.getValue(COV[i], s);
					if (ei != 0)
						fprintf(fout, "%u\t%d\t%.4f\n", genomes[i + 1]->taxID, ei, ci);
				}
			}
		}

		fclose(fout);
	} catch (IloException &ex) {}

	if (exist != NULL)
		delete [] exist;

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::high_resolution_clock::now() - start).count();
	fprintf(stderr, "Time for ilp: %lu ms.\n", duration);
}

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

