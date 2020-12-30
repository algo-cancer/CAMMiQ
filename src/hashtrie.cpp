#include <assert.h>
#include <chrono>

#include "hashtrie.hpp"
#include "util.hpp"
//#include "binaryio_temp.hpp"
#include "binaryio.hpp"

trieNode::trieNode() {
	isEnd = false;
	children = new trieNode*[4];
	for (int i = 0; i < 4; i++)
		children[i] = NULL;
}

trieNode::~trieNode() {
	if (children != NULL) {
		for (int i = 0; i < 4; i++)
			if (children[i] != NULL)
				delete children[i];
		delete [] children;
	}
}

leafNode::leafNode() {
	isEnd = true;
	/*children = new trieNode*[4];
	for (int i = 0; i < 4; i++)
		children[i] = NULL;*/
}

leafNode::leafNode(uint32_t refID_) {
	isEnd = true;
	refID = refID_;
	/*children = new trieNode*[4];
	for (int i = 0; i < 4; i++)
		children[i] = NULL;*/
}

dleafNode::dleafNode() {
	isEnd = true;
	refID = std::make_pair(0, 0);   
	ucount = std::make_pair(0, 0);
	/*children = new trieNode*[4];
	for (int i = 0; i < 4; i++)
		children[i] = NULL;*/   
}

void dleafNode::resetRefID() {
	refID.first = 0;
	refID.second = 0;
	ucount.first = 0;
	ucount.second = 0;
}

void dleafNode::setRefID_first(uint32_t refID_1, uint16_t uc1) {
	refID.first = refID_1;
	ucount.first = uc1;
}

void dleafNode::setRefID_second(uint32_t refID_2, uint16_t uc2) {
	refID.second = refID_2;
	ucount.second = uc2;
}

int dleafNode::compareRefID(uint32_t refID_1, uint16_t uc1, uint32_t refID_2, uint16_t uc2) {
	if (refID.first == refID_1) {
		if ((ucount.first != uc1) || (refID.second != refID_2) || (ucount.second != uc2))
			return 0;
		return 1;
	}
	if (refID.second == refID_1) {
		if ((ucount.second != uc1) || (refID.first != refID_2) || (ucount.first != uc2))
			return 0;
		return 2;
	}
	return 0;
}

pleafNode::pleafNode() {
	isEnd = true; 
}

void pleafNode::reset() {
	depth = 0;
	refID1 = 0;
	refID2 = 0;
	//next1 = NULL;
	//next2 = NULL;
	ucount1 = 0;
	ucount2 = 0;
	rcount = 0;
}

/*Hash::Hash(uint32_t hash_len, uint8_t d_) {
	hash_len_ = hash_len;
	uniq_ = d_;
}*/

Hash::Hash(uint32_t hash_len) {
	hash_len_ = hash_len;
}

Hash::~Hash() {
	if (!map32.empty()) {
		for (auto it : map32)
			delete it.second;
		map32.clear();
	}
	if (!map64.empty()) {
		for (auto it : map64)
			delete it.second;
		map64.clear();
	}
	//if (!smap32.empty())
	//	smap32.clear();
	//if (!smap64.empty())
	//	smap64.clear();
}

void Hash::clear() {
	if (!map32.empty()) {
		for (auto it : map32)
			delete it.second;
		map32.clear();
		fprintf(stderr, "HASHMAP: cleared 32-bits map.\n");
	}
	if (!map64.empty()) {
		for (auto it : map64)
			delete it.second;
		map64.clear();
		fprintf(stderr, "HASHMAP: cleared 64-bits map.\n");
	}
}

trieNode* Hash::getNode() {
	trieNode *pNode = new trieNode();
	return pNode;
}

leafNode* Hash::getLeaf() {
	leafNode *lNode = new leafNode();
	return lNode;
}

dleafNode* Hash::getLeaf_d() {
	dleafNode *dNode = new dleafNode();
	return dNode;
}

pleafNode* Hash::getLeaf_p() {
	pleafNode *pNode = new pleafNode();
	return pNode;
}

uint32_t Hash::computeHashVal(const uint8_t* key) {
	uint32_t res = 0;
	for (uint16_t i = 0; i < hash_len_; i++)
		res = ((res << 2) | symbolIdx[key[i]]);
	return res;
}

uint64_t Hash::computeHashVal64(const uint8_t* key) {
	uint64_t res = 0;
	for (uint16_t i = 0; i < hash_len_; i++)
		res = ((res << 2) | symbolIdx[key[i]]);
	return res;
}

uint32_t Hash::computeHashVal(const char* key) {
	uint32_t res = 0;
	for (uint16_t i = 0; i < hash_len_; i++)
		res = ((res << 2) | symbolIdx[(unsigned char) key[i]]);
	return res;
}

uint64_t Hash::computeHashVal64(const char* key) {
	uint64_t res = 0;
	for (uint16_t i = 0; i < hash_len_; i++)
		res = ((res << 2) | symbolIdx[(unsigned char) key[i]]);
	return res;
}

void abortInsert() {
	fprintf(stderr, "Illegal insertion, another key with the same prefix already exists.\n");
	abort();
}

/*
	Insert a unique substring (stroed in KEY) to hash table.
	KEY_LEN_ (>= HASH_LEN_) includes hashed part.
 */
//void Hash::insert32(uint32_t bucket_, const uint8_t* key, uint32_t key_len_, uint32_t refID_, uint8_t uc, uint32_t pos) {
void Hash::insert32(uint32_t bucket_, const uint8_t* key, uint32_t key_len_, uint32_t refID_, uint8_t uc) {
	HashMap32::iterator iter = map32.find(bucket_);
	leafNode *leaf = NULL;
	if (iter != map32.end()) {
		//fprintf(stderr, "111\n");
		if (key_len_ == hash_len_) {
			if (iter->second->isEnd) {
				leaf = (leafNode*)iter->second;
				if (refID_ == leaf->refID) {
					//leaf->ucount++;
					assert(uc == leaf->ucount);
					return;
				}
			}
			abortInsert();
		}
		//fprintf(stderr, "aa\n");
		leaf = insert(iter->second, key + hash_len_, key_len_ - hash_len_, refID_, uc);
	} else {
		//fprintf(stderr, "222\n");
		if (key_len_ == hash_len_) {
			leaf = getLeaf();
			//map32.insert(std::make_pair(bucket_, leaf));
			map32[bucket_] = leaf;
		} else {
			trieNode *root = getNode();
			//map32.insert(std::make_pair(bucket_, root));
			map32[bucket_] = root;
 			leaf = insert(root, key + hash_len_, key_len_ - hash_len_, refID_);
		}
	}
	if (leaf != NULL) {
		leaf->refID = refID_;
		leaf->ucount = uc;
		//leaf->pos = pos;
	}
}

void Hash::insert32_d(uint32_t bucket_, const uint8_t* key, 
			uint32_t key_len_, uint32_t refID_1, 
			uint8_t uc1, uint32_t refID_2, uint8_t uc2) {
	HashMap32::iterator iter = map32.find(bucket_);
	dleafNode *dleaf = NULL;
	if (iter != map32.end()) {
		if (key_len_ == hash_len_) {
			if (iter->second->isEnd) {
				dleaf = (dleafNode*)iter->second;
				if (dleaf->compareRefID(refID_1, uc1, refID_2, uc2) == 0) {
				//	abortInsert();
					fprintf(stderr, "Information not matching.\n");
					abort();
				}
				return;
			} else
				abortInsert();
		}
		dleaf = insert_d(iter->second, key + hash_len_, key_len_ - hash_len_, refID_1, uc1, refID_2, uc2);
	} else {
		if (key_len_ == hash_len_) {
			dleaf = getLeaf_d();
			dleaf->resetRefID();
			//map32.insert(std::make_pair(bucket_, dleaf));
			map32[bucket_] = dleaf;

		} else {
			trieNode *root = getNode();
			//map32.insert(std::make_pair(bucket_, root));
			map32[bucket_] = root;
 			dleaf = insert_d(root, key + hash_len_, key_len_ - hash_len_);
		}
	}
	if ((dleaf != NULL) && (dleaf->refID.first == 0)) { 
		dleaf->setRefID_first(refID_1, uc1);
		dleaf->setRefID_second(refID_2, uc2);
		//dleaf->pos = pos;
	}
}

//void Hash::insert64(uint64_t bucket_, const uint8_t* key, uint32_t key_len_, uint32_t refID_, uint8_t uc, uint32_t pos) {
void Hash::insert64(uint64_t bucket_, const uint8_t* key, uint32_t key_len_, uint32_t refID_, uint8_t uc) {
	HashMap64::iterator iter = map64.find(bucket_);
	leafNode *leaf = NULL;
	if (iter != map64.end()) {
		if (key_len_ == hash_len_) {
			if (iter->second->isEnd) {
				leaf = (leafNode*)iter->second;
				if (refID_ == leaf->refID) {
					assert(uc == leaf->ucount);
					return;
				}
			}
			abortInsert();
		}
		leaf = insert(iter->second, key + hash_len_, key_len_ - hash_len_, refID_, uc);
	} else {
		if (key_len_ == hash_len_) {
			leaf = getLeaf();
			//map64.insert(std::make_pair(bucket_, leaf));
			map64[bucket_] = leaf;
		} else {
			trieNode *root = getNode();
			//map64.insert(std::make_pair(bucket_, root));
			map64[bucket_] = root;
 			leaf = insert(root, key + hash_len_, key_len_ - hash_len_, refID_);
		}
	}
	if (leaf != NULL) {
		leaf->refID = refID_;
		leaf->ucount = uc;
		//leaf->pos = pos;
	}	
}

//void Hash::insert64_d(uint64_t bucket_, const uint8_t* key, uint32_t key_len_, uint32_t refID_1,
//			uint8_t uc1, uint32_t refID_2, uint8_t uc2, uint32_t pos) {
void Hash::insert64_d(uint64_t bucket_, const uint8_t* key, uint32_t key_len_, uint32_t refID_1,
			uint8_t uc1, uint32_t refID_2, uint8_t uc2) {
	HashMap64::iterator iter = map64.find(bucket_);
	dleafNode *dleaf = NULL;
	if (iter != map64.end()) {
		if (key_len_ == hash_len_) {
			if (iter->second->isEnd) {
				dleaf = (dleafNode*)iter->second;
				/*int cr = dleaf->compareRefID(refID_1, uc1, refID_2, uc2);
				switch (cr) {
					case 1:
						if (dleaf->pos.first == 0)
							dleaf->pos.first = pos;
						break;
					case 2:
						if (dleaf->pos.second == 0)
							dleaf->pos.second = pos;
						break;
					default:
						fprintf(stderr, "Information not matching.\n");
						abort();
				*/
				if (dleaf->compareRefID(refID_1, uc1, refID_2, uc2) == 0) {
					fprintf(stderr, "%u\t%u\t%u\t%u\t", dleaf->refID.first, dleaf->ucount.first, dleaf->refID.second, dleaf->ucount.second);
					fprintf(stderr, "%u\t%u\t%u\t%u\t", refID_1, uc1, refID_2, uc2);
					fprintf(stderr, "Information not matching.\n");
					abort();
				}
				//}
				return;
			} else
				abortInsert();
		}
		dleaf = insert_d(iter->second, key + hash_len_, key_len_ - hash_len_, refID_1, uc1, refID_2, uc2);
	} else {
		if (key_len_ == hash_len_) {
			dleaf = getLeaf_d();
			dleaf->resetRefID();
			//map64.insert(std::make_pair(bucket_, dleaf));
			map64[bucket_] = dleaf;
		} else {
			trieNode *root = getNode();
			//map64.insert(std::make_pair(bucket_, root));
			map64[bucket_] = root;
 			dleaf = insert_d(root, key + hash_len_, key_len_ - hash_len_);
		}
	}
	if (dleaf != NULL) {
		if (dleaf->refID.first == 0) {
			dleaf->setRefID_first(refID_1, uc1);
			dleaf->setRefID_second(refID_2, uc2);
			//dleaf->pos.first = 0;
			//dleaf->pos.second = 0;
		}
		/*int cr = dleaf->compareRefID(refID_1, uc1, refID_2, uc2);
		switch (cr) {
			case 1:
				if (dleaf->pos.first == 0)
					dleaf->pos.first = pos;
				break;
			case 2:
				if (dleaf->pos.second == 0)
					dleaf->pos.second = pos;
				break;
			default:
				break;
		}*/
	}
}

/*
	The actual insert (helper) to a trie structure.
	KEY_LEN_ here (> 0) does not include hashed part
	Observe that none of the "key" can be a prefix of any other "key". 
 */
leafNode* Hash::insert(trieNode *root, const uint8_t* key, uint32_t key_len_, uint32_t refID_, uint8_t uc) {
	//fprintf(stderr, "%u, %u\n", key_len_, refID_);
	trieNode *cur = root;
	int index;
	for (uint32_t i = 0; i < key_len_ - 1; i++) {
        	index = symbolIdx[key[i]];
		if (!cur->children[index])
			cur->children[index] = getNode();
		cur = cur->children[index];
		if (cur->isEnd)
			abortInsert();
	}
	//fprintf(stderr, "%u, %u\n", key_len_, refID_);
	/* Mark last node as leaf and set its label. */
	index = symbolIdx[key[key_len_ - 1]];
	if (!cur->children[index]) {
		cur->children[index] = getLeaf();
		return (leafNode*) cur->children[index];
	} else {
		if (cur->children[index]->isEnd) {	
			leafNode *leaf = (leafNode*)cur->children[index];
			if (refID_ == leaf->refID) {
				//leaf->ucount++;
				assert(leaf->ucount == uc);
				return NULL;
			}
		}
		fprintf(stderr, "INDEX %d\n", index);
		abortInsert();
	}
	return NULL;
}

leafNode* Hash::insert(trieNode *root, const uint8_t* key, uint32_t key_len_, uint32_t refID_) {
	//fprintf(stderr, "%u, %u\n", key_len_, refID_);
	trieNode *cur = root;
	int index;
	for (uint32_t i = 0; i < key_len_ - 1; i++) {
        	index = symbolIdx[key[i]];
		if (!cur->children[index])
			cur->children[index] = getNode();
		cur = cur->children[index];
		if (cur->isEnd)
			abortInsert();
	}
	//fprintf(stderr, "%u, %u\n", key_len_, refID_);
	/* Mark last node as leaf and set its label. */
	index = symbolIdx[key[key_len_ - 1]];
	if (!cur->children[index]) {
		cur->children[index] = getLeaf();
		return (leafNode*) cur->children[index];
	} else {
		if (cur->children[index]->isEnd) {	
			leafNode *leaf = (leafNode*)cur->children[index];
			if (refID_ == leaf->refID) {
				//leaf->ucount++;
				//assert(leaf->ucount == uc);
				return NULL;
			}
		}
		fprintf(stderr, "INDEX %d\n", index);
		abortInsert();
	}
	return NULL;
}

// Update 2019/07/27
// When first insert we do not need refID
dleafNode* Hash::insert_d(trieNode *root, const uint8_t* key, uint32_t key_len_) {
	trieNode *cur = root;
	int index;
	for (uint32_t i = 0; i < key_len_ - 1; i++) {
        	index = symbolIdx[key[i]];
		if (!cur->children[index])
			cur->children[index] = getNode();
		cur = cur->children[index];
		if (cur->isEnd) {
			abortInsert();
		}
	}
	/* Mark last node as leaf and set its label. */
	index = symbolIdx[key[key_len_ - 1]];
	dleafNode *dleaf;
	if (!cur->children[index]) {
		cur->children[index] = getLeaf_d();
		dleaf = ((dleafNode*) cur->children[index]);
		//dleaf->resetRefID_first(refID_);
		dleaf->resetRefID();
		return dleaf;
	} else {
		if (!cur->children[index]->isEnd) {
			abortInsert();
		}
	}
	return NULL;
}

dleafNode* Hash::insert_d(trieNode *root, const uint8_t* key, uint32_t key_len_,
	uint32_t refID_1, uint8_t uc1, uint32_t refID_2, uint8_t uc2) {
	trieNode *cur = root;
	int index;
	for (uint32_t i = 0; i < key_len_ - 1; i++) {
        	index = symbolIdx[key[i]];
		if (!cur->children[index])
			cur->children[index] = getNode();
		cur = cur->children[index];
		if (cur->isEnd) {
			abortInsert();
		}
	}
	/* Mark last node as leaf and set its label. */
	index = symbolIdx[key[key_len_ - 1]];
	dleafNode *dleaf;
	if (!cur->children[index]) {
		cur->children[index] = getLeaf_d();
		dleaf = ((dleafNode*) cur->children[index]);
		//dleaf->resetRefID_first(refID_);
		dleaf->resetRefID();
		return dleaf;
	} else {
		if (!cur->children[index]->isEnd) {
			dleaf = (dleafNode*)cur->children[index];
			if (dleaf->compareRefID(refID_1, uc1, refID_2, uc2))
				return NULL;
			else
				abortInsert();
		}
	}
	return NULL;
}


/*
	Matching substring CAND within BUCKET_ and maximum length LEN.
 */
uint32_t Hash::find32(uint32_t bucket_, const char* cand, size_t len_) {
	HashMap32::iterator iter = map32.find(bucket_);
	if (iter == map32.end())
		return 0;
	else {
		trieNode *cur = iter->second;
		int index;
		for (size_t i = 0; i < len_; i++) {
			index = symbolIdx[(unsigned char) cand[i]];
			if (cur->isEnd)
				return ((leafNode*) cur)->refID;
			if (cur->children[index] == NULL)
				return 0;
			cur = cur->children[index];
		}
		if (cur != NULL && cur->isEnd)
			return ((leafNode*) cur)->refID;
	}
	return 0;
}


uint32_t Hash::find32(uint32_t bucket_, const uint8_t* cand, size_t len_) {
	
	//fprintf(stderr, "%u; %lu\n", bucket_, len_);
	HashMap32::iterator iter = map32.find(bucket_);
	if (iter == map32.end())
		return 0;
	else {
		trieNode *cur = iter->second;
		int index;
		for (size_t i = 0; i < len_; i++) {
			index = symbolIdx[cand[i]];
			//fprintf(stderr, "A: %lu, %d.\n", len_, index);
			if (cur->isEnd) {
				//fprintf(stderr, "%u\n", bucket_);

				//for (size_t j = 0; j < len_ + hash_len_; j++)
				//	fprintf(stderr, "%c", (char) cand[j - hash_len_]);
				//fprintf(stderr, "\n");
				return ((leafNode*) cur)->refID;
			}
			if (cur->children[index] == NULL) {
				//fprintf(stderr, "C\n");
				return 0;
			}
			cur = cur->children[index];
		}
		if (cur != NULL && cur->isEnd) {
			//fprintf(stderr, "B\n");
			return ((leafNode*) cur)->refID;
		}
	}
	return 0;
}

pleafNode* Hash::find32_p(uint32_t bucket_, const uint8_t* cand, size_t len_) {
	HashMap32::iterator iter = map32.find(bucket_);
	if (iter == map32.end())
		return NULL;
	else {
		trieNode *cur = iter->second;
		int index;
		for (size_t i = 0; i < len_; i++) {
			index = symbolIdx[cand[i]];
			if (cur->isEnd)
				return ((pleafNode*) cur);
			if (cur->children[index] == NULL)
				return NULL;
			cur = cur->children[index];
		}
		if (cur != NULL && cur->isEnd)
			return ((pleafNode*) cur);
	}
	return 0;
}

pleafNode* Hash::find64_p(uint64_t bucket_, const uint8_t* cand, size_t len_) {
	HashMap64::iterator iter = map64.find(bucket_);
	if (iter == map64.end())
		return NULL;
	else {
		trieNode *cur = iter->second;
		int index;
		for (size_t i = 0; i < len_; i++) {
			index = symbolIdx[cand[i]];
			if (cur->isEnd)
				return ((pleafNode*) cur);
			if (cur->children[index] == NULL)
				return NULL;
			cur = cur->children[index];
		}
		if (cur != NULL && cur->isEnd)
			return ((pleafNode*) cur);
	}
	return 0;
}

uint32_t Hash::find32_test(uint32_t bucket_, const uint8_t* cand, size_t len_) {
	fprintf(stderr, "%u; %lu\n", bucket_, len_);
	fprintf(stderr, "AAAA\n");
	HashMap32::iterator iter = map32.find(bucket_);
	fprintf(stderr, "AAAA\n");
	if (iter == map32.end()) {
		fprintf(stderr, "BBBB\n");
		return 0;
	}
	else {
		fprintf(stderr, "CCCC\n");

		trieNode *cur = iter->second;
		fprintf(stderr, "DDDD\n");
		int index;
		for (size_t i = 0; i < len_; i++) {
			fprintf(stderr, "%u; %lu\n", bucket_, i);
			index = symbolIdx[cand[i]];
			if (cur->isEnd) {
				fprintf(stderr, "%u.\n", ((leafNode*) cur)->refID);
				return ((leafNode*) cur)->refID;
			}
			if (cur->children[index] == NULL)
				return 0;
			cur = cur->children[index];
		}
		if (cur != NULL && cur->isEnd)
			return ((leafNode*) cur)->refID;
	}
	return 0;
}

std::pair<uint32_t, uint32_t> Hash::find32_d(uint32_t bucket_, const char* cand, size_t len_) {
	HashMap32::iterator iter = map32.find(bucket_);
	if (iter == map32.end())
		return std::make_pair(0, 0);
	else {
		trieNode *cur = iter->second;
		int index;
		for (size_t i = 0; i < len_; i++) {
			index = symbolIdx[(unsigned char) cand[i]];
			if (cur->isEnd)
				return ((dleafNode*) cur)->refID;
			if (cur->children[index] == NULL)
				return std::make_pair(0, 0);
			cur = cur->children[index];
		}
		if (cur != NULL && cur->isEnd)
			return ((dleafNode*) cur)->refID;
	}
	return std::make_pair(0, 0);
}

std::pair<uint32_t, uint32_t> Hash::find32_d(uint32_t bucket_, const uint8_t* cand, size_t len_) {
	HashMap32::iterator iter = map32.find(bucket_);
	if (iter == map32.end())
		return std::make_pair(0, 0);
	else {
		trieNode *cur = iter->second;
		int index;
		for (size_t i = 0; i < len_; i++) {
			index = symbolIdx[cand[i]];
			if (cur->isEnd) {
				//fprintf(stderr, "%lu\n", i);
				return ((dleafNode*) cur)->refID;
			}
			if (cur->children[index] == NULL)
				return std::make_pair(0, 0);
			cur = cur->children[index];
		}
		if (cur != NULL && cur->isEnd) {
			//fprintf(stderr, "%lu\n", len_);

			return ((dleafNode*) cur)->refID;}
	}
	return std::make_pair(0, 0);
}

std::vector<uint32_t> Hash::find_rids(uint32_t bucket_, const uint8_t* cand, size_t len_) {
	std::vector<uint32_t> result;
	HashMap32::iterator iter = map32.find(bucket_);
	if (iter == map32.end())
		return result;
	else {
		trieNode *cur = iter->second;
		int index;
		for (size_t i = 0; i < len_ - hash_len_; i++) {
			index = symbolIdx[(cand + hash_len_)[i]];
			/* Illegal case: one (doubly) unique substring is a superstring of another. */
			if (cur->isEnd)
				abort();
			/* Possible case: two (doubly) unique substrings share the same prefix. */
			if (cur->children[index] == NULL)
				return result;
			cur = cur->children[index];
		}
		if (cur == NULL)
			return result;
		if (cur->isEnd) {
			if (((dleafNode*)cur)->refID.first != 0)
				result.push_back(((dleafNode*) cur)->refID.first);
			if (((dleafNode*)cur)->refID.second != 0)
				result.push_back(((dleafNode*)cur)->refID.second);
			return result;
		} else {
			abort();
		}
	}
	return result;
}

uint32_t Hash::find_rids_u(uint32_t bucket_, const uint8_t* cand, size_t len_) {
	HashMap32::iterator iter = map32.find(bucket_);
	if (iter == map32.end())
		return 0;
	else {
		trieNode *cur = iter->second;
		int index;
		for (size_t i = 0; i < len_ - hash_len_; i++) {
			index = symbolIdx[(cand + hash_len_)[i]];
			/* Illegal case: one (doubly) unique substring is a superstring of another. */
			if (cur->isEnd)
				abort();
			/* Possible case: two (doubly) unique substrings share the same prefix. */
			if (cur->children[index] == NULL)
				return 0;
			cur = cur->children[index];
		}
		if (cur == NULL)
			return 0;
		if (cur->isEnd) {
			return ((leafNode*)cur)->refID;
		} else {
			abort();
		}
	}
	return 0;
}

int32_t Hash::find_and_increase_u(uint32_t bucket_, const uint8_t* cand, size_t len_) {
	HashMap32::iterator iter = map32.find(bucket_);
	if (iter == map32.end())
		return 0;
	else {
		trieNode *cur = iter->second;
		int index;
		for (size_t i = 0; i < len_ - hash_len_; i++) {
			index = symbolIdx[(cand + hash_len_)[i]];
			/* Illegal case: one (doubly) unique substring is a superstring of another. */
			if (cur->isEnd) {
				//fprintf(stderr, "aaa.\n");
				return -1;
			}
			/* Possible case: two (doubly) unique substrings share the same prefix. */
			if (cur->children[index] == NULL)
				return 0;
			cur = cur->children[index];
		}
		if (cur == NULL)
			return 0;
		if (cur->isEnd) {
			((leafNode*) cur)->ucount++;
			return 1;
		} else
			return -1;
	}
	return 0;
}

uint32_t Hash::find64(uint64_t bucket_, const char* cand, size_t len_) {
	HashMap64::iterator iter = map64.find(bucket_);
	if (iter == map64.end())
		return 0;
	else {
		trieNode *cur = iter->second;
		int index;
		for (size_t i = 0; i < len_; i++) {
			index = symbolIdx[(unsigned char) cand[i]];
			if (cur->isEnd)
				return ((leafNode*) cur)->refID;
			if (cur->children[index] == NULL)
				return 0;
			cur = cur->children[index];
		}
		if (cur != NULL && cur->isEnd)
			return ((leafNode*) cur)->refID;
	}
	return 0;
}

/*
	Load encoded index file (*.IDX)
 */
trieNode* Hash::decodeTrie() {
	int b = idx_reader->readBit();
	if (b == 0)
		return NULL;
	else {
		trieNode *root = getNode();
		bool isleaf = true;
		for (int i = 0; i < 4; i++) {
			trieNode *child = decodeTrie();
			root->children[i] = child;
			if (child != NULL)
				isleaf = false;
		}
		if (isleaf) {
			leafNode *leaf = getLeaf();
			leaf->refID = idx_reader->readBits32();
			leaf->ucount = idx_reader->readBits16();
			delete root;
			return leaf;
		}
		return root;
	}
}

trieNode* Hash::decodeTrie_d() {
	int b = idx_reader->readBit();
	if (b == 0)
		return NULL;
	else {
		trieNode *root = getNode();
		bool isdleaf = true;
		for (int i = 0; i < 4; i++) {
			trieNode *child = decodeTrie_d();
			root->children[i] = child;
			if (child != NULL)
				isdleaf = false;
		}
		if (isdleaf) {
			dleafNode *dleaf = getLeaf_d();
			dleaf->refID.first = idx_reader->readBits32();
			dleaf->refID.second = idx_reader->readBits32();
			dleaf->ucount.first = idx_reader->readBits16();
			dleaf->ucount.first = idx_reader->readBits16();
			delete root;
			return dleaf;
		}
		return root;
	}
}

/*trieNode* Hash::decodeTrie_d(uint8_t depth) {
	int b = idx_reader->readBits32(1);
	if (b == 0)
		return NULL;
	else {
		trieNode *root = getNode();
		bool isdleaf = true;
		for (int i = 0; i < 4; i++) {
			trieNode *child = decodeTrie_d(depth + 1);
			root->children[i] = child;
			if (child != NULL)
				isdleaf = false;
		}
		if (isdleaf) {
			dleafNode *dleaf = getLeaf_d();
			dleaf->refID.first = idx_reader->readBits32(32);
			dleaf->refID.second = idx_reader->readBits32(32);
			dleaf->ucount.first = idx_reader->readBits32(16);
			dleaf->ucount.second = idx_reader->readBits32(16);
			if (dleaf->refID.second != 0)
				if (expected_u2.find(dleaf->refID.first) == expected_u2.end()) {
					std::unordered_map<uint32_t, float> eu2_s1;
					eu2_s1[dleaf->refID.second] = (100 - depth - hash_len_) * 0.1 * ;
					expected_u2[dleaf->refID.first] = eu2_s1;
				} else {
					if (expected_u2[dleaf->refID.first].find(dleaf->refID.second) == expected_u2[dleaf->refID.first].end())
						expected_u2[dleaf->refID.first][dleaf->refID.second] = (100 - depth - hash_len_) * 0.1;
					else
						expected_u2[dleaf->refID.first][dleaf->refID.second] += (100 - depth - hash_len_) * 0.1;
				}
			delete root;
			return dleaf;
		}
		return root;
	}
}*/

trieNode* Hash::decodeTrie_p(int d_flag, uint8_t depth) {
	//fprintf(stderr, "AAAA.\n");
	int b = idx_reader->readBit();
	if (b == 0)
		return NULL;
	else {
		trieNode *root = getNode();
		if (d_flag) {
			bool isdleaf = true;
			for (int i = 0; i < 4; i++) {
				trieNode *child = decodeTrie_p(d_flag, depth + 1);
				root->children[i] = child;
				if (child != NULL)
					isdleaf = false;
			}
			if (isdleaf) {
				leaf_cnt++;
				pleafNode *pleaf = getLeaf_p();
				pleaf->depth = depth + hash_len_;
				uint32_t rid1 = idx_reader->readBits32();
				uint32_t rid2 = idx_reader->readBits32();
				/* */
				assert(rid1 != 0 && rid2 != 0);
				pleaf->refID1 = rid1;
				pleaf->refID2 = rid2;
				pleaf->ucount1 = idx_reader->readBits16();
				pleaf->ucount2 = idx_reader->readBits16();
				pleaf->rcount = 0;
				//uint32_t pos1 = idx_reader->readBits32();
				//uint32_t pos2 = idx_reader->readBits32();
				//if (pos1 == 0)
					map_sp[rid1].push_back(pleaf);
				//else
				//	map_sp[rid1][pos1 - 1] = pleaf;
				//if (pos2 == 0)
					map_sp[rid2].push_back(pleaf);
				//else
				//	map_sp[rid2][pos2 - 1] = pleaf;

				delete root;
				return pleaf;
			}
			return root;
		} else {
			//fprintf(stderr, "BBBB.\n");
			bool isleaf = true;
			for (int i = 0; i < 4; i++) {
				trieNode *child = decodeTrie_p(d_flag, depth + 1);
				root->children[i] = child;
				if (child != NULL)
					isleaf = false;
			}
			if (isleaf) {
				leaf_cnt++;
				pleafNode *pleaf = getLeaf_p();
				pleaf->depth = depth + hash_len_;
				uint32_t rid = idx_reader->readBits32();
				pleaf->refID1 = rid;
				pleaf->refID2 = 0;
				pleaf->ucount1 = idx_reader->readBits16();
				pleaf->rcount = 0;
				//uint32_t pos = idx_reader->readBits32();
				//if (pos == 0)
					map_sp[rid].push_back(pleaf);
				//else
				//	map_sp[rid][pos - 1] = pleaf;

				delete root;
				return pleaf;
			}
			return root;
		}
	}
}

/*trieNode* Hash::mergeTrie_p(trieNode* t1, trieNode* t2) {
	if (t1 == NULL)
		return t2;
	if (t2 == NULL)
		return t1;
	for (int i = 0; i < 4; i++)
        	t1->children[i] = mergeTrie_p(t1->children[i], t2->children[i]);
	return t1;
}*/

void Hash::loadIdx32(std::string &fn) {
	idx_reader = new BitReader();
	idx_reader->openFile(fn);
	int doubly_unique = idx_reader->readBit();
	assert(doubly_unique == 0);
	int option = idx_reader->readBits(7);
	assert(option == 32);
	hash_len_ = idx_reader->readBits(8);
	uint32_t bucket = idx_reader->readBits32();
	while (bucket != END32) {
		trieNode *root = decodeTrie();
		map32[bucket] = root;
		bucket = idx_reader->readBits32();
	}
	idx_reader->closeFile();
	fprintf(stderr, "Loaded index from file.\n");
	delete idx_reader;
}

void Hash::loadIdx32_d(std::string &fn) {
	idx_reader = new BitReader();
	idx_reader->openFile(fn);
	int doubly_unique = idx_reader->readBit();
	assert(doubly_unique == 1);
	int option = idx_reader->readBits(7);
	assert(option == 32);
	hash_len_ = idx_reader->readBits(8);
	uint32_t bucket = idx_reader->readBits32();
	while (bucket != END32) {
		trieNode *root = decodeTrie_d();
		map32[bucket] = root;
		bucket = idx_reader->readBits32();
	}
	idx_reader->closeFile();
	fprintf(stderr, "Loaded index from file.\n");
	delete idx_reader;
}

void Hash::loadIdx32_p(std::string &fn) {
	idx_reader = new BitReader();
	idx_reader->openFile(fn);
	//idx_reader->openFile("index_u2.bin2");
	int doubly_unique = idx_reader->readBit();
	//print()
	//assert(doubly_unique == 1);
	int option = idx_reader->readBits(7);
	assert(option == 32);
	hash_len_ = idx_reader->readBits(8);
	uint32_t bucket = idx_reader->readBits32();
	//int fff = 0;
	fprintf(stderr, "N: %s\n", fn.c_str());
	fprintf(stderr, "L: %d\n", hash_len_);
	while (bucket != END32) {
		//if (fff < 2) {
		//	fprintf(stderr, "%u\n", bucket); 
		//} 
		//fff++;
		//if (fff >= 100) break;
		//break;
		trieNode *root = decodeTrie_p(doubly_unique, 0);
		//trieNode *root = decodeTrie();
		map32[bucket] = root;
		bucket = idx_reader->readBits32();
	}
	idx_reader->closeFile();
	fprintf(stderr, "Loaded index from file.\n");
	delete idx_reader;
}

void Hash::loadIdx64_p(std::string &fn) {
	//auto start = std::chrono::high_resolution_clock::now();

	idx_reader = new BitReader();
	idx_reader->openFile(fn);
	int doubly_unique = idx_reader->readBit();
	//size_t count1 = 0, count2 = 0;
	
	int option = idx_reader->readBits(7);
	assert(option == 64);
	hash_len_ = idx_reader->readBits(8);
	uint64_t bucket = idx_reader->readBits64();
	
	fprintf(stderr, "Index: %s\n", fn.c_str());
	fprintf(stderr, "Hash Length: %d\n", hash_len_);
	while (bucket != END64) {
		/*1168575096146679
		0
		4047904091797938
		1294735211906458
		1963562052986519
		2610839374132841
		1819508562565474
		4179598439307609
		1862942747953436
		1294438259600985*/

		trieNode *root = decodeTrie_p(doubly_unique, 0);	
		//fprintf(stderr, "bucket = %lu\n", bucket);
		
		
		//if (root->isEnd) {
		//	count1++;
			//fprintf(stderr, "RefID1 = %u; RefID2 = %u\n", ((pleafNode *) root) -> refID1, ((pleafNode *) root) -> refID2);
			//if (bucket > 4503599627370496)
			//abort();
		//} else
		//	count2++;
		map64[bucket] = root;
		bucket = idx_reader->readBits64();
	}
	idx_reader->closeFile();
	//fprintf(stderr, "Loaded index from file.\n");
	//fprintf(stderr, "Num buckets: %lu; %lu\n", count1, count2);
	delete idx_reader;
	//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
	//		(std::chrono::high_resolution_clock::now() - start).count();
	//fprintf(stderr, "Time for query: %lu ms.\n", duration);
}

void Hash::loadIdx64_test(std::string &fn) {
	idx_reader = new BitReader();

	size_t ext_pos = fn.find("bin");
	std::string fn_main, fn_extra;
	fn_main = fn.substr(0, ext_pos - 1) + "_main." + fn.substr(ext_pos);
	fn_extra = fn.substr(0, ext_pos - 1) + "_extra." + fn.substr(ext_pos);

	

	
	idx_reader->openFile(fn_main);
	int doubly_unique = idx_reader->readBit();
	int option = idx_reader->readBits(7);
	fprintf(stderr, "%d, %d\n", option, doubly_unique);
	assert(option == 64);
	fprintf(stderr, "%d, %d\n", idx_reader->curBits, idx_reader->curByte);
	uint64_t bucket = idx_reader->readBits64();
	int indicator = 0;

	/* correct vals:
	2482328851543654
	67053043997664
	4318113957499127
	703548341397820
	2443053312351933
	2747081960180816
	3477702105075717
	2671024462725259
	4227154730525565
	1707883066870428
	
	51392913212908
	3922311404835771
	2904333816863000
	1520583437555794
	2839361992221291
	807317021782211
	1747703315486937
	3436626172995151
	2421924265360233
	2648242704181452*/
	
	while (bucket != END64) {
		//if (indicator < 10) {
		//	fprintf(stderr, "%lu\n", bucket);
		//	indicator++;
		//}	
		pleafNode *pleaf = getLeaf_p();
		pleaf->depth = 26;
		if (doubly_unique) {
			//uint32_t rid1 = idx_reader->read_uint32_t();
			//uint32_t rid2 = idx_reader->read_uint32_t();
			uint32_t rid1 = idx_reader->readBits32();
			uint32_t rid2 = idx_reader->readBits32();
			assert(rid1 != 0 && rid2 != 0);
			pleaf->refID1 = rid1;
			pleaf->refID2 = rid2;
			//pleaf->ucount1 = idx_reader->read_uint16_t();
			//pleaf->ucount2 = idx_reader->read_uint16_t();
			pleaf->ucount1 = idx_reader->readBits16();
			pleaf->ucount2 = idx_reader->readBits16();
			pleaf->rcount = 0;
			map_sp[rid1].push_back(pleaf);
			map_sp[rid2].push_back(pleaf);	
		} else {
			//uint32_t rid = idx_reader->read_uint32_t();
			uint32_t rid = idx_reader->readBits32();
			pleaf->refID1 = rid;
			pleaf->refID2 = 0;
			//pleaf->ucount1 = idx_reader->read_uint16_t();
			pleaf->ucount1 = idx_reader->readBits16();
			pleaf->rcount = 0;
			map_sp[rid].push_back(pleaf);
		}
		map64[bucket] = (trieNode*) pleaf;
		bucket = idx_reader->readBits64();
	}
	idx_reader->closeFile();
	fprintf(stderr, "Loaded main index from file.\n");

	idx_reader->reset();

	idx_reader->openFile(fn_extra);
	doubly_unique = idx_reader->readBit();
	
	option = idx_reader->readBits(7);
	assert(option == 64);
	hash_len_ = idx_reader->readBits(8);
	bucket = idx_reader->readBits64();
	fprintf(stderr, "N: %s\n", fn_extra.c_str());
	fprintf(stderr, "L: %d\n", hash_len_);
	while (bucket != END64) {	
		trieNode *root = decodeTrie_p(doubly_unique, 0);
		map64[bucket] = root;
		bucket = idx_reader->readBits64();
	}
	idx_reader->closeFile();
	fprintf(stderr, "Loaded extra index from file.\n");

	//fprintf(stderr, "%lu; %lu", count1, count2);
	delete idx_reader;
}

/*void Hash::loadIdx64_b(std::string &fn_u, std::string &fn_d) {
	fprintf(stderr, "Size of node: %lu\n", sizeof(trieNode));
	//fprintf(stderr, "Size of node: %lu\n", sizeof(trieNodePtr));
	fprintf(stderr, "Size of leaf node: %lu\n", sizeof(leafNode));
	fprintf(stderr, "Size of pleaf node: %lu\n", sizeof(pleafNode));

	idx_reader = new BitReader();
	idx_reader->openFile(fn_u);
	int doubly_unique = idx_reader->readBits32(1);
	
	int option = idx_reader->readBits32(7);
	assert(option == 64);
	hash_len_ = idx_reader->readBits32(8);
	uint64_t bucket = idx_reader->readBits64(64);
	
	fprintf(stderr, "N: %s\n", fn_u.c_str());
	fprintf(stderr, "L: %d\n", hash_len_);
	while (bucket != END64) {	
		//fprintf(stderr, "bucket = %lu\n", bucket);
		trieNode *root = decodeTrie_p(doubly_unique, 0);
		map64[bucket] = root;
		bucket = idx_reader->readBits64(64);
	}
	idx_reader->closeFile();
	fprintf(stderr, "Loaded index from file.\n");

	idx_reader->reset();
	idx_reader->openFile(fn_d);
	doubly_unique = idx_reader->readBits32(1);
	fprintf(stderr, "N: %s\n", fn_d.c_str());

	fprintf(stderr, "N: %d\n", doubly_unique);

	option = idx_reader->readBits32(7);
	fprintf(stderr, "N: %d\n", option);
	assert(option == 64);
	hash_len_ = idx_reader->readBits32(8);
	bucket = idx_reader->readBits64(64);
	
	fprintf(stderr, "N: %s\n", fn_d.c_str());
	fprintf(stderr, "L: %d\n", hash_len_);
	while (bucket != END64) {	
		//fprintf(stderr, "bucket = %lu\n", bucket);
		trieNode *root = decodeTrie_p(doubly_unique, 0);
		HashMap64::iterator it64 = map64.find(bucket);
		if (it64 == map64.end())
			map64[bucket] = root;
		else {
			map64[bucket] = mergeTrie_p(map64[bucket], root);
			delete root;
		}
		bucket = idx_reader->readBits64(64);
	}
	idx_reader->closeFile();
	fprintf(stderr, "Loaded index from file.\n");
	delete idx_reader;
}*/

void Hash::loadIdx64(std::string &fn) {
	idx_reader = new BitReader();
	idx_reader->openFile(fn);
	int doubly_unique = idx_reader->readBit();
	assert(doubly_unique == 0);
	int option = idx_reader->readBits(7);
	assert(option == 64);
	hash_len_ = idx_reader->readBits(8);
	uint64_t bucket = idx_reader->readBits64();
	while (bucket != END64) {
		trieNode *root = decodeTrie();
		map64[bucket] = root;
		bucket = idx_reader->readBits64();
	}
	idx_reader->closeFile();
	fprintf(stderr, "Loaded index from file.\n");
	delete idx_reader;
}

void Hash::loadIdx64_d(std::string &fn) {
	idx_reader = new BitReader();
	idx_reader->openFile(fn);
	int doubly_unique = idx_reader->readBit();
	assert(doubly_unique == 1);
	int option = idx_reader->readBits(7);
	assert(option == 64);
	hash_len_ = idx_reader->readBits(8);
	uint64_t bucket = idx_reader->readBits64();
	while (bucket != END64) {
		trieNode *root = decodeTrie_d();
		map64[bucket] = root;
		bucket = idx_reader->readBits64();
	}
	idx_reader->closeFile();
	fprintf(stderr, "Loaded index from file.\n");
	delete idx_reader;
}

/*
	Count the number of unique substrings and trie nodes within each bucket.
 */
//std::vector<pleafNode*> traverse_p() {
//}


std::pair<uint32_t, uint32_t> Hash::countLeafNodes(const trieNode *root) {
	uint32_t n1 = 0, n2 = 0;
	if (root->isEnd) {
		n1 = ((leafNode*) root)->ucount;
		n2 = 1;
	}
	for (int i = 0; i < 4; i++) {
		if (root->children[i] != NULL) {
			std::pair<uint32_t, uint32_t> nn = countLeafNodes(root->children[i]);
			n1 += nn.first;
			n2 += nn.second;
		}
	}
	return std::make_pair(n1, n2);
}

std::pair<uint32_t, uint32_t> Hash::countdLeafNodes(const trieNode *root) {
	uint32_t n1 = 0, n2 = 0;
	if (root->isEnd) {
		std::pair<uint32_t, uint32_t> nn = ((dleafNode*) root)->ucount;
		n1 = nn.first + nn.second;
		n2 = 1;
	}
	for (int i = 0; i < 4; i++) {
		if (root->children[i] != NULL) {
			std::pair<uint32_t, uint32_t> nn = countLeafNodes(root->children[i]);
			n1 += nn.first;
			n2 += nn.second;
		}
	}
	return std::make_pair(n1, n2);
}

uint32_t Hash::countTrieNodes(const trieNode *root) {
	uint32_t r = 1;
	for (int i = 0; i < 4; i++)
		if (root->children[i] != NULL)
			r += countTrieNodes(root->children[i]);
	return r;
}

/*
	Encode index to a file (*.IDX)
 */
void Hash::encodeTrie(const trieNode *root) {
	if (root == NULL)
		idx_writer->writeBit(0);
	else {
		idx_writer->writeBit(1);
		for (int i = 0; i < 4; i++)
			encodeTrie(root->children[i]);
		if (root->isEnd) {
			idx_writer->writeBits32(((leafNode*) root)->refID);
			idx_writer->writeBits16(((leafNode*) root)->ucount);
			//idx_writer->writeBits32(32, ((leafNode*) root)->pos);
		}
	}
}

void Hash::encodeTrie_d(const trieNode *root) {
	if (root == NULL)
		idx_writer->writeBit(0);
	else {
		idx_writer->writeBit(1);
		for (int i = 0; i < 4; i++)
			encodeTrie_d(root->children[i]);
		if (root->isEnd) {
			idx_writer->writeBits32(((dleafNode*) root)->refID.first);
			idx_writer->writeBits32(((dleafNode*) root)->refID.second);
			idx_writer->writeBits16(((dleafNode*) root)->ucount.first);
			idx_writer->writeBits16(((dleafNode*) root)->ucount.second);
			//idx_writer->writeBits32(32, ((dleafNode*) root)->pos.first);
			//idx_writer->writeBits32(32, ((dleafNode*) root)->pos.second);
		}
	}
}

void Hash::encodeIdx32(std::string &ofn, int debug) {
	uint32_t num_buckets = 0;
	double avg_bucket_size = 0.0;
	uint32_t total_ = 0;
	
	if (debug) {
		for (auto it : map32) {
			std::pair<uint32_t, uint32_t> res_pair = countLeafNodes(it.second);
			//smap32[it.first] = res_pair;
			total_ += res_pair.first;
			num_buckets++;
			avg_bucket_size += res_pair.second;
		}
		fprintf(stderr, "Total number of unique substrings: %u.\n", total_);
		fprintf(stderr, "Total number of buckets: %u.\n", num_buckets);
		fprintf(stderr, "Average bucket size: %.4f.\n", (avg_bucket_size / num_buckets));
	}
	
	idx_writer = new BitWriter();
	idx_writer->openFile(ofn);
	
	/* Flag: doubly-unique = 0. */
	idx_writer->writeBit(0);

	/* Hash length. */
	idx_writer->writeBits(7, 32);
	idx_writer->writeBits(8, hash_len_);
	
	for (auto it : map32) {
		idx_writer->writeBits32(it.first);
		encodeTrie(it.second);
	}
	idx_writer->flush32();
	idx_writer->closeFile();
	fprintf(stderr, "Dumped index to file: %s.\n", ofn.c_str());

	delete idx_writer;
}

void Hash::encodeIdx32_d(std::string &ofn, int debug) {
	uint32_t num_buckets = 0;
	double avg_bucket_size = 0.0;
	uint32_t total_ = 0;

	if (debug) {
		for (auto it : map32) {
			std::pair<uint32_t, uint32_t> res_pair = countdLeafNodes(it.second);
			total_ += res_pair.first;
			num_buckets++;
			avg_bucket_size += res_pair.second;
		}
		fprintf(stderr, "Total number of unique substrings: %u.\n", total_);
		fprintf(stderr, "Total number of buckets: %u.\n", num_buckets);
		fprintf(stderr, "Average bucket size: %.4f.\n", (avg_bucket_size / num_buckets));
	}

	idx_writer = new BitWriter();
	idx_writer->openFile(ofn);
	
	/* Flag: doubly-unique = 1. */
	idx_writer->writeBit(1);

	/* Hash length. */
	idx_writer->writeBits(7, 32);
	idx_writer->writeBits(8, hash_len_);

	for (auto it : map32) {
		idx_writer->writeBits32(it.first);
		encodeTrie_d(it.second);
	}
	idx_writer->flush32();
	idx_writer->closeFile();
	fprintf(stderr, "Dumped index to file: %s.\n", ofn.c_str());

	delete idx_writer;
}

void Hash::encodeIdx64(std::string &ofn, int debug) {
	uint32_t num_buckets = 0;
	double avg_bucket_size = 0.0;
	uint32_t total_ = 0;
	
	if (debug) {
		for (auto it : map64) {
			std::pair<uint32_t, uint32_t> res_pair = countLeafNodes(it.second);
			//smap64[it.first] = res_pair;
			total_ += res_pair.first;
			num_buckets++;
			avg_bucket_size += res_pair.second;
		}
		fprintf(stderr, "Total number of unique substrings: %u.\n", total_);
		fprintf(stderr, "Total number of buckets: %u.\n", num_buckets);
		fprintf(stderr, "Average bucket size: %.4f.\n", (avg_bucket_size / num_buckets));
	}
	
	idx_writer = new BitWriter();
	idx_writer->openFile(ofn);

	/* Flag: doubly-unique = 0. */
	idx_writer->writeBit(0);

	/* Hash length. */
	idx_writer->writeBits(7, 64);
	idx_writer->writeBits(8, hash_len_);

	//int aaaaa = 0;
	for (auto it : map64) {
		//if (aaaaa < 10)
		//	fprintf(stderr, "%lu\n", it.first);
		//idx_writer->writeBits64(0);
		idx_writer->writeBits64(it.first);
		encodeTrie(it.second);
		//aaaaa++;
	}
	idx_writer->flush64();
	idx_writer->closeFile();
	fprintf(stderr, "Dumped index to file: %s.\n", ofn.c_str());

	delete idx_writer;
}

void Hash::encodeIdx64_d(std::string &ofn, int debug) {
	idx_writer = new BitWriter();
	idx_writer->openFile(ofn);

	idx_writer->writeBit(1);
	idx_writer->writeBits(7, 64);
	idx_writer->writeBits(8, hash_len_);
	
	for (auto it : map64) {
		idx_writer->writeBits64(it.first);
		encodeTrie_d(it.second);
	}
	idx_writer->flush64();
	idx_writer->closeFile();
	fprintf(stderr, "Dumped index to file: %s.\n", ofn.c_str());

	delete idx_writer;
}

//Updated 0331/2020
void Hash::encodeIdx64_temp(std::string &ofn_main, std::string &ofn_extra, int du) {
	idx_writer = new BitWriter();
	idx_writer->openFile(ofn_main);
	int indicator = 0;

	if (du)
		idx_writer->writeBit(1);
	else
		idx_writer->writeBit(0);
	idx_writer->writeBits(7, 64);
	//idx_writer->writeBits32(8, hash_len_);
	
	for (auto it : map64) {
		if (it.second->isEnd) {
			if (indicator < 10) {
				fprintf(stderr, "%lu\n", it.first);
				indicator++;
			}
			idx_writer->writeBits64(it.first);
			if (du) {
				idx_writer->writeBits32(((pleafNode*) it.second)->refID1);
				idx_writer->writeBits32(((pleafNode*) it.second)->refID2);
				idx_writer->writeBits16(((pleafNode*) it.second)->ucount1);
				idx_writer->writeBits16(((pleafNode*) it.second)->ucount2);
			} else {
				idx_writer->writeBits32(((pleafNode*) it.second)->refID1);
				idx_writer->writeBits16(((pleafNode*) it.second)->ucount1);
			}
		}
	}
	idx_writer->flush64();
	idx_writer->closeFile();
	fprintf(stderr, "Dumped main index to file: %s.\n", ofn_main.c_str());

	idx_writer->reset();
	idx_writer->openFile(ofn_extra);
	if (du)
		idx_writer->writeBit(1);
	else
		idx_writer->writeBit(0);
	idx_writer->writeBits(7, 64);
	idx_writer->writeBits(8, hash_len_);
	
	for (auto it : map64) {
		if (!it.second->isEnd) {
			idx_writer->writeBits64(it.first);
			encodeTrie_temp(it.second, du);
		}
	}
	idx_writer->flush64();
	idx_writer->closeFile();
	fprintf(stderr, "Dumped extra index to file: %s.\n", ofn_extra.c_str());

	delete idx_writer;
}

void Hash::encodeTrie_temp(const trieNode *root, int du) {
	if (root == NULL)
		idx_writer->writeBit(0);
	else {
		idx_writer->writeBit(1);
		for (int i = 0; i < 4; i++)
			encodeTrie_temp(root->children[i], du);
		if (root->isEnd) {
			if (du) {
				idx_writer->writeBits32(((pleafNode*) root)->refID1);
				idx_writer->writeBits32(((pleafNode*) root)->refID2);
				idx_writer->writeBits16(((pleafNode*) root)->ucount1);
				idx_writer->writeBits16(((pleafNode*) root)->ucount2);
			} else {
				idx_writer->writeBits32(((pleafNode*) root)->refID1);
				idx_writer->writeBits16(((pleafNode*) root)->ucount1);
			}
		}
	}
}

int Hash::symbolIdx[256] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
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

