#include <assert.h>
#include <chrono>

#include "hashtrie.hpp"
#include "util.hpp"
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
	//ucount = 0;
}

leafNode::leafNode(uint32_t refID_) {
	isEnd = true;
	refID = refID_;
	//ucount = 0;
}

dleafNode::dleafNode() {
	isEnd = true;
	refID = std::make_pair(0, 0);   
	ucount = std::make_pair(0, 0);   
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
	ucount1 = 0;
	ucount2 = 0;
	rcount = 0;
}

Hash::Hash(uint32_t hash_len) {
	hash_len_ = hash_len;
}

Hash::~Hash() {
	if (!map64.empty()) {
		for (auto it : map64)
			delete it.second;
		map64.clear();
	}
}

void Hash::clear() {
	/*if (!map32.empty()) {
		for (auto it : map32)
			delete it.second;
		map32.clear();
		fprintf(stderr, "HASHMAP: cleared 32-bits map.\n");
	}*/
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

uint64_t Hash::computeHashVal64(const uint8_t* key) {
	uint64_t res = 0;
	for (uint16_t i = 0; i < hash_len_; i++)
		res = ((res << 2) | symbolIdx[key[i]]);
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
			map64[bucket_] = leaf;
		} else {
			trieNode *root = getNode();
			map64[bucket_] = root;
 			leaf = insert(root, key + hash_len_, key_len_ - hash_len_, refID_);
		}
	}
	if (leaf != NULL) {
		leaf->refID = refID_;
		leaf->ucount = uc;
	}	
}

void Hash::insert64_d(uint64_t bucket_, const uint8_t* key, uint32_t key_len_, uint32_t refID_1,
			uint8_t uc1, uint32_t refID_2, uint8_t uc2) {
	HashMap64::iterator iter = map64.find(bucket_);
	dleafNode *dleaf = NULL;
	if (iter != map64.end()) {
		if (key_len_ == hash_len_) {
			if (iter->second->isEnd) {
				dleaf = (dleafNode*)iter->second;
				if (dleaf->compareRefID(refID_1, uc1, refID_2, uc2) == 0) {
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
			map64[bucket_] = dleaf;
		} else {
			trieNode *root = getNode();
			map64[bucket_] = root;
 			dleaf = insert_d(root, key + hash_len_, key_len_ - hash_len_);
		}
	}
	if (dleaf != NULL) {
		if (dleaf->refID.first == 0) {
			dleaf->setRefID_first(refID_1, uc1);
			dleaf->setRefID_second(refID_2, uc2);
		}
	}
}

/*
	The actual insert (helper) to a trie structure.
	KEY_LEN_ here (> 0) does not include hashed part
	Observe that none of the "key" can be a prefix of any other "key". 
 */
leafNode* Hash::insert(trieNode *root, const uint8_t* key, uint32_t key_len_, uint32_t refID_, uint8_t uc) {
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
	/* Mark last node as leaf and set its label. */
	index = symbolIdx[key[key_len_ - 1]];
	if (!cur->children[index]) {
		cur->children[index] = getLeaf();
		return (leafNode*) cur->children[index];
	} else {
		if (cur->children[index]->isEnd) {	
			leafNode *leaf = (leafNode*)cur->children[index];
			if (refID_ == leaf->refID) {
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
	/* Mark last node as leaf and set its label. */
	index = symbolIdx[key[key_len_ - 1]];
	if (!cur->children[index]) {
		cur->children[index] = getLeaf();
		return (leafNode*) cur->children[index];
	} else {
		if (cur->children[index]->isEnd) {	
			leafNode *leaf = (leafNode*)cur->children[index];
			if (refID_ == leaf->refID)
				return NULL;
		}
		fprintf(stderr, "INDEX %d\n", index);
		abortInsert();
	}
	return NULL;
}

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

trieNode* Hash::decodeTrie_p(int d_flag, uint8_t depth) {
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
				/* Doubly-unique substring must have two RIDs. */
				assert(rid1 != 0 && rid2 != 0);
				pleaf->refID1 = rid1;
				pleaf->refID2 = rid2;
				pleaf->ucount1 = idx_reader->readBits16();
				pleaf->ucount2 = idx_reader->readBits16();
				pleaf->rcount = 0;
				map_sp[rid1].push_back(pleaf);
				map_sp[rid2].push_back(pleaf);

				delete root;
				return pleaf;
			}
			return root;
		} else {
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
				map_sp[rid].push_back(pleaf);

				delete root;
				return pleaf;
			}
			return root;
		}
	}
}

void Hash::loadIdx64_p(std::string &fn) {
	idx_reader = new BitReader();
	idx_reader->openFile(fn);
	int doubly_unique = idx_reader->readBit();
	
	int option = idx_reader->readBits(7);
	assert(option == 64);
	hash_len_ = idx_reader->readBits(8);
	uint64_t bucket = idx_reader->readBits64();
	
	fprintf(stderr, "Index: %s\n", fn.c_str());
	fprintf(stderr, "Hash Length: %d\n", hash_len_);
	while (bucket != END64) {
		trieNode *root = decodeTrie_p(doubly_unique, 0);	
		map64[bucket] = root;
		bucket = idx_reader->readBits64();
	}
	idx_reader->closeFile();
	//fprintf(stderr, "Loaded index from file.\n");
	//fprintf(stderr, "Num buckets: %lu; %lu\n", count1, count2);
	delete idx_reader;
}


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
			std::pair<uint32_t, uint32_t> nn = countdLeafNodes(root->children[i]);
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
		}
	}
}

void Hash::encodeIdx64(std::string &ofn, int debug) {
	uint32_t num_buckets = 0;
	double avg_bucket_size = 0.0;
	uint32_t total_ = 0;
	
	if (debug) {
		for (auto it : map64) {
			std::pair<uint32_t, uint32_t> res_pair = countLeafNodes(it.second);
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

	for (auto it : map64) {
		idx_writer->writeBits64(it.first);
		encodeTrie(it.second);
	}
	idx_writer->flush64();
	idx_writer->closeFile();
	fprintf(stderr, "Dumped index to file: %s.\n", ofn.c_str());

	delete idx_writer;
}

void Hash::encodeIdx64_d(std::string &ofn, int debug) {
	uint32_t num_buckets = 0;
	double avg_bucket_size = 0.0;
	uint32_t total_ = 0;
	
	if (debug) {
		for (auto it : map64) {
			std::pair<uint32_t, uint32_t> res_pair = countdLeafNodes(it.second);
			total_ += res_pair.first;
			num_buckets++;
			avg_bucket_size += res_pair.second;
		}
		fprintf(stderr, "Total number of doubly unique substrings: %u.\n", total_);
		fprintf(stderr, "Total number of buckets: %u.\n", num_buckets);
		fprintf(stderr, "Average bucket size: %.4f.\n", (avg_bucket_size / num_buckets));
	}

	idx_writer = new BitWriter();
	idx_writer->openFile(ofn);

	/* Flag: doubly-unique = 1. */
	idx_writer->writeBit(1);

	/* Hash length. */
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

