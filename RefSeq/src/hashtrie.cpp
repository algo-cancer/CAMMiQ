#include "hashtrie.hpp"
#include "util.hpp"
#include "binaryio.hpp"

trieNode::trieNode() {
	isEnd = false;
}

trieNode::~trieNode() {
	for (int i = 0; i < 4; i++)
		if (children[i] != NULL)
			delete children[i];
}

leafNode::leafNode() {
	isEnd = true;
}

leafNode::leafNode(uint32_t refID_) {
	isEnd = true;
	refID = refID_;
}

dleafNode::dleafNode() {
	isEnd = true;
	refID = std::make_pair(0, 0);   
	ucount = std::make_pair(0, 0);   
}

void dleafNode::setRefID_first(uint32_t refID_) {
	refID.first = refID_;
	ucount.first = 1;
}

void dleafNode::setRefID_second(uint32_t refID_) {
	refID.second = refID_;
	ucount.second = 1;
}

int dleafNode::insertRefID(uint32_t refID_) {
	if (refID.first == refID_) {
		ucount.first++;
		return 0;
	}
	if (refID.second == 0) {
		refID.second = refID_;
		ucount.second = 1;
		return 1;
	}
	if (refID.second == refID_) {
		ucount.second++;
		return 2;
	}
	fprintf(stderr, "%u\t%u\n", refID.first, refID_);
	return -1;
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
	if (!smap32.empty())
		smap32.clear();
	if (!smap64.empty())
		smap64.clear();
}

trieNode* Hash::getNode() {
	trieNode *pNode = new trieNode();
	for (int i = 0; i < 4; i++)
        	pNode->children.push_back(NULL);
	return pNode;
}

leafNode* Hash::getLeaf() {
	leafNode *lNode = new leafNode();
	lNode->ucount = 1;
	for (int i = 0; i < 4; i++)
        	lNode->children.push_back(NULL);
	return lNode;
}

dleafNode* Hash::getLeaf_d() {
	dleafNode *dNode = new dleafNode();
	for (int i = 0; i < 4; i++)
        	dNode->children.push_back(NULL);
	return dNode;
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
void Hash::insert32(uint32_t bucket_, const uint8_t* key, uint32_t key_len_, uint32_t refID_) {
	HashMap32::iterator iter = map32.find(bucket_);
	leafNode *leaf = NULL;
	if (iter != map32.end()) {
		if (key_len_ == hash_len_) {
			if (iter->second->isEnd) {
				leaf = (leafNode*)iter->second;
				if (refID_ == leaf->refID) {
					leaf->ucount++;
					return;
				}
			}
			abortInsert();
		}
		leaf = insert(iter->second, key + hash_len_, key_len_ - hash_len_, refID_);
	} else {
		if (key_len_ == hash_len_) {
			leaf = getLeaf();
			map32.insert(std::make_pair(bucket_, leaf));
		} else {
			trieNode *root = getNode();
			map32.insert(std::make_pair(bucket_, root));
 			leaf = insert(root, key + hash_len_, key_len_ - hash_len_, refID_);
		}
	}
	if (leaf != NULL)
		leaf->refID = refID_;
}

void Hash::insert32_d(uint32_t bucket_, const uint8_t* key, uint32_t key_len_, uint32_t refID_) {
	HashMap32::iterator iter = map32.find(bucket_);
	dleafNode *dleaf = NULL;
	if (iter != map32.end()) {
		if (key_len_ == hash_len_) {
			if (iter->second->isEnd) {
				dleaf = (dleafNode*)iter->second;
				if (dleaf->insertRefID(refID_) < 0) {
					abortInsert();
				}
				return;
			} else
				abortInsert();
		}
		insert_d(iter->second, key + hash_len_, key_len_ - hash_len_, refID_);
	} else {
		if (key_len_ == hash_len_) {
			dleaf = getLeaf_d();
			dleaf->setRefID_first(refID_);
			map32.insert(std::make_pair(bucket_, dleaf));
		} else {
			trieNode *root = getNode();
			map32.insert(std::make_pair(bucket_, root));
 			insert_d(root, key + hash_len_, key_len_ - hash_len_, refID_);
		}
	}
}

void Hash::insert64(uint64_t bucket_, const uint8_t* key, uint32_t key_len_, uint32_t refID_) {
	//uint64_t bucket_ = computeHashVal64(key);
	HashMap64::iterator iter = map64.find(bucket_);
	leafNode *leaf = NULL;
	if (iter != map64.end()) {
		if (key_len_ == hash_len_) {
			if (iter->second->isEnd) {
				leaf = (leafNode*)iter->second;
				if (refID_ == leaf->refID) {
					leaf->ucount++;
					return;
				}
			}
			abortInsert();
		}
		leaf = insert(iter->second, key + hash_len_, key_len_ - hash_len_, refID_);
	} else {
		if (key_len_ == hash_len_) {
			leaf = getLeaf();
			map64.insert(std::make_pair(bucket_, leaf));
		} else {
			trieNode *root = getNode();
			map64.insert(std::make_pair(bucket_, root));
 			leaf = insert(root, key + hash_len_, key_len_ - hash_len_, refID_);
		}
	}
	if (leaf != NULL)
		leaf->refID = refID_;	
}

void Hash::insert64_d(uint64_t bucket_, const uint8_t* key, uint32_t key_len_, uint32_t refID_) {
	HashMap64::iterator iter = map64.find(bucket_);
	dleafNode *dleaf = NULL;
	if (iter != map64.end()) {
		if (key_len_ == hash_len_) {
			if (iter->second->isEnd) {
				dleaf = (dleafNode*)iter->second;
				if (dleaf->insertRefID(refID_) < 0)
					abortInsert();
			}
			abortInsert();
		}
		insert_d(iter->second, key + hash_len_, key_len_ - hash_len_, refID_);
	} else {
		if (key_len_ == hash_len_) {
			dleaf = getLeaf_d();
			dleaf->setRefID_first(refID_);
			map64.insert(std::make_pair(bucket_, dleaf));
		} else {
			trieNode *root = getNode();
			map64.insert(std::make_pair(bucket_, root));
 			insert_d(root, key + hash_len_, key_len_ - hash_len_, refID_);
		}
	}
}

/*
	The actual insert (helper) to a trie structure.
	KEY_LEN_ here (> 0) does not include hashed part
	Observe that none of the "key" can be a prefix of any other "key". 
 */
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
			if (refID_ == leaf->refID) {
				leaf->ucount++;
				return NULL;
			}
		}
		abortInsert();
	}
	return NULL;
}

void Hash::insert_d(trieNode *root, const uint8_t* key, uint32_t key_len_, uint32_t refID_) {
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
	dleafNode *dleaf;
	if (!cur->children[index]) {
		cur->children[index] = getLeaf_d();
		dleaf = ((dleafNode*) cur->children[index]);
		dleaf->setRefID_first(refID_);
	} else {
		if (cur->children[index]->isEnd) {	
			dleaf = ((dleafNode*) cur->children[index]);
			if (dleaf->insertRefID(refID_) < 0)
				abortInsert();
		} else 
			abortInsert();
	}
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
	int b = idx_reader->readBits32(1);
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
			leaf->refID = idx_reader->readBits32(32);
			delete root;
			return leaf;
		}
		return root;
	}
}

void Hash::loadIdx32(std::string &fn) {
	idx_reader = new BitReader();
	idx_reader->openFile(fn);
	uint32_t bucket = idx_reader->readBits32(32);
	while (bucket != END32) {
		trieNode *root = decodeTrie();
		map32[bucket] = root;
		bucket = idx_reader->readBits32(32);
	}
	idx_reader->closeFile();
	fprintf(stderr, "Loaded index from file.\n");
	delete idx_reader;
}

void Hash::loadIdx64(std::string &fn) {
	idx_reader = new BitReader();
	idx_reader->openFile(fn);
	uint64_t bucket = idx_reader->readBits64(64);
	while (bucket != END64) {
		trieNode *root = decodeTrie();
		map64[bucket] = root;
		bucket = idx_reader->readBits64(64);
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
		idx_writer->writeBits32(1, 0);
	else {
		idx_writer->writeBits32(1, 1);
		for (int i = 0; i < 4; i++)
			encodeTrie(root->children[i]);
		if (root->isEnd)
			idx_writer->writeBits32(32, ((leafNode*) root)->refID);
	}
}

// Just for test 2019/07/18: need to update for doubly-unique later.
This needs to be modified since we need to print also the counts. 
void Hash::encodeTrie_d(const trieNode *root, int depth) {
	if (root != NULL) {
		for (int i = 0; i < 4; i++)
			encodeTrie_d(root->children[i], depth + 1);
		if (root->isEnd) {
			uint32_t r1 = ((dleafNode*) root)->refID.first;
			uint32_t r2 = ((dleafNode*) root)->refID.second;
			fprintf(stderr, "%d\t%u\t%u\n", depth, r1, r2);
		}
	}
}

void Hash::encodeIdx32_d() {
	for (auto it : map32) {
		encodeTrie_d(it.second, 11);
	}
}

void Hash::encodeIdx32(std::string &ofn) {
	uint32_t num_buckets = 0;
	double avg_bucket_size = 0.0;
	uint32_t total_ = 0;
	
	//fprintf(stderr, "\nNum buckets: %lu\n", map32.size());
	for (auto it : map32) {
		std::pair<uint32_t, uint32_t> res_pair = countLeafNodes(it.second);
		smap32[it.first] = res_pair;
		total_ += res_pair.first;
		num_buckets++;
		avg_bucket_size += res_pair.second;
	}
	fprintf(stderr, "Total number of unique substrings: %u.\n", total_);
	fprintf(stderr, "Total number of buckets: %u.\n", num_buckets);
	fprintf(stderr, "Average bucket size: %.4f.\n", (avg_bucket_size / num_buckets));
	
	idx_writer = new BitWriter();
	idx_writer->openFile(ofn);
	for (auto it : map32) {
		idx_writer->writeBits32(32, it.first);
		encodeTrie(it.second);
	}
	idx_writer->flush32();
	idx_writer->closeFile();
	fprintf(stderr, "Dumped index to file: %s.\n", ofn.c_str());

	delete idx_writer;
}

void Hash::encodeIdx64(std::string &ofn) {
	uint64_t num_buckets = 0;
	double avg_bucket_size = 0.0;
	//uint32_t test = 0;
	for (auto it : map64) {
		std::pair<uint32_t, uint32_t> res_pair = countLeafNodes(it.second);
		smap64[it.first] = res_pair;
	//	test += countTrieNodesC(it.second);
		num_buckets++;
		avg_bucket_size += res_pair.second;
	}
	//fprintf(stderr, "CCCC: %u.\n", test);
	fprintf(stderr, "Total number of buckets: %lu.\n", num_buckets);
	fprintf(stderr, "Average bucket size: %.4f.\n", (avg_bucket_size / num_buckets));
	
	idx_writer = new BitWriter();
	idx_writer->openFile(ofn);
	for (auto it : map64) {
		idx_writer->writeBits64(64, it.first);
		encodeTrie(it.second);
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

