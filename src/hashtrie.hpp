#ifndef HASHTRIE_HPP_
#define HASHTRIE_HPP_

#include <vector>
#include <robin_hood.h>
#include "binaryio.hpp"

struct trieNode {
	trieNode** children = NULL;
	bool isEnd;
	trieNode();
	~trieNode();
};

/* Leaf Node: with one refID label. */
struct leafNode: public trieNode {
	uint32_t refID;
	uint32_t pos;
	uint16_t ucount;
	leafNode();
	leafNode(uint32_t);
};

/* Leaf Node: with two refID labels. */
struct dleafNode: public trieNode {
	std::pair<uint32_t, uint32_t> refID;
	std::pair<uint32_t, uint32_t> pos;
	std::pair<uint16_t, uint16_t> ucount;
	dleafNode();
	void setRefID_first(uint32_t, uint16_t);
	void setRefID_second(uint32_t, uint16_t);
	void resetRefID();
	int compareRefID(uint32_t, uint16_t, uint32_t, uint16_t);
};

/* Leaf Node: with one refID label and read count. */
struct pleafNode: public trieNode {
	uint32_t refID1 = 0;
	uint32_t refID2 = 0;
	uint8_t depth = 0;
	uint16_t ucount1 = 0;
	uint16_t ucount2 = 0;
	uint32_t rcount = 0;

	pleafNode();
	void reset();
};

typedef robin_hood::unordered_map<uint64_t, trieNode*> HashMap64;
typedef robin_hood::unordered_map<uint32_t, std::vector<pleafNode*>> HashMapSP;

class Hash {
	private:
		uint32_t hash_len_ = 26;
		HashMap64 map64;
		BitWriter *idx_writer = NULL;
		BitReader *idx_reader = NULL;

	public:
		HashMapSP map_sp;
		uint64_t leaf_cnt = 0;
		Hash(uint32_t);
		~Hash();
		void clear();

		/* Insert a unique substring to hash table. */
		trieNode* getNode();
		leafNode* getLeaf();
		dleafNode* getLeaf_d();
		pleafNode* getLeaf_p();

		/* Compute hash value for a given string. */
		uint64_t computeHashVal64(const uint8_t*);
		uint64_t computeHashVal64(const char*);

		/* Insert unique or doubly-unique substrings into hash table. */
		void insert64(uint64_t, const uint8_t*, uint32_t, uint32_t, uint8_t);
		leafNode* insert(trieNode*, const uint8_t*, uint32_t, uint32_t);
		leafNode* insert(trieNode*, const uint8_t*, uint32_t, uint32_t, uint8_t);
		void insert64_d(uint64_t, const uint8_t*, uint32_t, uint32_t, uint8_t, uint32_t, uint8_t);
		dleafNode* insert_d(trieNode*, const uint8_t*, uint32_t);
		dleafNode* insert_d(trieNode*, const uint8_t*, uint32_t, uint32_t, uint8_t, uint32_t, uint8_t);
		//void insert64_p(uint64_t, const uint8_t*, uint32_t, uint32_t);
		//void insert_p(trieNode*, const uint8_t*, uint32_t);

		/* Find a substring matching in the hash table. */
		pleafNode* find64_p(uint64_t, const uint8_t*, size_t);
		
		/* Dump the index into a binary file. */
		std::pair<uint32_t, uint32_t> countLeafNodes(const trieNode*);
		std::pair<uint32_t, uint32_t> countdLeafNodes(const trieNode*);
		uint32_t countTrieNodes(const trieNode*);
		void encodeTrie(const trieNode*);
		void encodeIdx64(std::string &ofn, int);
		void encodeTrie_d(const trieNode*);
		void encodeIdx64_d(std::string &ofn, int);

		/* Load the index from a binary file. */
		void loadIdx64(std::string &fn);
		void loadIdx64_d(std::string &fn);
		trieNode* decodeTrie();
		trieNode* decodeTrie_d();
		trieNode* decodeTrie_p(int, uint8_t);
		void loadIdx64_p(std::string &fn);

		HashMap64& getmap() {
			return map64;
		}
		uint32_t getHashLength() {
			return hash_len_;
		}

		static int symbolIdx[256];
};

#endif
