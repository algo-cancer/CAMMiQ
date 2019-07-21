#ifndef HASHTRIE_HPP_
#define HASHTRIE_HPP_

#include <vector>
#include <unordered_map>
#include "binaryio.hpp"

struct trieNode {
	std::vector<trieNode*> children;
	bool isEnd;
	trieNode();
	~trieNode();
};

struct leafNode: public trieNode {
	uint32_t refID;
	uint32_t ucount;
	leafNode();
	leafNode(uint32_t);
};

struct dleafNode: public trieNode {
	std::pair<uint32_t, uint32_t> refID;
	std::pair<uint32_t, uint32_t> ucount;
	dleafNode();
	void setRefID_first(uint32_t);
	void setRefID_second(uint32_t);
	int insertRefID(uint32_t);
};

typedef std::unordered_map<uint32_t, trieNode*> HashMap32;
typedef std::unordered_map<uint64_t, trieNode*> HashMap64;
typedef std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>> HashMapSum32;
typedef std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>> HashMapSum64;

class Hash {
	private:
		uint32_t hash_len_ = 12;
		//uint8_t uniq_ = 1;
		HashMap32 map32;
		HashMap64 map64;
		HashMapSum32 smap32;
		HashMapSum64 smap64;
		BitWriter *idx_writer = NULL;
		BitReader *idx_reader = NULL;

	public:
		//Hash(uint32_t, uint8_t);
		Hash(uint32_t);
		~Hash();

		/* Insert a unique substring to hash table. */
		trieNode* getNode();
		leafNode* getLeaf();
		dleafNode* getLeaf_d();
		uint32_t computeHashVal(const uint8_t*);
		uint32_t computeHashVal(const char*);
		uint64_t computeHashVal64(const uint8_t*);
		uint64_t computeHashVal64(const char*);
		void insert32(uint32_t, const uint8_t*, uint32_t, uint32_t);
		void insert64(uint64_t, const uint8_t*, uint32_t, uint32_t);
		leafNode* insert(trieNode*, const uint8_t*, uint32_t, uint32_t);
		void insert32_d(uint32_t, const uint8_t*, uint32_t, uint32_t);
		void insert64_d(uint64_t, const uint8_t*, uint32_t, uint32_t);
		void insert_d(trieNode*, const uint8_t*, uint32_t, uint32_t);

		/* Find a substring matching in the hash table. */
		uint32_t find32(uint32_t, const char*, size_t);
		uint32_t find64(uint64_t, const char*, size_t);
		std::vector<uint32_t> find32_all(uint32_t, const char*, size_t);
		std::vector<uint32_t> find64_all(uint64_t, const char*, size_t);
		
		/* Dump the index into a binary file. */
		std::pair<uint32_t, uint32_t> countLeafNodes(const trieNode*);
		uint32_t countTrieNodes(const trieNode*);

		void encodeTrie(const trieNode*);
		void encodeIdx32(std::string &ofn);
		void encodeIdx64(std::string &ofn);
		void encodeTrie_d(const trieNode*, int);
		void encodeIdx32_d();

		/* Load the index from a binary file. */
		void loadIdx32(std::string &fn);
		void loadIdx64(std::string &fn);
		trieNode* decodeTrie();

		static int symbolIdx[256];
};

#endif