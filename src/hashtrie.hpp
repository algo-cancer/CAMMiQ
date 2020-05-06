#ifndef HASHTRIE_HPP_
#define HASHTRIE_HPP_

#include <vector>
//#include <unordered_map>
#include "robin_hood.h"
#include "binaryio.hpp"

/*struct trieNode {
	std::vector<trieNode*> children;
	bool isEnd;
	trieNode();
	~trieNode();
};*/

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
	//int insertRefID(uint32_t);
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
	//int insertRefID(uint32_t);
	void reset();
};

typedef robin_hood::unordered_map<uint32_t, trieNode*> HashMap32;
typedef robin_hood::unordered_map<uint64_t, trieNode*> HashMap64;
//typedef std::unordered_map<uint32_t, pleafNode*> HashMapSP;
typedef robin_hood::unordered_map<uint32_t, std::vector<pleafNode*>> HashMapSP;

//typedef std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>> HashMapSum32;
//typedef std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>> HashMapSum64;
//typedef std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> Graph;

class Hash {
	private:
		uint32_t hash_len_ = 12;
		//uint8_t uniq_ = 1;
		HashMap32 map32;
		HashMap64 map64;
		//HashMapSum32 smap32;
		//HashMapSum64 smap64;
		
		BitWriter *idx_writer = NULL;
		BitReader *idx_reader = NULL;

	public:
		HashMapSP map_sp;
		//HashMapCnt map_cnt;
		uint64_t leaf_cnt = 0;
		//Graph us_graph;
		//std::unordered_map<uint32_t, float> expected_u1;
		//std::unordered_map<uint32_t, std::unordered_map<uint32_t, float>> expected_u2;
		//Hash(uint32_t, uint8_t);
		Hash(uint32_t);
		~Hash();

		/* Insert a unique substring to hash table. */
		trieNode* getNode();
		leafNode* getLeaf();
		dleafNode* getLeaf_d();
		pleafNode* getLeaf_p();

		/**/
		uint32_t computeHashVal(const uint8_t*);
		uint32_t computeHashVal(const char*);
		uint64_t computeHashVal64(const uint8_t*);
		uint64_t computeHashVal64(const char*);

		/**/
		//void insert32(uint32_t, const uint8_t*, uint32_t, uint32_t, uint8_t, uint32_t);
		//void insert64(uint64_t, const uint8_t*, uint32_t, uint32_t, uint8_t, uint32_t);
		void insert32(uint32_t, const uint8_t*, uint32_t, uint32_t, uint8_t);
		void insert64(uint64_t, const uint8_t*, uint32_t, uint32_t, uint8_t);
		leafNode* insert(trieNode*, const uint8_t*, uint32_t, uint32_t);
		leafNode* insert(trieNode*, const uint8_t*, uint32_t, uint32_t, uint8_t);
		void insert32_d(uint32_t, const uint8_t*, uint32_t, uint32_t, uint8_t, uint32_t, uint8_t);
		//void insert64_d(uint64_t, const uint8_t*, uint32_t, uint32_t, uint8_t, uint32_t, uint8_t, uint32_t);
		void insert64_d(uint64_t, const uint8_t*, uint32_t, uint32_t, uint8_t, uint32_t, uint8_t);
		dleafNode* insert_d(trieNode*, const uint8_t*, uint32_t);
		dleafNode* insert_d(trieNode*, const uint8_t*, uint32_t, uint32_t, uint8_t, uint32_t, uint8_t);
		void insert32_p(uint32_t, const uint8_t*, uint32_t);
		void insert64_p(uint64_t, const uint8_t*, uint32_t, uint32_t);
		void insert_p(trieNode*, const uint8_t*, uint32_t);

		/* Find a substring matching in the hash table. */
		uint32_t find32(uint32_t, const char*, size_t);
		uint32_t find32(uint32_t, const uint8_t*, size_t);
		uint32_t find32_test(uint32_t, const uint8_t*, size_t);
		std::pair<uint32_t, uint32_t> find32_d(uint32_t, const char*, size_t);
		std::pair<uint32_t, uint32_t> find32_d(uint32_t, const uint8_t*, size_t);
		pleafNode* find32_p(uint32_t, const uint8_t*, size_t);
		pleafNode* find64_p(uint64_t, const uint8_t*, size_t);
		uint32_t find64(uint64_t, const char*, size_t);
		std::vector<uint32_t> find32_all(uint32_t, const char*, size_t);
		std::vector<uint32_t> find64_all(uint64_t, const char*, size_t);
		int32_t find_and_increase(uint32_t, const uint8_t*, size_t, uint32_t, bool);
		int32_t find_and_increase_u(uint32_t, const uint8_t*, size_t);
		std::vector<uint32_t> find_rids(uint32_t, const uint8_t*, size_t);
		uint32_t find_rids_u(uint32_t, const uint8_t*, size_t);
		
		/* Dump the index into a binary file. */
		std::pair<uint32_t, uint32_t> countLeafNodes(const trieNode*);
		std::pair<uint32_t, uint32_t> countdLeafNodes(const trieNode*);
		uint32_t countTrieNodes(const trieNode*);
		//std::vector<pleafNode*> traverse_p();

		//trieNode* mergeTrie_p(trieNode*, trieNode*);

		void encodeTrie(const trieNode*);
		void encodeIdx32(std::string &ofn, int);
		void encodeIdx64(std::string &ofn, int);
		void encodeTrie_d(const trieNode*);
		void encodeIdx32_d(std::string &ofn, int);
		void encodeIdx64_d(std::string &ofn, int);

		/* Load the index from a binary file. */
		void loadIdx32(std::string &fn);
		void loadIdx32_d(std::string &fn);
		void loadIdx32_p(std::string &fn);
		void loadIdx64(std::string &fn);
		void loadIdx64_d(std::string &fn);
		trieNode* decodeTrie();
		trieNode* decodeTrie_d();
		trieNode* decodeTrie_p(int, uint8_t);
		void loadIdx64_p(std::string &fn);
		void loadIdx64_test(std::string &fn);
		//void loadIdx64_b(std::string&, std::string&);
		//trieNode* decodeTrie_g(int);
		//void loadIdx64_g(std::string &fn);
		void encodeIdx64_temp(std::string &, std::string &, int);
		void encodeTrie_temp(const trieNode*, int);

		/*void testbucket(uint32_t buc) {
			HashMap32::iterator iter = map32.find(buc);
			if (iter == map32.end()) {
				fprintf(stderr, "BBBB\n");
			}
			else 
				fprintf(stderr, "%d\n", iter->second);
		}*/

		HashMap32& getmap() {
			return map32;
		}

		static int symbolIdx[256];
};

#endif
