#include <string>
#include <cstring>
#include "binaryio.hpp"
#include "hashtrie.hpp"

/*void FqReader::readFastq(std::string &INFILE) { 
 	std::string bases;
	std::ifstream inputFile;

	/* Read in fasta file. / 
	inputFile.open(INFILE);
	while (std::getline(inputFile, bases)) {
		std::getline(inputFile, bases);
		read_count[(bases)]++;
		std::getline(inputFile, bases);
		std::getline(inputFile, bases);
	}
	
	inputFile.close();
}

void assign(std::string &INFILE)*/

void copyuint(uint8_t *dest, char * src, int begin, int len) {
	for (int i = 0; i < len; i++)
		dest[begin + i] = src[i];
}

int main(int argc, char** argv) {
	char* k1 = "ACAGT"; 
	char* k2 = "CACTAAT";
	char* k3 = "CACTACTCT";
	char* k4 = "GCCCCTC";
	char* k5 = "GCCCCTAC";

	uint8_t *a = NULL;

	a = new uint8_t[100];
	copyuint(a, k1, 0, 5);
	copyuint(a, k2, 5, 7);
	copyuint(a, k3, 12, 9);
	copyuint(a, k4, 21, 7);
	copyuint(a, k5, 28, 8);

	Hash *h = new Hash(5);

	uint32_t h1 = h->computeHashVal(k1);
	uint32_t h2 = h->computeHashVal(k2);
	uint32_t h3 = h->computeHashVal(k3);
	uint32_t h4 = h->computeHashVal(k4);
	uint32_t h5 = h->computeHashVal(k5);

	//fprintf(stderr, "%d, %d, %d, %d, %d.\n", h1, h2, h3, h4, h5);
	h->insert32(h1, a, 5, 1);
	h->insert32(h2, a + 5, 7, 2);
	h->insert32(h3, a + 12, 9, 3);
	h->insert32(h4, a + 21, 7, 4);
	h->insert32(h5, a + 28, 8, 5);

	std::string sss = "text.idx";
	h->encodeIdx32(sss);

	delete h;

	Hash *hh = new Hash(5);
	hh->loadIdx32(sss);
	delete hh;

	Hash *hhh = new Hash(5);

	uint64_t h1_ = hhh->computeHashVal64(k1);
	uint64_t h2_ = hhh->computeHashVal64(k2);
	uint64_t h3_ = hhh->computeHashVal64(k3);
	uint64_t h4_ = hhh->computeHashVal64(k4);
	uint64_t h5_ = hhh->computeHashVal64(k5);

	hhh->insert64(h1_, a, 5, 1);
	hhh->insert64(h2_, a + 5, 7, 2);
	hhh->insert64(h3_, a + 12, 9, 3);
	hhh->insert64(h4_, a + 21, 7, 4);
	hhh->insert64(h5_, a + 28, 8, 5);

	std::string ssss = "text2.idx";
	hhh->encodeIdx64(ssss);

	delete []a;
	delete hhh;

	Hash *hhhh = new Hash(5);
	hhhh->loadIdx64(ssss);
	delete hhhh;

	return 0;
}