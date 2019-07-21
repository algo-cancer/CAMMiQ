#ifndef BINARYIO_HPP_
#define BINARYIO_HPP_

#include <fstream>
#include <vector>

#define BYTE 8
#define END32 0xFFFFFFFF
#define END64 0xFFFFFFFFFFFFFFFF

class BitWriter {
	public:
		std::ofstream stream;
		int curByte = 0;
		int curBits = 0;
		BitWriter();
		~BitWriter() {
		}

		void openFile(const std::string &ofn);
		void writeBits32(int, uint32_t);
		void writeBits64(int, uint64_t);	
		void flush32();
		void flush64();
		void closeFile();
};

class BitReader {
	public:
		std::ifstream stream;
		int curBits = 0;
		int curByte = 0;
		//int test = 0;

		/* Contructors */
		BitReader();
		~BitReader() {
		}

		void openFile(const std::string &fn);
		uint32_t readBits32(int);
		uint64_t readBits64(int);
		void closeFile();
};

#endif

