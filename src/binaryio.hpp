#ifndef BINARYIO_HPP_
#define BINARYIO_HPP_

#include <fstream>
#include <vector>

#define BYTE 8
#define END32 0xFFFFFFFF
#define END64 0xFFFFFFFFFFFFFFFFul

class BitWriter {
	public:
		int test = 0;
		std::ofstream stream_INT;
		std::ofstream stream_AUX;
		int curByte = 0;
		int curBits = 0;
		BitWriter();
		~BitWriter() {
		}

		void openFile(const std::string &ofn);

		void writeBit(bool);
		void writeBits(int, uint32_t);
		void writeBits16(uint32_t);
		void writeBits32(uint32_t);
		void writeBits64(uint64_t);

		void flush32();
		void flush64();
		void flush32i();
		void flush32a();
		void flush64i();
		void flush64a();
		void closeFile();
		void reset() {
			curBits = 0;
			curByte = 0;
		}

};

class BitReader {
	public:
		std::ifstream stream_INT;
		std::ifstream stream_AUX;
		int curBits = 0;
		int curByte = 0;
		size_t cur_INT = 0, cur_AUX = 0;
		size_t fsize_INT = 0, fsize_AUX = 0;
		char* buffer_INT = NULL;
		char* buffer_AUX = NULL;

		/* Contructors */
		BitReader();
		~BitReader() {
			delete []buffer_INT;
			delete []buffer_AUX;
		}

		void openFile(const std::string &fn);

		uint32_t readBit();
		uint32_t readBits(int);
		uint16_t readBits16();
		uint32_t readBits32();
		uint64_t readBits64();

		void closeFile();
		void reset() {
			curBits = 0;
			curByte = 0;
		}
};

#endif

