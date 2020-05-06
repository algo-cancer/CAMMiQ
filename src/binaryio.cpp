/* 
	Extracted (and modified) from AssemblTire.
	see https://github.com/kyzhu/assembltrie/
 */
#include "binaryio.hpp"

BitWriter::BitWriter() {
	curBits = 0;
	curByte = 0;
}

uint32_t lowBits32(int count, uint32_t value) {
	return (uint32_t) ((value & END32) & ((1 << count) - 1));
}

int highBits32(int begin, int count, uint32_t value) {
	return (int) (value >> (begin - count));
}

uint64_t lowBits64(int count, uint64_t value) {
	return (uint64_t) ((value & END64) & ((1L << count) - 1));
}

int highBits64(int begin, int count, uint64_t value) {
	return (int) (value >> (begin - count));
}

/* 
	Writes the lowest COUNT bits of VALUE to file. 
 */
void BitWriter::writeBits32(int count, uint32_t value) {
	/* Write next byte when available */
	try {
		while (curBits + count >= BYTE) {
			int toWrite = BYTE - curBits;
			curByte = (curByte << toWrite) + highBits32(count, toWrite, value);
			count -= toWrite;
			value = lowBits32(count, value);

			stream.write((char*) &curByte, 1);
			curByte = 0;
			curBits = 0;
		}
	} catch (std::ofstream::failure e) {
		fprintf(stderr, "Exception writing to file.\n");
	}

	/* Add to current byte */
	curByte = (curByte << count) + lowBits32(count, value);
	curBits += count;
}

void BitWriter::writeBits64(int count, uint64_t value) {
	/* Write next byte when available */
	try {
		while (curBits + count >= BYTE) {
			int toWrite = BYTE - curBits;
			curByte = (curByte << toWrite) + highBits64(count, toWrite, value);
			count -= toWrite;
			value = lowBits64(count, value);

			stream.write((char*) &curByte, 1);
			curByte = 0;
			curBits = 0;
		}
	} catch (std::ofstream::failure e) {
		fprintf(stderr, "Exception writing to file.\n");
	}

	/* Add to current byte */
	curByte = (curByte << count) + lowBits64(count, value);
	curBits += count;
}

void BitWriter::openFile(const std::string &ofn) {
	try { 
		stream.open(ofn.c_str(), std::ios::out | std::ios::binary); 
	} catch (std::ofstream::failure e) { 
		fprintf(stderr, "Cannot open file: %s.\n", ofn.c_str()); 
	}
}

void BitWriter::flush32() {
	writeBits32(32, END32);
	writeBits32(8, 255);
}

void BitWriter::flush64() {
	writeBits64(64, END64);
	writeBits64(8, 255);
}

void BitWriter::closeFile() {
	/* Close file. */ 
	try { 
		stream.close(); 
	} catch (std::ofstream::failure e) {
		fprintf(stderr, "Exception closing file.\n");
	} 
}

BitReader::BitReader() {
	curBits = 0;
	curByte = 0;
}

uint32_t BitReader::readBits32(int count) {
	int curShift = 0;
	uint32_t value = 0;
	int need = count;
	try {
		while (need > curBits) {
			value = (value << curShift) + lowBits32(curBits, curByte);
			need -= curBits;
			curShift = BYTE;
			curBits = BYTE;
			if (!stream.eof())
			//if (cur < fsize) {
				stream.read((char*) &curByte, 1);
			/*	curByte = file[cur++];
			} else
				curByte = -1;*/
		}
	} catch (std::ifstream::failure e) {
		fprintf(stderr, "Exception reading file.\n");
	}
	value = (value << need) + (lowBits32(curBits, curByte) >> (curBits - need));
	curBits -= need;
	return value;
}

uint64_t BitReader::readBits64(int count) {
	int curShift = 0;
	uint64_t value = 0;
	int need = count;
	try {
		while (need > curBits) {
			value = (value << curShift) + lowBits64(curBits, curByte);
			need -= curBits;
			curShift = BYTE;
			curBits = BYTE;
			if (!stream.eof())
				stream.read((char*) &curByte, 1);
			/*if (cur < fsize) 
				curByte = file[cur++];
			else
				curByte = -1;*/
		}
	} catch (std::ifstream::failure e) {
		fprintf(stderr, "Exception reading file.\n");
	}
	value = (value << need) + (lowBits64(curBits, curByte) >> (curBits - need));
	curBits -= need;
	return value;
}

uint16_t BitReader::read_uint16_t() {
	uint16_t value;
	try {
		stream.read((char*) &value, 2);
	} catch (std::ifstream::failure e) {
		fprintf(stderr, "Exception reading file.\n");
	}
	return value;
}

uint32_t BitReader::read_uint32_t() {
	uint32_t value;
	try {
		stream.read((char*) &value, 4);
	} catch (std::ifstream::failure e) {
		fprintf(stderr, "Exception reading file.\n");
	}
	return value;
}

uint64_t BitReader::read_uint64_t() {
	uint64_t value = 0;
	/*try {
		stream.read(buffer, 8);
	} catch (std::ifstream::failure e) {
		fprintf(stderr, "Exception reading file.\n");
	}*/
	return value;
}

/*void init_buffer(int len) {
	buffer = new char[len];
}

void read_to_buffer(int len) {
	try {
		stream.read(buffer, len);
	} catch (std::ifstream::failure e) {
		fprintf(stderr, "Exception reading file.\n");
	}
}*/

void BitReader::openFile(const std::string &fn) {
	try {
		stream.open(fn.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
		
		stream.seekg(0, std::ios::beg);
		//char nextByte = 0;
		/*while(!stream.eof()) {
			stream.read((char*) &nextByte, 1);
			file.push_back(nextByte);
			fsize++;
		}*/
	} catch (std::ifstream::failure e) {
		fprintf(stderr, "File does not exist.\n");
	}
}

void BitReader::closeFile() {
	try {
		stream.close();
		//file.clear();
	} catch (std::ifstream::failure e) {
		fprintf(stderr, "Exception closing file.\n");
	}
}

