#include "binaryio.hpp"

BitWriter::BitWriter() {
	curBits = 0;
	curByte = 0;
}

/* 
	Writes the 1/16/32/64 bits VALUE to file. 
 */
void BitWriter::writeBit(bool value) {
	int count = 1;
	/* Write next byte when available */
	try {
		if (curBits + count >= BYTE) {
			curByte = (curByte << 1) + value;
			count -= 1;
			stream_AUX.write((char*) &curByte, 1);
			curByte = 0;
			curBits = 0;
		}
	} catch (std::ofstream::failure e) {
		fprintf(stderr, "Exception writing to file.\n");
	}

	/* Add to current byte */
	curByte = (curByte << count) + ((count == 1) ? value : 0);
	curBits += count;
}

void BitWriter::writeBits(int count, uint32_t value) {
	bool value_;
	for (int i = 0; i < count; i++) {
		value_ = ((value >> (count - 1 - i)) & 0x1);
		writeBit(value_);
	}
}

void BitWriter::writeBits16(uint32_t value) {
	try {
		int value_ = (value >> 8) & 0xFF;
		stream_INT.write((char*) &value_, 1);
		value_ = value & 0xFF;
		stream_INT.write((char*) &value_, 1);
		/*uint32_t value_ = (value << 16);
		stream_INT.write((char*) &value_, 2);*/
	} catch (std::ofstream::failure e) {
		fprintf(stderr, "Exception writing to file.\n");
	}
}

void BitWriter::writeBits32(uint32_t value) {
	try {
		for (int i = 0; i < 4; i++) {
			int value_ = (value >> (8 * (3 - i))) & 0xFF;
			stream_INT.write((char*) &value_, 1);
		}
		//uint32_t value_ = value;
		//stream_INT.write((char*) &value_, 4);
	} catch (std::ofstream::failure e) {
		fprintf(stderr, "Exception writing to file.\n");
	}
}

void BitWriter::writeBits64(uint64_t value) {
	try {
		/*if (test == 0)
			fprintf(stderr, "--------%lu; %lu\n", value, sizeof(value));
		test++;
		uint64_t value_ = value;
		stream_INT.write((char*) &value_, 8);*/
		for (int i = 0; i < 8; i++) {
			int value_ = (value >> (8 * (7 - i))) & 0xFF;
			stream_INT.write((char*) &value_, 1);
		}
		//uint64_t value_ = value;
		//stream_INT.write((char*) &value_, 8);

	} catch (std::ofstream::failure e) {
		fprintf(stderr, "Exception writing to file.\n");
	}
}

void BitWriter::openFile(const std::string &ofn) {
	try {
		std::string aux_fn = ofn + ".aux";
		stream_INT.open(ofn.c_str(), std::ios::out | std::ios::binary);
		stream_AUX.open(aux_fn.c_str(), std::ios::out | std::ios::binary); 
	} catch (std::ofstream::failure e) { 
		fprintf(stderr, "Cannot open file: %s.\n", ofn.c_str()); 
	}
}

void BitWriter::flush32() {
	flush32i();
	flush32a();
}

void BitWriter::flush32i() {
	for (int i = 0; i < 40; i++)
		writeBit(1);
}

void BitWriter::flush32a() {
	writeBits32(END32);
	writeBits16(END32);
}

void BitWriter::flush64() {
	flush64i();
	flush64a();
}

void BitWriter::flush64i() {
	for (int i = 0; i < 72; i++)
		writeBit(1);
}

void BitWriter::flush64a() {
	writeBits64(END64);
	writeBits16(END32);
}

void BitWriter::closeFile() {
	/* Close file. */ 
	try {
		stream_INT.close();
		stream_AUX.close(); 
	} catch (std::ofstream::failure e) {
		fprintf(stderr, "Exception closing file.\n");
	} 
}

BitReader::BitReader() {
	curBits = 0;
	curByte = 0;
}

uint32_t BitReader::readBit() {
	uint32_t value = 0;
	
	if (curBits == 0) {
		curBits = BYTE;
		if (cur_AUX < fsize_AUX)
			curByte = buffer_AUX[cur_AUX++];
		else
			curByte = -1;
	}
	
	value = ((curByte >> (curBits - 1)) & 0x1u);
	curBits -= 1;
	return value;
}

uint32_t BitReader::readBits(int count) {
	uint32_t value = 0;
	for (int i = 0; i < count; i++)
		value = (value << 1) + readBit();
	return value;
}

uint16_t BitReader::readBits16() {
	//fprintf(stderr, "buffer_INT[8] = %d", buffer_INT[8]);
	uint16_t value = buffer_INT[cur_INT++] & 0xFF;
	value = ((value << 8) | (buffer_INT[cur_INT++] & 0xFF));
	return value;
}

uint32_t BitReader::readBits32() {
	uint32_t value = (buffer_INT[cur_INT++] & 0xFF);
	for (int i = 0; i < 3; i++)
		value = ((value << 8) | (buffer_INT[cur_INT++] & 0xFF));
	return value;
}

uint64_t BitReader::readBits64() {
	//fprintf(stderr, "value: %d, %d\n",  buffer_INT[0], buffer_INT[1]);
	uint64_t value = buffer_INT[cur_INT++] & 0xFF;
	for (int i = 0; i < 7; i++) {
		//fprintf(stderr, "value: %lu, %d\n", value, buffer_INT[cur_INT]);
		value = ((value << 8) | (buffer_INT[cur_INT++] & 0xFF));
	}
	return value;
}

void BitReader::openFile(const std::string &fn) {
	try {
		std::string aux_fn = fn + ".aux";
		//stream_INT.open(fn.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
		stream_INT.open(fn.c_str(), std::ios::in | std::ios::binary);
		stream_AUX.open(aux_fn.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
		
		stream_INT.seekg(0, std::ios::end);
		fsize_INT = stream_INT.tellg();
		//fprintf(stderr, "%lu\n", fsize_INT);
		buffer_INT = new char[fsize_INT + 100];
		//buffer_INT = new char[10000];
		//fprintf(stderr, "%d\n", buffer_INT[0]);
		stream_INT.seekg(0, std::ios::beg);
		stream_INT.read(buffer_INT, fsize_INT);
		//stream_INT.read(buffer_INT, 1000);
		//fprintf(stderr, "%d\n", buffer_INT[0]);

		stream_AUX.seekg(0, std::ios::end);
		fsize_AUX = stream_AUX.tellg();
		buffer_AUX = new char[fsize_AUX + 100];
		stream_AUX.seekg(0, std::ios::beg);
		stream_AUX.read(buffer_AUX, fsize_AUX);
	} catch (std::ifstream::failure e) {
		fprintf(stderr, "File does not exist.\n");
	}
}

void BitReader::closeFile() {
	try {
		stream_INT.close();
		stream_AUX.close();
	} catch (std::ifstream::failure e) {
		fprintf(stderr, "Exception closing file.\n");
	}
}

