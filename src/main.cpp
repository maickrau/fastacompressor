#include <cassert>
#include <iostream>
#include <fstream>
#include <phmap.h>
#include "CompressedStringIndex.h"
#include "MinimizerIterator.h"
#include "fastqloader.h"

int main(int argc, char** argv)
{
	size_t k = std::stoull(argv[1]);
	size_t w = std::stoull(argv[2]);
	std::string filename { argv[3] };
	std::cerr << "k " << k << " w " << w << std::endl;
	FastaCompressor::CompressedStringIndex index { k, w };
	size_t countBases = 0;
	std::vector<std::string> realStrings;
	FastQ::streamFastqFromFile(filename, false, [&index, &countBases, &realStrings](const FastQ& read)
	{
		countBases += read.sequence.size();
		index.addString(read.seq_id, read.sequence);
		realStrings.emplace_back(read.sequence);
	});
	index.removeConstructionVariables();
	index.printSizeInformation();
	for (size_t i = 0; i < realStrings.size(); i++)
	{
		if (!(realStrings[i] == index.getSequence(i)))
		{
			std::cerr << i << std::endl;
			std::cerr << realStrings[i] << std::endl;
			std::cerr << index.getSequence(i) << std::endl;
		}
		assert(realStrings[i] == index.getSequence(i));
		if (realStrings[i].size() > 20000)
		{
			assert(realStrings[i].substr(10000, 10000) == index.getSubstring(i, 10000, 10000));
		}
		if (realStrings[i].size() > 40000)
		{
			assert(realStrings[i].substr(30000, 10000) == index.getSubstring(i, 30000, 10000));
		}
	}
}
