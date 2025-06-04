#include <cassert>
#include <iostream>
#include <fstream>
#include <phmap.h>
#include "CompressedStringIndex.h"
#include "MinimizerIterator.h"
#include "fastqloader.h"

auto getTime()
{
	return std::chrono::steady_clock::now();
}

std::string formatTime(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end)
{
	size_t milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
	return std::to_string(milliseconds / 1000) + "," + std::to_string(milliseconds % 1000) + " s";
}

int main(int argc, char** argv)
{
	size_t k = std::stoull(argv[1]);
	size_t w = std::stoull(argv[2]);
	std::string filename { argv[3] };
	std::cerr << "k " << k << " w " << w << std::endl;
	auto programStartTime = getTime();
	FastaCompressor::CompressedStringIndex index { k, w };
	size_t countBases = 0;
	std::vector<std::string> realStrings;
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::cerr << "begin reading" << std::endl;
	FastQ::streamFastqFromFile(filename, false, [&index, &countBases, &realStrings](const FastQ& read)
	{
		countBases += read.sequence.size();
		index.addString(read.seq_id, read.sequence);
		realStrings.emplace_back(read.sequence);
	});
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::cerr << "read done, remove construction variables" << std::endl;
	index.removeConstructionVariables(1);
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	index.printSizeInformation();
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
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
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
}
