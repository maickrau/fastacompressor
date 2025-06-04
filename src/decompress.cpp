#include <cassert>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <phmap.h>
#include "CompressedStringIndex.h"
#include "MinimizerIterator.h"
#include "fastqloader.h"

void writeFasta(const FastaCompressor::CompressedStringIndex& index, std::ostream& stream)
{
	for (size_t i = 0; i < index.size(); i++)
	{
		stream << ">" << index.getName(i) << std::endl;
		stream << index.getSequence(i) << std::endl;
	}
}

int main(int argc, char** argv)
{
	std::string inputFile;
	int opt;
	while ((opt = getopt(argc, argv, "i:h")) != -1)
	{
		switch(opt)
		{
		case 'i':
			if (inputFile.size() > 0)
			{
				std::cerr << "Only one input file allowed" << std::endl;
				std::abort();
			}
			inputFile = optarg;
			break;
		case 'h':
			std::cerr << "Parameters:" << std::endl;
			std::cerr << "-i arg Input file (default stdin)" << std::endl;
			std::cerr << "-h     Help" << std::endl;
			std::cerr << "Decompressed fasta without qualities is written to stdout" << std::endl;
			std::exit(0);
			break;
		}
	}
	if (inputFile.size() == 0)
	{
		inputFile = "-";
	}
	{
		std::istream* stream;
		if (inputFile == "-")
		{
			FastaCompressor::CompressedStringIndex index = FastaCompressor::CompressedStringIndex::loadFromStream(std::cin);
			writeFasta(index, std::cout);
		}
		else
		{
			std::ifstream file { inputFile, std::ios_base::binary };
			FastaCompressor::CompressedStringIndex index = FastaCompressor::CompressedStringIndex::loadFromStream(file);
			writeFasta(index, std::cout);
		}
	}
}
