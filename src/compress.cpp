#include <thread>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <phmap.h>
#include "CompressedStringIndex.h"
#include "MinimizerIterator.h"
#include "fastqloader.h"

int main(int argc, char** argv)
{
	size_t numThreads = 1;
	size_t k = 5;
	size_t w = 100;
	std::vector<std::string> inputFiles;
	int opt;
	while ((opt = getopt(argc, argv, "j:i:k:w:h")) != -1)
	{
		switch(opt)
		{
		case 'j':
			numThreads = std::stoull(optarg);
			break;
		case 'i':
			inputFiles.emplace_back(optarg);
			break;
		case 'k':
			k = std::stoull(optarg);
			break;
		case 'w':
			w = std::stoull(optarg);
			break;
		case 'h':
			std::cerr << "Parameters:" << std::endl;
			std::cerr << "-j arg Number of compression threads (default 1)" << std::endl;
			std::cerr << "-i arg Input files (default stdin)" << std::endl;
			std::cerr << "-k arg K-mer size (default 5)" << std::endl;
			std::cerr << "-w arg Window size (default 100)" << std::endl;
			std::cerr << "-h     Help" << std::endl;
			std::cerr << "Compressed fasta without qualities is written to stdout" << std::endl;
			std::exit(0);
			break;
		}
	}
	if (inputFiles.size() == 0)
	{
		inputFiles.emplace_back("-");
	}
	FastaCompressor::CompressedStringIndex index { k, w };
	std::vector<std::tuple<size_t, std::string*, std::string*>> readStack;
	std::mutex stackMutex;
	std::atomic<bool> readDone;
	readDone = false;
	std::vector<std::thread> threads;
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&readDone, &index, &readStack, &stackMutex]
		{
			while (true)
			{
				std::tuple<size_t, std::string*, std::string*> read { std::numeric_limits<size_t>::max(), nullptr, nullptr };
				{
					std::lock_guard<std::mutex> lock { stackMutex };
					if (readStack.size() == 0)
					{
						if (readDone) break;
					}
					else
					{
						read = readStack.back();
						readStack.pop_back();
					}
				}
				if (std::get<0>(read) == std::numeric_limits<size_t>::max())
				{
					assert(std::get<1>(read) == nullptr);
					assert(std::get<2>(read) == nullptr);
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
					continue;
				}
				index.addString(std::get<0>(read), *std::get<1>(read), *std::get<2>(read));
				delete std::get<1>(read);
				delete std::get<2>(read);
			}
		});
	}
	size_t nextNum = 0;
	for (std::string filename : inputFiles)
	{
		if (filename == "-")
		{
			FastQ::streamFastqFastaFromStream(std::cin, false, [&readStack, &stackMutex, &nextNum](FastQ& read)
			{
				std::string* seq = new std::string;
				std::string* name = new std::string;
				std::swap(*seq, read.sequence);
				std::swap(*name, read.seq_id);
				while (true)
				{
					{
						std::lock_guard<std::mutex> lock { stackMutex };
						if (readStack.size() < 100)
						{
							readStack.emplace_back(nextNum, name, seq);
							break;
						}
					}
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
				}
				nextNum += 1;
			});
		}
		else
		{
			FastQ::streamFastqFromFile(filename, false, [&readStack, &stackMutex, &nextNum](FastQ& read)
			{
				std::string* seq = new std::string;
				std::string* name = new std::string;
				std::swap(*seq, read.sequence);
				std::swap(*name, read.seq_id);
				while (true)
				{
					{
						std::lock_guard<std::mutex> lock { stackMutex };
						if (readStack.size() < 100)
						{
							readStack.emplace_back(nextNum, name, seq);
							break;
						}
					}
					std::this_thread::sleep_for(std::chrono::milliseconds(10));
				}
				nextNum += 1;
			});
		}
	}
	readDone = true;
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
	index.removeConstructionVariables(numThreads);
	index.writeToStream(std::cout);
}
