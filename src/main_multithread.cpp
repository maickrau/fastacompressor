#include <thread>
#include <chrono>
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
	size_t numThreads = std::stoull(argv[3]);
	std::string filename { argv[4] };
	std::cerr << "k " << k << " w " << w << std::endl;
	auto programStartTime = getTime();
	FastaCompressor::CompressedStringIndex index { k, w };
	std::vector<std::string> realStrings;
	std::vector<std::tuple<size_t, std::string*, std::string*>> readStack;
	std::mutex stackMutex;
	size_t nextNum = 0;
	std::atomic<bool> readDone;
	readDone = false;
	std::vector<std::thread> threads;
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::cerr << "spawn threads" << std::endl;
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
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::cerr << "begin reading" << std::endl;
	std::vector<size_t> realReadLengths;
	FastQ::streamFastqFromFile(filename, false, [&readStack, &nextNum, &realReadLengths, &stackMutex](FastQ& read)
	{
		realReadLengths.emplace_back(read.sequence.size());
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
	readDone = true;
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::cerr << "read done" << std::endl;
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::cerr << "threads done, remove construction variables" << std::endl;
	index.removeConstructionVariables();
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	std::cerr << "construction done" << std::endl;
	index.printSizeInformation();
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
	for (size_t i = 0; i < realReadLengths.size(); i++)
	{
		assert(index.getSequence(i).size() == realReadLengths[i]);
	}
	std::cerr << "elapsed time " << formatTime(programStartTime, getTime()) << std::endl;
}
