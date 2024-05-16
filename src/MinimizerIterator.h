#ifndef MinimizerIterator_h
#define MinimizerIterator_h

#include <cassert>
#include <vector>
#include <string>
#include <tuple>
#include <cstdint>
#include <queue>
#include <cmath>

namespace MinimizerIterator
{
	size_t charToInt(char c);
	std::vector<bool> getValidChars();
	uint64_t hash(uint64_t key);

	extern std::vector<bool> validChar;

	template <typename CallbackF>
	void iterateKmers(const std::string& str, size_t kmerLength, size_t windowSize, CallbackF callback)
	{
		const size_t realWindow = windowSize - kmerLength + 1;
		assert(kmerLength * 2 <= sizeof(size_t) * 8);
		if (str.size() < kmerLength) return;
		const size_t mask = ~(0xFFFFFFFFFFFFFFFF << (kmerLength * 2));
		assert(mask == pow(4, kmerLength)-1);
		size_t offset = 0;
	start:
		while (offset < str.size() && !validChar[str[offset]]) offset++;
		if (offset + kmerLength > str.size()) return;
		size_t kmer = 0;
		for (size_t i = 0; i < kmerLength; i++)
		{
			if (!validChar[str[offset+i]])
			{
				offset += i;
				goto start;
			}
			kmer <<= 2;
			kmer |= charToInt(str[offset+i]);
		}
		callback(offset + kmerLength-1, kmer);
		size_t lastKmer = kmer;
		size_t lastPos = offset + kmerLength-1;
		for (size_t i = kmerLength; offset+i < str.size(); i++)
		{
			if (!validChar[str[offset+i]])
			{
				offset += i;
				goto start;
			}
			kmer <<= 2;
			kmer &= mask;
			kmer |= charToInt(str[offset + i]);
			if (lastKmer != kmer || lastPos <= offset + i - realWindow)
			{
				callback(offset + i, kmer);
				lastKmer = kmer;
				lastPos = offset + i;
			}
		}
	}

	template <typename CallbackF>
	void iterateMinimizers(const std::string& str, size_t minimizerLength, size_t windowSize, CallbackF callback)
	{
		assert(minimizerLength * 2 <= sizeof(size_t) * 8);
		assert(minimizerLength <= windowSize);
		if (str.size() < minimizerLength) return;
		const size_t realWindow = windowSize - minimizerLength + 1;
		const size_t mask = ~(0xFFFFFFFFFFFFFFFF << (minimizerLength * 2));
		assert(mask == pow(4, minimizerLength)-1);
		size_t offset = 0;
		thread_local std::vector<std::tuple<size_t, size_t, size_t>> window;
	start:
		while (offset < str.size() && !validChar[str[offset]]) offset++;
		if (offset + windowSize > str.size()) return;
		size_t kmer = 0;
		for (size_t i = 0; i < minimizerLength; i++)
		{
			if (!validChar[str[offset+i]])
			{
				offset += i;
				goto start;
			}
			kmer <<= 2;
			kmer |= charToInt(str[offset+i]);
		}
		window.clear();
		window.emplace_back(offset+minimizerLength-1, kmer, hash(kmer));
		for (size_t i = minimizerLength; i < minimizerLength + realWindow; i++)
		{
			if (!validChar[str[offset + i]])
			{
				offset += i;
				goto start;
			}
			kmer <<= 2;
			kmer &= mask;
			kmer |= charToInt(str[offset + i]);
			auto hashed = hash(kmer);
			while (!window.empty() && std::get<2>(window.back()) > hashed) window.pop_back();
			window.emplace_back(offset+i, kmer, hashed);
		}
		auto iter = window.begin();
		while (iter != window.end() && std::get<2>(*iter) == std::get<2>(window.front()))
		{
			callback(std::get<0>(*iter), std::get<1>(*iter));
			++iter;
		}
		for (size_t i = minimizerLength + realWindow; offset+i < str.size(); i++)
		{
			if (!validChar[str[offset+i]])
			{
				offset += i;
				goto start;
			}
			kmer <<= 2;
			kmer &= mask;
			kmer |= charToInt(str[offset + i]);
			auto hashed = hash(kmer);
			size_t oldMinimum = std::get<2>(window[0]);
			bool frontPopped = false;
			while (!window.empty() && std::get<0>(window[0]) <= offset + i - realWindow)
			{
				frontPopped = true;
				window.erase(window.begin());
			}
			if (frontPopped)
			{
				while (window.size() >= 2 && std::get<2>(window[0]) == std::get<2>(*(window.begin()+1))) window.erase(window.begin());
			}
			while (!window.empty() && std::get<2>(window.back()) > hashed) window.pop_back();
			window.emplace_back(offset+i, kmer, hashed);
			if (std::get<2>(window[0]) != oldMinimum)
			{
				auto iter = window.begin();
				while (iter != window.end() && std::get<2>(*iter) == std::get<2>(window[0]))
				{
					callback(std::get<0>(*iter), std::get<1>(*iter));
					++iter;
				}
			}
			else if (std::get<2>(window.back()) == std::get<2>(window[0]))
			{
				callback(std::get<0>(window.back()), std::get<1>(window.back()));
			}
		}
	}
}

#endif
