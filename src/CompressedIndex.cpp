#include <iostream>
#include <chrono>
#include <thread>
#include "CompressedIndex.h"
#include "MinimizerIterator.h"
#include "RankBitvector.h"

namespace FastaCompressor
{
	CompressedIndex::CompressedIndex(size_t k, size_t w) :
		hierarchyIndex(),
		pieceIndex(),
		hierarchyTopDownFirst(),
		hierarchyTopDownSecond(),
		pieces(),
		firstHierarchicalIndex(std::numeric_limits<size_t>::max()),
		isfrozen(false),
		bitsPerIndex(0),
		k(k),
		w(w)
	{
		pieceReaderCount[0] = 0;
		pieceReaderCount[1] = 0;
		pieceReaderCount[2] = 0;
		pieceReaderCount[3] = 0;
		hierarchyReaderCount = 0;
	}
	std::vector<size_t> CompressedIndex::segmentFastaToPieces(const std::string& seq)
	{
		std::vector<size_t> minimizerPositions;
		MinimizerIterator::iterateMinimizers(seq, k, w+k-1, [&minimizerPositions](size_t pos, size_t kmer)
		{
			minimizerPositions.emplace_back(pos);
		});
		if (minimizerPositions.size() == 0)
		{
			minimizerPositions.emplace_back(0);
			minimizerPositions.emplace_back(seq.size());
		}
		else
		{
			if (minimizerPositions[0] != 0) minimizerPositions.insert(minimizerPositions.begin(), 0);
			minimizerPositions.emplace_back(seq.size());
		}
		std::vector<__uint128_t> shortPieces;
		std::vector<std::string> longPieces;
		shortPieces.resize(minimizerPositions.size()-1);
		longPieces.resize(minimizerPositions.size()-1);
		std::vector<std::vector<size_t>> piecesPerType;
		piecesPerType.resize(4);
		for (size_t i = 1; i < minimizerPositions.size(); i++)
		{
			size_t start = minimizerPositions[i-1];
			size_t end = minimizerPositions[i];
			std::string kmer = seq.substr(start, end-start);
			if (kmer.size() <= 15)
			{
				piecesPerType[0].emplace_back(i-1);
				shortPieces[i-1] = pieceIndex.encodeString(kmer);
			}
			else if (kmer.size() <= 31)
			{
				piecesPerType[1].emplace_back(i-1);
				shortPieces[i-1] = pieceIndex.encodeString(kmer);
			}
			else if (kmer.size() <= 63)
			{
				piecesPerType[2].emplace_back(i-1);
				shortPieces[i-1] = pieceIndex.encodeString(kmer);
			}
			else
			{
				piecesPerType[3].emplace_back(i-1);
				longPieces[i-1] = pieceIndex.encodeStringToString(kmer);
			}
		}
		std::vector<size_t> result;
		result.resize(minimizerPositions.size()-1, std::numeric_limits<size_t>::max());
		for (size_t type = 0; type < 4; type++)
		{
			if (piecesPerType[type].size() == 0) continue;
			std::vector<size_t> indicesNeedWriting;
			{
				std::lock_guard lock { pieceMutex[type] };
				pieceReaderCount[type] += 1;
			}
			for (size_t i : piecesPerType[type])
			{
				size_t index = 0;
				index = pieceIndex.getIndexOrNull(shortPieces[i], longPieces[i], type);
				result[i] = index;
				if (index == std::numeric_limits<size_t>::max())
				{
					indicesNeedWriting.emplace_back(i);
				}
			}
			pieceReaderCount[type] -= 1;
			if (pieceReaderCount[type] == 0)
			{
				pieceConditionVariable[type].notify_one();
			}
			if (indicesNeedWriting.size() >= 1)
			{
				{
					std::unique_lock<std::mutex> lock { pieceMutex[type] };
					while (pieceReaderCount[type] > 0)
					{
						pieceConditionVariable[type].wait_for(lock, std::chrono::milliseconds(1));
					}
					assert(pieceReaderCount[type] == 0);
					for (size_t ii = 0; ii < indicesNeedWriting.size(); ii++)
					{
						size_t i = indicesNeedWriting[ii];
						size_t index = 0;
						size_t addedIndex = pieceIndex.typeCount(type)*8 + (type*2);
						index = pieceIndex.getIndexOrSet(shortPieces[i], longPieces[i], type, addedIndex);
						result[i] = index;
					}
				}
				pieceConditionVariable[type].notify_one();
			}
		}
		for (size_t i = 0; i < result.size(); i++)
		{
			assert(result[i] != std::numeric_limits<size_t>::max());
		}
		return result;
	}
	std::vector<size_t> CompressedIndex::addString(const std::string& str)
	{
		auto indices = segmentFastaToPieces(str);
		{
			std::lock_guard<std::mutex> lock { hierarchyMutex };
			hierarchyReaderCount += 1;
		}
		auto fixedIndices = hierarchizeIndices(indices);
		hierarchyReaderCount -= 1;
		if (hierarchyReaderCount == 0)
		{
			hierarchyConditionVariable.notify_one();
		}
		std::vector<size_t> result;
		result.reserve(fixedIndices.size());
		assert(fixedIndices.size() >= 1);
		if (fixedIndices.size() >= 2)
		{
			std::unique_lock<std::mutex> lock { hierarchyMutex };
			while (hierarchyReaderCount > 0)
			{
				hierarchyConditionVariable.wait_for(lock, std::chrono::milliseconds(1));
			}
			assert(hierarchyReaderCount == 0);
			for (size_t i = 0; i+1 < fixedIndices.size(); i += 2)
			{
				if (hierarchyIndex.count(std::make_pair(fixedIndices[i], fixedIndices[i+1])) == 1)
				{
					result.emplace_back(hierarchyIndex.at(std::make_pair(fixedIndices[i], fixedIndices[i+1])));
				}
				else
				{
					assert(hierarchyIndex.count(std::make_pair(fixedIndices[i], fixedIndices[i+1])) == 0);
					size_t index = hierarchyIndex.size()*2+1;
					hierarchyIndex.set(std::make_pair(fixedIndices[i], fixedIndices[i+1]), index);
					result.emplace_back(index);
				}
			}
		}
		hierarchyConditionVariable.notify_one();
		if (fixedIndices.size() % 2 == 1) result.emplace_back(fixedIndices.back());
		return result;
	}
	VariableWidthIntVector CompressedIndex::hierarchizeIndices(const VariableWidthIntVector& indices) const
	{
		std::vector<size_t> tmp;
		for (size_t i = 0; i < indices.size(); i++)
		{
			tmp.emplace_back(i);
		}
		auto tmp2 = hierarchizeIndices(tmp);
		size_t maxHere = 0;
		VariableWidthIntVector result;
		for (size_t val : tmp2)
		{
			maxHere = std::max(maxHere, val);
		}
		result.setWidth(ceil(log2(maxHere+1)));
		result.resize(tmp2.size());
		for (size_t i = 0; i < tmp2.size(); i++)
		{
			result.set(i, tmp2[i]);
		}
		return result;
	}
	std::vector<size_t> CompressedIndex::hierarchizeIndices(const std::vector<size_t>& indices) const
	{
		std::vector<size_t> result;
		for (size_t i = 0; i < indices.size(); i++)
		{
			result.emplace_back(indices[i]);
			while (result.size() >= 2 && hierarchyIndex.count(std::make_pair(result[result.size()-2], result[result.size()-1])) == 1)
			{
				size_t replacement = hierarchyIndex.at(std::make_pair(result[result.size()-2], result[result.size()-1]));
				// assert(replacement < hierarchyIndex.size() + pieceIndex.size()); // true but can't check because no lock on pieceIndex
				result.pop_back();
				result.pop_back();
				result.emplace_back(replacement);
			}
		}
		return result;
	}
	std::string CompressedIndex::getSubstring(const VariableWidthIntVector& indices, const size_t startPos, const size_t length) const
	{
		std::vector<size_t> fixed;
		fixed.reserve(indices.size());
		for (size_t i = 0; i < indices.size(); i++)
		{
			fixed.emplace_back(indices.get(i));
		}
		return getSubstring(fixed, startPos, length);
	}
	std::string CompressedIndex::getString(const VariableWidthIntVector& indices) const
	{
		std::vector<size_t> fixed;
		fixed.reserve(indices.size());
		for (size_t i = 0; i < indices.size(); i++)
		{
			fixed.emplace_back(indices.get(i));
		}
		return getString(fixed);
	}
	std::string CompressedIndex::getSubstring(const std::vector<size_t>& indices, const size_t startPos, const size_t length) const
	{
		assert(frozen());
		std::string result;
		std::vector<size_t> baseIndices;
		std::vector<size_t> hierarchyIndices = indices;
		while (hierarchyIndices.size() > 0)
		{
			while (hierarchyIndices.back() >= firstHierarchicalIndex)
			{
				size_t index = hierarchyIndices.back() - firstHierarchicalIndex;
				assert(index < hierarchyTopDownFirst.size());
				hierarchyIndices.pop_back();
				hierarchyIndices.emplace_back(hierarchyTopDownFirst.get(index));
				hierarchyIndices.emplace_back(hierarchyTopDownSecond.get(index));
			}
			assert(hierarchyIndices.back() < firstHierarchicalIndex);
			baseIndices.push_back(hierarchyIndices.back());
			hierarchyIndices.pop_back();
		}
		size_t currentPos = 0;
		for (size_t i = baseIndices.size()-1; i < baseIndices.size(); i--)
		{
			assert(baseIndices[i] < pieces.size());
			size_t pieceSize = pieces.pieceSize(baseIndices[i]);
			if (currentPos + pieceSize <= startPos)
			{
				currentPos += pieceSize;
				continue;
			}
			result += pieces.get(baseIndices[i]);
			if (currentPos < startPos)
			{
				assert(result.size() > startPos-currentPos);
				result.erase(result.begin(), result.begin()+(startPos-currentPos));
			}
			currentPos += pieceSize;
			if (currentPos >= startPos+length)
			{
				if (currentPos > startPos+length)
				{
					assert(result.size() > currentPos-(startPos+length));
					result.erase(result.end()-(currentPos-(startPos+length)), result.end());
				}
				break;
			}
		}
		return result;
	}
	std::string CompressedIndex::getString(const std::vector<size_t>& indices) const
	{
		assert(frozen());
		std::string result;
		std::vector<size_t> baseIndices;
		std::vector<size_t> hierarchyIndices = indices;
		while (hierarchyIndices.size() > 0)
		{
			while (hierarchyIndices.back() >= firstHierarchicalIndex)
			{
				size_t index = hierarchyIndices.back() - firstHierarchicalIndex;
				assert(index < hierarchyTopDownFirst.size());
				hierarchyIndices.pop_back();
				hierarchyIndices.emplace_back(hierarchyTopDownFirst.get(index));
				hierarchyIndices.emplace_back(hierarchyTopDownSecond.get(index));
			}
			assert(hierarchyIndices.back() < firstHierarchicalIndex);
			baseIndices.push_back(hierarchyIndices.back());
			hierarchyIndices.pop_back();
		}
		for (size_t i = baseIndices.size()-1; i < baseIndices.size(); i--)
		{
			assert(baseIndices[i] < pieces.size());
			result += pieces.get(baseIndices[i]);
		}
		return result;
	}
	size_t CompressedIndex::baseCount() const
	{
		assert(frozen());
		return pieces.baseCount();
	}
	size_t CompressedIndex::maxIndex() const
	{
		if (frozen()) return hierarchyTopDownFirst.size() + firstHierarchicalIndex;
		return pieceIndex.size() + hierarchyIndex.size();
	}
	size_t CompressedIndex::pieceCount() const
	{
		if (frozen()) return pieces.size();
		return pieceIndex.size();
	}
	bool CompressedIndex::frozen() const
	{
		return isfrozen;
	}
	void CompressedIndex::removeConstructionVariables(std::vector<VariableWidthIntVector>& indices, const size_t numThreads)
	{
		assert(!isfrozen);
		bitsPerIndex = ceil(log2(pieceIndex.size() + hierarchyIndex.size()));
		firstHierarchicalIndex = pieceIndex.size();
		size_t countBases = pieceIndex.countBases();
		pieces.setMaxStringLength(k+w);
		VariableWidthIntVector pieceReordering;
		pieceReordering.setWidth(ceil(log2(pieceIndex.size()*4+4)));
		pieceReordering.resize(pieceIndex.size()*4);
		pieces.reserve(pieceIndex.size(), countBases);
		{
			pieceIndex.iterateKeyValues([&pieceReordering, this](const std::string key, const size_t value)
			{
				size_t index = value / 2;
				assert(pieceReordering.get(index) == 0);
				pieceReordering.set(index, pieces.size());
				pieces.push_back(key);
			});
		}
		//assert(pieceReordering.size() == firstHierarchicalIndex);
		{
			decltype(pieceIndex) tmp;
			std::swap(tmp, pieceIndex);
		}
		std::vector<std::thread> threads;
		std::atomic<size_t> indexIndex;
		indexIndex = 0;
		std::mutex indexMutex;
		for (size_t threadi = 0; threadi < std::max((size_t)1, (size_t)numThreads-1); threadi++)
		{
			threads.emplace_back([this, &indexIndex, &indexMutex, &indices, &pieceReordering]()
			{
				while (true)
				{
					size_t i = indices.size();
					{
						std::lock_guard<std::mutex> lock { indexMutex };
						i = indexIndex;
						indexIndex += 1;
					}
					if (i >= indices.size()) break;
					VariableWidthIntVector replacement;
					replacement.setWidth(bitsPerIndex);
					std::vector<size_t> tmp;
					for (size_t j = 0; j < indices[i].size(); j++)
					{
						tmp.emplace_back(indices[i].get(j));
					}
					tmp = hierarchizeIndices(tmp);
					replacement.resize(tmp.size());
					for (size_t j = 0; j < tmp.size(); j++)
					{
						if (tmp[j] % 2 == 0)
						{
							replacement.set(j, pieceReordering.get(tmp[j]/2));
						}
						else
						{
							replacement.set(j, (tmp[j]-1) / 2 + firstHierarchicalIndex);
						}
					}
					std::swap(replacement, indices[i]);
				}
			});
		}
		if (numThreads == 1) threads[0].join();
		assert(bitsPerIndex > 0);
		hierarchyTopDownFirst.setWidth(bitsPerIndex);
		hierarchyTopDownSecond.setWidth(bitsPerIndex);
		hierarchyTopDownFirst.resize(hierarchyIndex.size());
		hierarchyTopDownSecond.resize(hierarchyIndex.size());
		hierarchyIndex.iterateKeyValues([this, &pieceReordering](std::pair<uint64_t, uint64_t> key, size_t value)
		{
			size_t index = (value-1) / 2;
			assert(index < hierarchyTopDownFirst.size());
			assert(hierarchyTopDownFirst.get(index) == 0);
			assert(hierarchyTopDownSecond.get(index) == 0);
			if (key.first % 2 == 0)
			{
				hierarchyTopDownFirst.set(index, pieceReordering.get(key.first / 2));
			}
			else
			{
				hierarchyTopDownFirst.set(index, (key.first-1) / 2 + firstHierarchicalIndex);
			}
			if (key.second % 2 == 0)
			{
				hierarchyTopDownSecond.set(index, pieceReordering.get(key.second / 2));
			}
			else
			{
				hierarchyTopDownSecond.set(index, (key.second-1) / 2 + firstHierarchicalIndex);
			}
		});
		if (numThreads > 1)
		{
			for (size_t threadi = 0; threadi < threads.size(); threadi++)
			{
				threads[threadi].join();
			}
		}
		{
			decltype(hierarchyIndex) tmp;
			std::swap(tmp, hierarchyIndex);
		}
		isfrozen = true;
	}
}
