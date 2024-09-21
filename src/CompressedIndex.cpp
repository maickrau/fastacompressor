#include <iostream>
#include <chrono>
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
		std::vector<uint8_t> pieceType;
		shortPieces.resize(minimizerPositions.size()-1);
		longPieces.resize(minimizerPositions.size()-1);
		pieceType.resize(minimizerPositions.size()-1, std::numeric_limits<uint8_t>::max());
		for (size_t i = 1; i < minimizerPositions.size(); i++)
		{
			size_t start = minimizerPositions[i-1];
			size_t end = minimizerPositions[i];
			std::string kmer = seq.substr(start, end-start);
			if (kmer.size() <= 15)
			{
				pieceType[i-1] = 0;
				shortPieces[i-1] = pieceIndex.encodeString(kmer);
			}
			else if (kmer.size() <= 31)
			{
				pieceType[i-1] = 1;
				shortPieces[i-1] = pieceIndex.encodeString(kmer);
			}
			else if (kmer.size() <= 63)
			{
				pieceType[i-1] = 2;
				shortPieces[i-1] = pieceIndex.encodeString(kmer);
			}
			else
			{
				pieceType[i-1] = 3;
				longPieces[i-1] = pieceIndex.encodeStringToString(kmer);
			}
		}
		std::vector<size_t> result;
		result.reserve(minimizerPositions.size()-1);
		{
			std::lock_guard lock { pieceMutex };
			for (size_t i = 0; i < shortPieces.size(); i++)
			{
				size_t index = 0;
				size_t addedIndex = pieceIndex.size()*2;
				index = pieceIndex.getIndexOrSet(shortPieces[i], longPieces[i], pieceType[i], addedIndex);
				result.emplace_back(index);
			}
		}
		assert(result.size() == minimizerPositions.size()-1);
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
	void CompressedIndex::removeConstructionVariables(std::vector<VariableWidthIntVector>& indices)
	{
		assert(!isfrozen);
		bitsPerIndex = ceil(log2(pieceIndex.size() + hierarchyIndex.size()));
		firstHierarchicalIndex = pieceIndex.size();
		for (size_t i = 0; i < indices.size(); i++)
		{
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
					replacement.set(j, tmp[j]/2);
				}
				else
				{
					replacement.set(j, (tmp[j]-1) / 2 + firstHierarchicalIndex);
				}
			}
			std::swap(replacement, indices[i]);
		}
		assert(bitsPerIndex > 0);
		hierarchyTopDownFirst.setWidth(bitsPerIndex);
		hierarchyTopDownSecond.setWidth(bitsPerIndex);
		hierarchyTopDownFirst.resize(hierarchyIndex.size());
		hierarchyTopDownSecond.resize(hierarchyIndex.size());
		hierarchyIndex.iterateKeyValues([this](std::pair<uint64_t, uint64_t> key, size_t value)
		{
			size_t index = (value-1) / 2;
			assert(index < hierarchyTopDownFirst.size());
			assert(hierarchyTopDownFirst.get(index) == 0);
			assert(hierarchyTopDownSecond.get(index) == 0);
			if (key.first % 2 == 0)
			{
				hierarchyTopDownFirst.set(index, key.first / 2);
			}
			else
			{
				hierarchyTopDownFirst.set(index, (key.first-1) / 2 + firstHierarchicalIndex);
			}
			if (key.second % 2 == 0)
			{
				hierarchyTopDownSecond.set(index, key.second / 2);
			}
			else
			{
				hierarchyTopDownSecond.set(index, (key.second-1) / 2 + firstHierarchicalIndex);
			}
		});
		{
			decltype(hierarchyIndex) tmp;
			std::swap(tmp, hierarchyIndex);
		}
		size_t countBases = 0;
		pieceIndex.iterateKeyValues([&countBases](const std::string key, const size_t value)
		{
			countBases += key.size();
		});
		pieces.setMaxStringLength(k+w);
		VariableWidthIntVector pieceReordering;
		pieceReordering.setWidth(ceil(log2(pieceIndex.size()*2+1)));
		pieceReordering.resize(pieceIndex.size());
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
		std::vector<bool> valueFound;
		valueFound.resize(pieceIndex.size(), false);
		for (size_t i = 0; i < pieceReordering.size(); i++)
		{
			assert(!valueFound[pieceReordering.get(i)]);
			assert(pieceReordering.get(i) < valueFound.size());
			valueFound[pieceReordering.get(i)] = true;
		}
		assert(pieceReordering.size() == firstHierarchicalIndex);
		for (size_t i = 0; i < indices.size(); i++)
		{
			for (size_t j = 0; j < indices[i].size(); j++)
			{
				if (indices[i].get(j) < pieceReordering.size())
				{
					size_t newValue = pieceReordering.get(indices[i].get(j));
					indices[i].set(j, newValue);
					assert(indices[i].get(j) == newValue);
				}
			}
		}
		for (size_t i = 0; i < hierarchyTopDownFirst.size(); i++)
		{
			if (hierarchyTopDownFirst.get(i) < pieceReordering.size())
			{
				hierarchyTopDownFirst.set(i, pieceReordering.get(hierarchyTopDownFirst.get(i)));
			}
		}
		for (size_t i = 0; i < hierarchyTopDownSecond.size(); i++)
		{
			if (hierarchyTopDownSecond.get(i) < pieceReordering.size())
			{
				hierarchyTopDownSecond.set(i, pieceReordering.get(hierarchyTopDownSecond.get(i)));
			}
		}
		{
			decltype(pieceIndex) tmp;
			std::swap(tmp, pieceIndex);
		}
		isfrozen = true;
	}
}
