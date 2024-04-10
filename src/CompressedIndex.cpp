#include "CompressedIndex.h"
#include "MinimizerIterator.h"

namespace FastaCompressor
{
	CompressedIndex::CompressedIndex(size_t k, size_t w) :
		seenOnce(),
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
		std::vector<size_t> result;
		for (size_t i = 1; i < minimizerPositions.size(); i++)
		{
			size_t start = minimizerPositions[i-1];
			size_t end = minimizerPositions[i];
			std::string kmer = seq.substr(start, end-start);
			size_t index = 0;
			if (pieceIndex.count(kmer) == 0)
			{
				index = pieceIndex.size() + (1ull << 63ull);
				pieceIndex[kmer] = index;
			}
			else
			{
				index = pieceIndex.at(kmer);
			}
			assert(index & (1ull << 63ull));
			result.emplace_back(index);
		}
		return result;
	}
	std::vector<size_t> CompressedIndex::addString(const std::string& str)
	{
		auto indices = segmentFastaToPieces(str);
		auto fixedIndices = hierarchizeIndices(indices);
		std::vector<size_t> result;
		for (size_t i = 0; i+1 < fixedIndices.size(); i += 2)
		{
			if (i > 0 && hierarchyIndex.count(std::make_pair(fixedIndices[i], fixedIndices[i+1])) == 1)
			{
				result.emplace_back(hierarchyIndex.at(std::make_pair(fixedIndices[i], fixedIndices[i+1])));
			}
			else if (seenOnce.count(fixedIndices[i]) == 1 && seenOnce.count(fixedIndices[i+1]) == 1)
			{
				assert(hierarchyIndex.count(std::make_pair(fixedIndices[i], fixedIndices[i+1])) == 0);
				size_t index = hierarchyIndex.size();
				hierarchyIndex[std::make_pair(fixedIndices[i], fixedIndices[i+1])] = index;
				result.emplace_back(index);
			}
			else
			{
				seenOnce.insert(fixedIndices[i]);
				seenOnce.insert(fixedIndices[i+1]);
				result.emplace_back(fixedIndices[i]);
				result.emplace_back(fixedIndices[i+1]);
			}
		}
		if (fixedIndices.size() % 2 == 1) result.emplace_back(fixedIndices.back());
		return result;
	}
	std::vector<size_t> CompressedIndex::hierarchizeIndices(const std::vector<size_t>& indices) const
	{
		std::vector<size_t> result;
		for (size_t i = 0; i < indices.size(); i++)
		{
			result.emplace_back(indices[i]);
			assert(indices[i] < hierarchyIndex.size() || (indices[i] & (1ull << 63ull)));
			while (result.size() >= 2 && hierarchyIndex.count(std::make_pair(result[result.size()-2], result[result.size()-1])) == 1)
			{
				size_t replacement = hierarchyIndex.at(std::make_pair(result[result.size()-2], result[result.size()-1]));
				assert(replacement < hierarchyIndex.size());
				result.pop_back();
				result.pop_back();
				result.emplace_back(replacement);
			}
		}
		return result;
	}
	VariableWidthIntVector CompressedIndex::convertToPostConstructionFormat(const std::vector<size_t>& indices) const
	{
		assert(!frozen());
		VariableWidthIntVector result;
		result.setWidth(ceil(log2(pieceIndex.size() + hierarchyIndex.size())));
		result.resize(indices.size());
		for (size_t i = 0; i < result.size(); i++)
		{
			if (indices[i] & (1ull << 63ull))
			{
				result.set(i, indices[i] ^ (1ull << 63ull));
			}
			else
			{
				assert(indices[i] < hierarchyIndex.size());
				result.set(i, indices[i] + pieceIndex.size());
			}
		}
		return result;
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
		size_t result = 0;
		if (frozen())
		{
			return pieces.baseCount();
		}
		else
		{
			for (const auto& pair : pieceIndex)
			{
				result += pair.first.size();
			}
		}
		return result;
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
	void CompressedIndex::removeConstructionVariables()
	{
		assert(!isfrozen);
		bitsPerIndex = ceil(log2(pieceIndex.size() + hierarchyIndex.size()));
		assert(bitsPerIndex > 0);
		{
			decltype(seenOnce) tmp;
			std::swap(tmp, seenOnce);
		}
		size_t countBases = 0;
		for (const auto& pair : pieceIndex)
		{
			countBases += pair.first.size();
		}
		pieces.setMaxStringLength(k+w);
		pieces.reserve(pieceIndex.size(), countBases);
		{
			std::vector<std::string> tmp;
			tmp.resize(pieceIndex.size());
			for (const auto& pair : pieceIndex)
			{
				assert(pair.second & (1ull << 63ull));
				size_t index = (pair.second ^ (1ull << 63ull));
				assert(index < tmp.size());
				assert(tmp[index] == "");
				tmp[index] = pair.first;
			}
			for (size_t i = 0; i < tmp.size(); i++)
			{
				assert(tmp[i] != "");
				pieces.push_back(tmp[i]);
			}
		}
		{
			decltype(pieceIndex) tmp;
			std::swap(tmp, pieceIndex);
		}
		firstHierarchicalIndex = pieces.size();
		hierarchyTopDownFirst.setWidth(bitsPerIndex);
		hierarchyTopDownSecond.setWidth(bitsPerIndex);
		hierarchyTopDownFirst.resize(hierarchyIndex.size());
		hierarchyTopDownSecond.resize(hierarchyIndex.size());
		for (auto pair : hierarchyIndex)
		{
			assert(pair.second < hierarchyTopDownFirst.size());
			assert(hierarchyTopDownFirst.get(pair.second) == 0);
			assert(hierarchyTopDownSecond.get(pair.second) == 0);
			if (pair.first.first & (1ull << 63ull))
			{
				hierarchyTopDownFirst.set(pair.second, pair.first.first ^ (1ull << 63ull));
			}
			else
			{
				hierarchyTopDownFirst.set(pair.second, pair.first.first + firstHierarchicalIndex);
			}
			if (pair.first.second & (1ull << 63ull))
			{
				hierarchyTopDownSecond.set(pair.second, pair.first.second ^ (1ull << 63ull));
			}
			else
			{
				hierarchyTopDownSecond.set(pair.second, pair.first.second + firstHierarchicalIndex);
			}
		}
		{
			decltype(hierarchyIndex) tmp;
			std::swap(tmp, hierarchyIndex);
		}
		isfrozen = true;
	}
}
