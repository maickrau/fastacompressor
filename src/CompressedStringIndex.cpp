#include <iostream>
#include "CompressedStringIndex.h"

namespace FastaCompressor
{
	CompressedStringIndex::CompressedStringIndex(size_t k, size_t w) :
		index(k, w),
		readNames(),
		readIndices()
	{
	}
	size_t CompressedStringIndex::addString(const std::string& readName, const std::string& readSequence)
	{
		size_t result = readNames.size();
		readNames.emplace_back(readName);
		readIndices.emplace_back(index.addString(readSequence));
		return result;
	}
	std::string CompressedStringIndex::getName(const size_t i) const
	{
		assert(i < readNames.size());
		return readNames[i];
	}
	std::string CompressedStringIndex::getSequence(const size_t i) const
	{
		assert(i < readIndices.size());
		assert(index.frozen());
		return index.getString(readIndices[i]);
	}
	void CompressedStringIndex::removeConstructionVariables()
	{
		for (size_t i = 0; i < readIndices.size(); i++)
		{
			readIndices[i] = index.convertToPostConstructionFormat(index.hierarchizeIndices(readIndices[i]));
		}
		index.removeConstructionVariables();
	}
	size_t CompressedStringIndex::size() const
	{
		return readNames.size();
	}
	void CompressedStringIndex::printSizeInformation() const
	{
		size_t bitsPerIndex = ceil(log2(index.maxIndex()));
		size_t countIndices = 0;
		for (size_t i = 0; i < readIndices.size(); i++)
		{
			countIndices += readIndices[i].size();
		}
		size_t bitsForPieces = countIndices * bitsPerIndex;
		size_t bitsForBases = index.baseCount() * 2;
		size_t bitsForHierarchy = bitsPerIndex * index.maxIndex();
		size_t totalBases = 0;
		for (size_t i = 0; i < readIndices.size(); i++)
		{
			totalBases += index.getString(readIndices[i]).size();
		}
		std::cerr << "bits per index: " << bitsPerIndex << std::endl;
		std::cerr << "remade bases " << totalBases << std::endl;
		std::cerr << "bits per bp: " << (double)(bitsForBases + bitsForPieces + bitsForHierarchy)/(double)totalBases << " (bases " << (double)bitsForBases/(double)totalBases << " pieces " << (double)bitsForPieces/(double)totalBases << " hierarchy " << (double)bitsForHierarchy/(double)totalBases << ")" << std::endl;
	}
}
