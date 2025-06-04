#include <iostream>
#include "CompressedStringIndex.h"
#include "Serializer.h"

namespace FastaCompressor
{
	CompressedStringIndex::CompressedStringIndex(size_t k, size_t w) :
		index(k, w),
		readNames(),
		readIndices()
	{
		indexMutex = std::unique_ptr<std::mutex> { new std::mutex };
	}
	size_t CompressedStringIndex::addString(const std::string& readName, const std::string& readSequence)
	{
		size_t result;
		{
			std::lock_guard lock { *indexMutex };
			result = readNames.size();
			readNames.emplace_back();
			readIndices.emplace_back();
			readNames[result] = readName;
		}
		std::vector<size_t> indices = index.addString(readSequence);
		assert(indices.size() >= 1);
		size_t maxIndex = 1;
		for (size_t val : indices)
		{
			maxIndex = std::max(maxIndex, val);
		}
		VariableWidthIntVector tmp;
		tmp.setWidth(ceil(log2(maxIndex+1)));
		tmp.resize(indices.size());
		for (size_t i = 0; i < indices.size(); i++)
		{
			tmp.set(i, indices[i]);
		}
		{
			std::lock_guard lock { *indexMutex };
			std::swap(tmp, readIndices[result]);
		}
		return result;
	}
	size_t CompressedStringIndex::addString(const size_t readID, const std::string& readName, const std::string& readSequence)
	{
		size_t result;
		result = readID;
		{
			std::lock_guard lock { *indexMutex };
			while (readID >= readNames.size())
			{
				readNames.emplace_back();
				readIndices.emplace_back();
			}
			assert(readNames[readID].size() == 0);
			assert(readIndices[readID].size() == 0);
			readNames[result] = readName;
		}
		std::vector<size_t> indices = index.addString(readSequence);
		assert(indices.size() >= 1);
		size_t maxIndex = 1;
		for (size_t val : indices)
		{
			maxIndex = std::max(maxIndex, val);
		}
		{
			std::lock_guard lock { *indexMutex };
			readIndices[result].setWidth(ceil(log2(maxIndex+1)));
			readIndices[result].resize(indices.size());
			for (size_t i = 0; i < indices.size(); i++)
			{
				readIndices[result].set(i, indices[i]);
			}
		}
		return result;
	}
	std::string CompressedStringIndex::getName(const size_t i) const
	{
		assert(i < readNames.size());
		return readNames[i];
	}
	std::string CompressedStringIndex::getSubstring(const size_t i, const size_t startPos, const size_t length) const
	{
		assert(i < readIndices.size());
		assert(index.frozen());
		return index.getSubstring(readIndices[i], startPos, length);
	}
	std::string CompressedStringIndex::getSequence(const size_t i) const
	{
		assert(i < readIndices.size());
		assert(index.frozen());
		return index.getString(readIndices[i]);
	}
	void CompressedStringIndex::removeConstructionVariables(const size_t numThreads)
	{
		index.removeConstructionVariables(readIndices, numThreads);
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
		size_t bitsForBases = index.baseCount() * 2 + index.pieceCount() * 17;
		size_t bitsForHierarchy = bitsPerIndex * index.maxIndex();
		size_t totalBases = 0;
		for (size_t i = 0; i < readIndices.size(); i++)
		{
			totalBases += index.getString(readIndices[i]).size();
		}
		std::cerr << "bits per index: " << bitsPerIndex << std::endl;
		std::cerr << "remade bases " << totalBases << std::endl;
		std::cerr << "bits per bp: " << (double)(bitsForBases + bitsForPieces + bitsForHierarchy)/(double)totalBases << " (bases " << (double)bitsForBases/(double)totalBases << " pieces " << (double)bitsForPieces/(double)totalBases << " hierarchy " << (double)bitsForHierarchy/(double)totalBases << ")" << std::endl;
		std::cerr << "expected size: " << (bitsForBases + bitsForPieces + bitsForHierarchy) / 8.0 / 1024.0 / 1024.0 / 1024.0 << " gb" << std::endl;
		std::cerr << "bases in piece index " << index.baseCount() << " construction bytes ~" << (index.baseCount())/1024.0/1024.0/1024.0 << " gb" << std::endl;
		std::cerr << "count in piece index " << index.pieceCount() << " construction bytes ~" <<  (index.pieceCount()*64)/1024.0/1024.0/1024.0 << " gb" << std::endl;
		std::cerr << "hierarchy index size " << index.maxIndex() << " construction bytes ~" << (index.maxIndex()*48)/1024.0/1024.0/1024.0 << " gb" << std::endl;
	}
	void CompressedStringIndex::writeToStream(std::ostream& stream) const
	{
		assert(index.frozen());
		index.writeToStream(stream);
		Serializer::writeUint64_t(readNames.size(), stream);
		for (size_t i = 0; i < readNames.size(); i++)
		{
			Serializer::writeString(readNames[i], stream);
		}
		Serializer::writeUint64_t(readIndices.size(), stream);
		size_t width = 0;
		if (readIndices.size() > 0) width = readIndices[0].width();
		Serializer::writeUint64_t(width, stream);
		if (readIndices.size() > 0)
		{
			for (size_t i = 0; i < readIndices.size(); i++)
			{
				assert(readIndices[i].width() == width);
				readIndices[i].writeToStream(stream);
			}
		}
	}
	CompressedStringIndex CompressedStringIndex::loadFromStream(std::istream& stream)
	{
		CompressedStringIndex result { 0, 0 };
		result.index = CompressedIndex::loadFromStream(stream);
		size_t readNameCount = Serializer::readUint64_t(stream);
		result.readNames.reserve(readNameCount);
		for (size_t i = 0; i < readNameCount; i++)
		{
			result.readNames.emplace_back(Serializer::readString(stream));
		}
		size_t readIndicesCount = Serializer::readUint64_t(stream);
		size_t readIndicesWidth = Serializer::readUint64_t(stream);
		result.readIndices.reserve(readIndicesCount);
		for (size_t i = 0; i < readIndicesCount; i++)
		{
			result.readIndices.emplace_back(VariableWidthIntVector::loadFromStream(readIndicesWidth, stream));
		}
		return result;
	}
}
