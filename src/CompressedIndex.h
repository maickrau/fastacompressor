#ifndef CompressedIndex_h
#define CompressedIndex_h

#include <mutex>
#include <vector>
#include <string>
#include <condition_variable>
#include <phmap.h>
#include "VariableWidthIntVector.h"
#include "StringContainer.h"
#include "StringHashIndex.h"
#include "HierarchyIndex.h"

namespace FastaCompressor
{
	class CompressedIndex
	{
	public:
		CompressedIndex(size_t k, size_t w);
		std::vector<size_t> addString(const std::string& str);
		std::vector<size_t> hierarchizeIndices(const std::vector<size_t>& indices) const;
		VariableWidthIntVector hierarchizeIndices(const VariableWidthIntVector& indices) const;
		std::string getString(const std::vector<size_t>& indices) const;
		std::string getString(const VariableWidthIntVector& indices) const;
		std::string getSubstring(const std::vector<size_t>& indices, const size_t startPos, const size_t length) const;
		std::string getSubstring(const VariableWidthIntVector& indices, const size_t startPos, const size_t length) const;
		size_t maxIndex() const;
		size_t pieceCount() const;
		size_t baseCount() const;
		void removeConstructionVariables(std::vector<VariableWidthIntVector>& indices);
		bool frozen() const;
	private:
		std::vector<size_t> segmentFastaToPieces(const std::string& sequence);
		HierarchyIndex hierarchyIndex;
		StringHashIndex pieceIndex;
		VariableWidthIntVector hierarchyTopDownFirst;
		VariableWidthIntVector hierarchyTopDownSecond;
		StringContainer pieces;
		size_t firstHierarchicalIndex;
		bool isfrozen;
		size_t bitsPerIndex;
		size_t k;
		size_t w;
		std::mutex pieceMutex;
		std::atomic<size_t> pieceReaderCount;
		std::condition_variable pieceConditionVariable;
		std::mutex hierarchyMutex;
		std::atomic<size_t> hierarchyReaderCount;
		std::condition_variable hierarchyConditionVariable;
	};
}

#endif
