#ifndef CompressedIndex_h
#define CompressedIndex_h

#include <fstream>
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
		struct TemporaryConstructionInfo
		{
		public:
			TemporaryConstructionInfo();
			HierarchyIndex hierarchyIndex;
			StringHashIndex pieceIndex;
			std::array<std::mutex, 4> pieceMutex;
			std::array<std::atomic<size_t>, 4> pieceReaderCount;
			std::array<std::condition_variable, 4> pieceConditionVariable;
			std::mutex hierarchyMutex;
			std::atomic<size_t> hierarchyReaderCount;
			std::condition_variable hierarchyConditionVariable;
		};
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
		void removeConstructionVariables(std::vector<VariableWidthIntVector>& indices, const size_t numThreads);
		bool frozen() const;
		void writeToStream(std::ostream& stream) const;
		static CompressedIndex loadFromStream(std::istream& stream);
	private:
		std::vector<size_t> segmentFastaToPieces(const std::string& sequence);
		VariableWidthIntVector hierarchyTopDownFirst;
		VariableWidthIntVector hierarchyTopDownSecond;
		StringContainer pieces;
		size_t firstHierarchicalIndex;
		bool isfrozen;
		size_t bitsPerIndex;
		size_t k;
		size_t w;
		std::unique_ptr<TemporaryConstructionInfo> temps;
	};
}

#endif
