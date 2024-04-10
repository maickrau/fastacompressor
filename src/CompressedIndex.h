#ifndef CompressedIndex_h
#define CompressedIndex_h

#include <vector>
#include <string>
#include <phmap.h>
#include "VariableWidthIntVector.h"

namespace FastaCompressor
{
	class CompressedIndex
	{
	public:
		CompressedIndex(size_t k, size_t w);
		std::vector<size_t> addString(const std::string& str);
		std::vector<size_t> hierarchizeIndices(const std::vector<size_t>& indices) const;
		VariableWidthIntVector convertToPostConstructionFormat(const std::vector<size_t>& indices) const;
		std::string getString(const std::vector<size_t>& indices) const;
		std::string getString(const VariableWidthIntVector& indices) const;
		size_t maxIndex() const;
		size_t pieceCount() const;
		size_t baseCount() const;
		void removeConstructionVariables();
		bool frozen() const;
	private:
		std::vector<size_t> segmentFastaToPieces(const std::string& sequence);
		phmap::flat_hash_set<size_t> seenOnce;
		phmap::flat_hash_map<std::pair<size_t, size_t>, size_t> hierarchyIndex;
		phmap::flat_hash_map<std::string, size_t> pieceIndex;
		VariableWidthIntVector hierarchyTopDownFirst;
		VariableWidthIntVector hierarchyTopDownSecond;
		std::vector<std::string> pieces;
		size_t firstHierarchicalIndex;
		bool isfrozen;
		size_t bitsPerIndex;
		size_t k;
		size_t w;
	};
}

#endif
