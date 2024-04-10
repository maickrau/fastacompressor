#ifndef CompressedStringIndex_h
#define CompressedStringIndex_h

#include "CompressedIndex.h"

namespace FastaCompressor
{
	class CompressedStringIndex
	{
	public:
		CompressedStringIndex(size_t k, size_t w);
		size_t addString(const std::string& readName, const std::string& readSequence);
		std::string getName(const size_t index) const;
		std::string getSequence(const size_t index) const;
		void removeConstructionVariables();
		size_t size() const;
		void printSizeInformation() const;
	private:
		CompressedIndex index;
		std::vector<std::string> readNames;
		std::vector<std::vector<size_t>> readIndicesUnbuilt;
		std::vector<VariableWidthIntVector> readIndicesFinished;
	};
}

#endif
