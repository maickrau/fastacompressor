#ifndef CompressedStringIndex_h
#define CompressedStringIndex_h

#include <fstream>
#include <mutex>
#include <fstream>
#include "CompressedIndex.h"

namespace FastaCompressor
{
	class CompressedStringIndex
	{
	public:
		CompressedStringIndex(size_t k, size_t w);
		size_t addString(const std::string& readName, const std::string& readSequence);
		size_t addString(const size_t index, const std::string& readName, const std::string& readSequence);
		std::string getName(const size_t index) const;
		std::string getSequence(const size_t index) const;
		std::string getSubstring(const size_t index, const size_t startPos, const size_t length) const;
		void removeConstructionVariables(const size_t numThreads);
		size_t size() const;
		void printSizeInformation() const;
		void writeToStream(std::ostream& stream) const;
		static CompressedStringIndex loadFromStream(std::istream& stream);
	private:
		CompressedIndex index;
		std::vector<std::string> readNames;
		std::vector<VariableWidthIntVector> readIndices;
		std::unique_ptr<std::mutex> indexMutex;
	};
}

#endif
