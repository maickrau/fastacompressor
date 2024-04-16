#ifndef StringContainer_h
#define StringContainer_h

#include <cstddef>
#include <cstdint>
#include <vector>
#include <string>

namespace FastaCompressor
{
	class StringContainer
	{
	public:
		StringContainer();
		void setMaxStringLength(size_t length);
		void reserve(size_t countStrings, size_t countBps);
		void push_back(const std::string& str);
		std::string get(size_t index) const;
		size_t pieceSize(size_t index) const;
		size_t size() const;
		size_t baseCount() const;
	private:
		size_t bigOffsetEveryNSmallOffsets() const;
		std::vector<uint64_t> bits;
		std::vector<uint16_t> smallOffsets; // end position of string
		std::vector<size_t> bigOffsets; // end position of string
		size_t maxStringLength;
		size_t realSize;
	};
}

#endif
