#ifndef VariableWidthIntVector_h
#define VariableWidthIntVector_h

#include <cstddef>
#include <cstdint>
#include <vector>

namespace FastaCompressor
{
	class VariableWidthIntVector
	{
	public:
		VariableWidthIntVector();
		void setWidth(size_t bits);
		void resize(size_t newSize);
		size_t get(size_t index) const;
		void set(size_t index, size_t value);
		size_t size() const;
		size_t width() const;
	private:
		std::vector<uint64_t> data;
		size_t realSize;
		size_t bitWidth;
		size_t mask;
	};
}

#endif
