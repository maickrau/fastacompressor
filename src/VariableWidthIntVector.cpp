#include <cassert>
#include "VariableWidthIntVector.h"

namespace FastaCompressor
{
	VariableWidthIntVector::VariableWidthIntVector() :
		data(),
		realSize(0),
		bitWidth(0),
		mask(0)
	{
	}
	void VariableWidthIntVector::setWidth(size_t bits)
	{
		assert(bits >= 1);
		assert(bits <= 64);
		assert(size() == 0);
		bitWidth = bits;
		if (bits == 64)
		{
			mask = -1;
		}
		else
		{
			mask = (1ull << bits) - 1;
		}
	}
	void VariableWidthIntVector::resize(size_t newSize)
	{
		if (size() >= newSize)
		{
			realSize = newSize;
		}
		else
		{
			data.resize((newSize*width()+63)/64, 0);
			realSize = newSize;
		}
	}
	size_t VariableWidthIntVector::get(size_t i) const
	{
		assert(i < realSize);
		size_t bitPos = i * width();
		size_t index = bitPos / 64;
		assert(index < data.size());
		size_t offset = bitPos % 64;
		if (64 - offset >= width())
		{
			return (data[index] >> offset) & mask;
		}
		size_t result = 0;
		size_t bitsLeft = width() - (64 - offset);
		result = (data[index] >> offset) << bitsLeft;
		assert(index+1 < data.size());
		result += data[index+1] & (mask >> (width()-bitsLeft));
		return result;
	}
	void VariableWidthIntVector::set(size_t i, size_t value)
	{
		assert(i < realSize);
		assert((value & mask) == value);
		size_t bitPos = i * width();
		size_t index = bitPos / 64;
		assert(index < data.size());
		size_t offset = bitPos % 64;
		if (64 - offset >= width())
		{
			data[index] &= ~(mask << offset);
			data[index] += value << offset;
			return;
		}
		data[index] &= ~(mask << offset);
		data[index] += (value << (64 - width())) & (mask << offset);
		size_t bitsLeft = width() - (64 - offset);
		assert(index+1 < data.size());
		data[index+1] &= ~(mask >> (width()-bitsLeft));
		data[index+1] += value & (mask >> (width()-bitsLeft));
	}
	size_t VariableWidthIntVector::size() const
	{
		return realSize;
	}
	size_t VariableWidthIntVector::width() const
	{
		return bitWidth;
	}
}
