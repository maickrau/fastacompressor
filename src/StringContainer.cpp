#include <cassert>
#include <limits>
#include "StringContainer.h"

namespace FastaCompressor
{
	StringContainer::StringContainer() :
		bits(),
		smallOffsets(),
		bigOffsets(),
		maxStringLength(0),
		realSize(0)
	{
	}
	void StringContainer::setMaxStringLength(size_t length)
	{
		assert(size() == 0);
		maxStringLength = length;
	}
	size_t StringContainer::bigOffsetEveryNSmallOffsets() const
	{
		return (size_t)std::numeric_limits<uint16_t>::max() / maxStringLength;
	}
	void StringContainer::reserve(size_t countStrings, size_t countBps)
	{
		assert(size() == 0);
		assert(maxStringLength > 0);
		bits.resize((countBps+31) / 32, 0);
		smallOffsets.resize(countStrings);
		bigOffsets.resize((countStrings + bigOffsetEveryNSmallOffsets()-1)/bigOffsetEveryNSmallOffsets());
	}
	void StringContainer::push_back(const std::string& str)
	{
		size_t index = 0;
		if (size() > 0)
		{
			index = smallOffsets[size()-1] + bigOffsets[(size()-1) / bigOffsetEveryNSmallOffsets()];
		}
		for (size_t i = 0; i < str.size(); i++)
		{
			switch(str[i])
			{
			case 'a':
			case 'A':
				break;
			case 'c':
			case 'C':
				bits[index/32] += 1ull << ((size_t)(index%32)*2);
				break;
			case 'g':
			case 'G':
				bits[index/32] += 2ull << ((size_t)(index%32)*2);
				break;
			case 't':
			case 'T':
				bits[index/32] += 3ull << ((size_t)(index%32)*2);
				break;
			default:
				assert(false);
			}
			index += 1;
			assert(index/32 < bits.size());
		}
		if (size() % bigOffsetEveryNSmallOffsets() == 0)
		{
			bigOffsets[size() / bigOffsetEveryNSmallOffsets()] = index;
		}
		size_t smallpos = index - bigOffsets[size() / bigOffsetEveryNSmallOffsets()];
		assert(smallpos < (size_t)std::numeric_limits<uint16_t>::max());
		smallOffsets[size()] = smallpos;
		realSize += 1;
	}
	size_t StringContainer::baseCount() const
	{
		return smallOffsets.back() + bigOffsets.back();
	}
	std::string StringContainer::get(size_t index) const
	{
		assert(index < size());
		size_t startpos = 0;
		if (index > 0)
		{
			startpos = (size_t)smallOffsets[index-1] + bigOffsets[(index-1)/bigOffsetEveryNSmallOffsets()];
		}
		size_t endpos = (size_t)smallOffsets[index] + bigOffsets[(index)/bigOffsetEveryNSmallOffsets()];
		assert(endpos - startpos < maxStringLength);
		std::string result;
		result.resize(endpos - startpos);
		for (size_t i = 0; i < endpos-startpos; i++)
		{
			result[i] = "ACGT"[(bits[(startpos+i)/32] >> (((startpos+i)%32)*2)) & 3];
		}
		return result;
	}
	size_t StringContainer::size() const
	{
		return realSize;
	}
}