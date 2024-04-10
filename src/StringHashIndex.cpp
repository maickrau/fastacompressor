#include "StringHashIndex.h"

namespace FastaCompressor
{
	bool StringHashIndex::count(const std::string& str) const
	{
		if (str.size() <= 31)
		{
			uint64_t encodedString = encodeString(str);
			if (str.size() <= 15)
			{
				assert(encodedString < (size_t)std::numeric_limits<uint32_t>::max());
				if (len16Strings.count(encodedString) == 1) return 1;
			}
			if (len32Strings.count(encodedString) == 1) return 1;
		}
		if (biglenStrings.count(str) == 1) return 1;
		return 0;
	}
	size_t StringHashIndex::at(const std::string& str) const
	{
		if (str.size() <= 31)
		{
			uint64_t encodedString = encodeString(str);
			if (str.size() <= 15)
			{
				assert(encodedString < (size_t)std::numeric_limits<uint32_t>::max());
				if (len16Strings.count(encodedString) == 1) return len16Strings.at(encodedString);
			}
			if (len32Strings.count(encodedString) == 1) return len32Strings.at(encodedString);
		}
		return biglenStrings.at(str);
	}
	void StringHashIndex::setIndex(const std::string& str, const size_t value)
	{
		if (value < (size_t)std::numeric_limits<uint32_t>::max() && str.size() <= 15)
		{
			uint64_t encodedString = encodeString(str);
			assert(encodedString < (size_t)std::numeric_limits<uint32_t>::max());
			assert(len16Strings.count(encodedString) == 0);
			len16Strings[encodedString] = value;
			return;
		}
		if (str.size() <= 31)
		{
			uint64_t encodedString = encodeString(str);
			assert(len32Strings.count(encodedString) == 0);
			len32Strings[encodedString] = value;
			return;
		}
		assert(biglenStrings.count(str) == 0);
		biglenStrings[str] = value;
	}
	size_t StringHashIndex::size() const
	{
		return len16Strings.size() + len32Strings.size() + biglenStrings.size();
	}
	uint64_t StringHashIndex::encodeString(const std::string& str) const
	{
		uint64_t result = 1;
		for (size_t i = str.size()-1; i < str.size(); i--)
		{
			result <<= 2;
			switch(str[i])
			{
			case 'A':
				break;
			case 'C':
				result += 1;
				break;
			case 'G':
				result += 2;
				break;
			case 'T':
				result += 3;
				break;
			}
		}
		return result;
	}
	std::string StringHashIndex::decodeString(uint64_t val) const
	{
		std::string result;
		while (val != 1)
		{
			result += "ACGT"[val & 3];
			val >>= 2;
		}
		return result;
	}
}