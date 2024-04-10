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
		auto fixstr = encodeStringToString(str);
		if (biglenStrings.count(fixstr) == 1) return 1;
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
		auto fixstr = encodeStringToString(str);
		return biglenStrings.at(fixstr);
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
		auto fixstr = encodeStringToString(str);
		assert(biglenStrings.count(fixstr) == 0);
		biglenStrings[fixstr] = value;
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
	std::string StringHashIndex::encodeStringToString(const std::string& str) const
	{
		std::string result;
		result.push_back(1);
		for (size_t i = 0; i < str.size(); i++)
		{
			size_t c = 0;
			switch(str[i])
			{
			case 'A':
				c = 0;
				break;
			case 'C':
				c = 1;
				break;
			case 'G':
				c = 2;
				break;
			case 'T':
				c = 3;
				break;
			default:
				assert(false);
			}
			result.back() <<= 2;
			result.back() += c;
			if (i % 4 == 3) result.push_back(1);
		}
		return result;
	}
	std::string StringHashIndex::decodeStringToString(std::string str) const
	{
		std::string result;
		for (size_t i = 0; i < str.size()-1; i++)
		{
			result += "ACGT"[(str[i] >> 6) & 3];
			result += "ACGT"[(str[i] >> 4) & 3];
			result += "ACGT"[(str[i] >> 2) & 3];
			result += "ACGT"[str[i] & 3];
		}
		uint8_t last = str.back();
		size_t pickedInLast = 0;
		while (last != 1)
		{
			result += "ACGT"[last & 3];
			pickedInLast += 1;
			last >>= 2;
			assert(last != 0);
		}
		assert(pickedInLast <= 3);
		std::reverse(result.end() - pickedInLast, result.end());
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