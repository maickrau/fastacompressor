#include "StringHashIndex.h"

namespace std
{
	template <> struct hash<__uint128_t>
	{
		size_t operator()(__uint128_t x) const
		{
			return std::hash<uint64_t>{}((uint64_t)x ^ (uint64_t)(x >> 64));
		}
	};
}

namespace FastaCompressor
{
	size_t StringHashIndex::countBases() const
	{
		size_t result = 0;
		for (auto pair : len16Strings)
		{
			result += decodeStringLength(pair.first);
		}
		for (auto pair : len32Strings)
		{
			result += decodeStringLength(pair.first);
		}
		for (auto pair : len64Strings)
		{
			result += decodeStringLength(pair.first);
		}
		for (const auto& pair : biglenStrings)
		{
			result += decodeStringToStringLength(pair.first);
		}
		return result;
	}
	bool StringHashIndex::count(const std::string& str) const
	{
		if (str.size() <= 15)
		{
			__uint128_t encodedString = encodeString(str);
			assert(encodedString < (__uint128_t)std::numeric_limits<uint32_t>::max());
			return len16Strings.count(encodedString);
		}
		if (str.size() <= 31)
		{
			__uint128_t encodedString = encodeString(str);
			assert(encodedString < (__uint128_t)std::numeric_limits<uint64_t>::max());
			return len32Strings.count(encodedString);
		}
		if (str.size() <= 63)
		{
			__uint128_t encodedString = encodeString(str);
			return len64Strings.count(encodedString);
		}
		auto fixstr = encodeStringToString(str);
		if (biglenStrings.count(fixstr) == 1) return 1;
		return 0;
	}
	size_t StringHashIndex::at(const std::string& str) const
	{
		if (str.size() <= 15)
		{
			__uint128_t encodedString = encodeString(str);
			assert(encodedString < (__uint128_t)std::numeric_limits<uint32_t>::max());
			return len16Strings.at(encodedString);
		}
		if (str.size() <= 31)
		{
			__uint128_t encodedString = encodeString(str);
			assert(encodedString < (__uint128_t)std::numeric_limits<uint64_t>::max());
			return len32Strings.at(encodedString);
		}
		if (str.size() <= 63)
		{
			__uint128_t encodedString = encodeString(str);
			return len64Strings.at(encodedString);
		}
		auto fixstr = encodeStringToString(str);
		return biglenStrings.at(fixstr);
	}
	size_t StringHashIndex::getIndexOrNull(const __uint128_t shortPiece, const std::string& longPiece, const uint8_t pieceType) const
	{
		switch(pieceType)
		{
		case 0:
			assert(shortPiece < (__uint128_t)std::numeric_limits<uint32_t>::max());
			if (len16Strings.count(shortPiece) == 1) return len16Strings.at(shortPiece);
			return std::numeric_limits<size_t>::max();
		case 1:
			assert(shortPiece < (__uint128_t)std::numeric_limits<uint64_t>::max());
			if (len32Strings.count(shortPiece) == 1) return len32Strings.at(shortPiece);
			return std::numeric_limits<size_t>::max();
		case 2:
			if (len64Strings.count(shortPiece) == 1) return len64Strings.at(shortPiece);
			return std::numeric_limits<size_t>::max();
		case 3:
			{
				auto found = biglenStrings.find(longPiece);
				if (found != biglenStrings.end()) return found->second;
				return std::numeric_limits<size_t>::max();
			}
		default:
			assert(false);
		}
	}
	size_t StringHashIndex::getIndexOrSet(const __uint128_t shortPiece, const std::string& longPiece, const uint8_t pieceType, const size_t value)
	{
		switch(pieceType)
		{
		case 0:
			assert(shortPiece < (__uint128_t)std::numeric_limits<uint32_t>::max());
			if (len16Strings.count(shortPiece) == 1) return len16Strings.at(shortPiece);
			len16Strings[shortPiece] = value;
			return value;
		case 1:
			assert(shortPiece < (__uint128_t)std::numeric_limits<uint64_t>::max());
			if (len32Strings.count(shortPiece) == 1) return len32Strings.at(shortPiece);
			len32Strings[shortPiece] = value;
			return value;
		case 2:
			if (len64Strings.count(shortPiece) == 1) return len64Strings.at(shortPiece);
			len64Strings[shortPiece] = value;
			return value;
		case 3:
			{
				auto found = biglenStrings.find(longPiece);
				if (found != biglenStrings.end()) return found->second;
				biglenStrings[longPiece] = value;
				return value;
			}
		default:
			assert(false);
		}
	}
	size_t StringHashIndex::getIndexOrSet(const std::string& str, const size_t value)
	{
		if (str.size() <= 15)
		{
			__uint128_t encodedString = encodeString(str);
			assert(encodedString < (__uint128_t)std::numeric_limits<uint32_t>::max());
			if (len16Strings.count(encodedString) == 1) return len16Strings.at(encodedString);
			len16Strings[encodedString] = value;
			return value;
		}
		if (str.size() <= 31)
		{
			__uint128_t encodedString = encodeString(str);
			assert(encodedString < (__uint128_t)std::numeric_limits<uint64_t>::max());
			if (len32Strings.count(encodedString) == 1) return len32Strings.at(encodedString);
			len32Strings[encodedString] = value;
			return value;
		}
		if (str.size() <= 63)
		{
			__uint128_t encodedString = encodeString(str);
			if (len64Strings.count(encodedString) == 1) return len64Strings.at(encodedString);
			len64Strings[encodedString] = value;
			return value;
		}
		auto fixstr = encodeStringToString(str);
		auto found = biglenStrings.find(fixstr);
		if (found != biglenStrings.end()) return found->second;
		biglenStrings[fixstr] = value;
		return value;
	}
	void StringHashIndex::setIndex(const std::string& str, const size_t value)
	{
		if (str.size() <= 15)
		{
			__uint128_t encodedString = encodeString(str);
			assert(encodedString < (__uint128_t)std::numeric_limits<uint32_t>::max());
			assert(len16Strings.count(encodedString) == 0);
			len16Strings[encodedString] = value;
			return;
		}
		if (str.size() <= 31)
		{
			__uint128_t encodedString = encodeString(str);
			assert(encodedString < (__uint128_t)std::numeric_limits<uint64_t>::max());
			assert(len32Strings.count(encodedString) == 0);
			len32Strings[encodedString] = value;
			return;
		}
		if (str.size() <= 63)
		{
			__uint128_t encodedString = encodeString(str);
			assert(len64Strings.count(encodedString) == 0);
			len64Strings[encodedString] = value;
			return;
		}
		auto fixstr = encodeStringToString(str);
		assert(biglenStrings.count(fixstr) == 0);
		biglenStrings[fixstr] = value;
	}
	size_t StringHashIndex::size() const
	{
		return len16Strings.size() + len32Strings.size() + len64Strings.size() + biglenStrings.size();
	}
	__uint128_t StringHashIndex::encodeString(const std::string& str) const
	{
		__uint128_t result = 1;
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
	size_t StringHashIndex::decodeStringLength(__uint128_t val) const
	{
		size_t result = 0;
		while (val != 1)
		{
			result += 1;
			val >>= 2;
		}
		return result;
	}
	size_t StringHashIndex::decodeStringToStringLength(const std::string& str) const
	{
		size_t result = (str.size()-1) * 4;
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
		return result + pickedInLast;
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
	std::string StringHashIndex::decodeString(__uint128_t val) const
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