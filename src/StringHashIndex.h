#ifndef StringHashIndex_h
#define StringContainer_h

#include <cstdint>
#include <string>
#include <phmap.h>

namespace FastaCompressor
{
	class StringHashIndex
	{
	public:
		bool count(const std::string& str) const;
		void setIndex(const std::string& str, const size_t value);
		size_t getIndexOrNull(const __uint128_t shortPiece, const std::string& longPiece, const uint8_t pieceType) const;
		size_t getIndexOrSet(const __uint128_t shortPiece, const std::string& longPiece, const uint8_t pieceType, const size_t value);
		size_t getIndexOrSet(const std::string& str, const size_t value);
		size_t at(const std::string& str) const;
		size_t size() const;
		size_t typeCount(size_t type) const;
		size_t countBases() const;
		template <typename F>
		void iterateValues(F callback) const
		{
			for (auto pair : len16Strings) callback(pair.second);
			for (auto pair : len32Strings) callback(pair.second);
			for (auto pair : len64Strings) callback(pair.second);
			for (const auto& pair : biglenStrings) callback(pair.second);
		}
		template <typename F>
		void iterateKeyValues(F callback) const
		{
			for (auto pair : len16Strings)
			{
				std::string key = decodeString(pair.first);
				callback(key, pair.second);
			}
			for (auto pair : len32Strings)
			{
				std::string key = decodeString(pair.first);
				callback(key, pair.second);
			}
			for (auto pair : len64Strings)
			{
				std::string key = decodeString(pair.first);
				callback(key, pair.second);
			}
			for (const auto& pair : biglenStrings)
			{
				callback(decodeStringToString(pair.first), pair.second);
			}
		}
		__uint128_t encodeString(const std::string& str) const;
		std::string encodeStringToString(const std::string& str) const;
	private:
		size_t decodeStringToStringLength(const std::string& str) const;
		size_t decodeStringLength(__uint128_t val) const;
		std::string decodeString(__uint128_t str) const;
		std::string decodeStringToString(std::string str) const;
		phmap::flat_hash_map<uint32_t, size_t> len16Strings;
		phmap::flat_hash_map<uint64_t, size_t> len32Strings;
		phmap::flat_hash_map<__uint128_t, size_t> len64Strings;
		phmap::flat_hash_map<std::string, size_t> biglenStrings;
	};
}

#endif
