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
		size_t at(const std::string& str) const;
		size_t size() const;
		template <typename F>
		void iterateValues(F callback) const
		{
			for (auto pair : len16Strings) callback(pair.second);
			for (auto pair : len32Strings) callback(pair.second);
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
			for (const auto& pair : biglenStrings)
			{
				callback(decodeStringToString(pair.first), pair.second);
			}
		}
	private:
		uint64_t encodeString(const std::string& str) const;
		std::string decodeString(uint64_t str) const;
		std::string encodeStringToString(const std::string& str) const;
		std::string decodeStringToString(std::string str) const;
		phmap::flat_hash_map<uint32_t, uint32_t> len16Strings;
		phmap::flat_hash_map<uint64_t, size_t> len32Strings;
		phmap::flat_hash_map<std::string, size_t> biglenStrings;
	};
}

#endif
