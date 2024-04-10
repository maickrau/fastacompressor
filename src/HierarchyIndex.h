#ifndef HierarchyIndex_h
#define HierarchyIndex_h

#include <cstdint>
#include <tuple>
#include <phmap.h>

namespace FastaCompressor
{
	class HierarchyIndex
	{
	public:
		size_t count(std::pair<uint64_t, uint64_t> key) const;
		size_t at(std::pair<uint64_t, uint64_t> key) const;
		void set(std::pair<uint64_t, uint64_t> key, size_t value);
		size_t size() const;
		template <typename F>
		void iterateKeyValues(F callback) const
		{
			for (auto pair : smalls)
			{
				callback(std::pair<uint64_t, uint64_t> { pair.first.first, pair.first.second }, (uint64_t)pair.second);
			}
			for (auto pair : bigs)
			{
				callback(pair.first, pair.second);
			}
		}
	private:
		phmap::flat_hash_map<std::pair<uint32_t, uint32_t>, uint32_t> smalls;
		phmap::flat_hash_map<std::pair<uint64_t, uint64_t>, uint64_t> bigs;
	};
}

#endif
