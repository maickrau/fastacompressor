#include "HierarchyIndex.h"

namespace FastaCompressor
{
	size_t HierarchyIndex::count(std::pair<uint64_t, uint64_t> key) const
	{
		if (key.first < (size_t)std::numeric_limits<uint32_t>::max() && key.second < (size_t)std::numeric_limits<uint32_t>::max())
		{
			std::pair<uint32_t, uint32_t> smallkey { key.first, key.second };
			if (smalls.count(smallkey) == 1) return 1;
		}
		return bigs.count(key);
	}
	size_t HierarchyIndex::at(std::pair<uint64_t, uint64_t> key) const
	{
		if (key.first < (size_t)std::numeric_limits<uint32_t>::max() && key.second < (size_t)std::numeric_limits<uint32_t>::max())
		{
			std::pair<uint32_t, uint32_t> smallkey { key.first, key.second };
			if (smalls.count(smallkey) == 1) return smalls.at(key);
		}
		return bigs.at(key);
	}
	void HierarchyIndex::set(std::pair<uint64_t, uint64_t> key, size_t value)
	{
		if (key.first < (size_t)std::numeric_limits<uint32_t>::max() && key.second < (size_t)std::numeric_limits<uint32_t>::max() && value < (size_t)std::numeric_limits<uint32_t>::max())
		{
			std::pair<uint32_t, uint32_t> smallkey { key.first, key.second };
			assert(smalls.count(smallkey) == 0);
			smalls[smallkey] = value;
			return;
		}
		assert(bigs.count(key) == 0);
		bigs[key] = value;
	}
	size_t HierarchyIndex::size() const
	{
		return smalls.size() + bigs.size();
	}
}