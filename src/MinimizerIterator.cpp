#include "MinimizerIterator.h"

namespace MinimizerIterator
{
	size_t charToInt(char c)
	{
		switch(c)
		{
			case 'a':
			case 'A':
				return 0;
			case 'c':
			case 'C':
				return 1;
			case 'g':
			case 'G':
				return 2;
			case 't':
			case 'T':
				return 3;
		}
		assert(false);
		return 0;
	}

	std::vector<bool> getValidChars()
	{
		std::vector<bool> result;
		result.resize(256, false);
		result['a'] = true;
		result['A'] = true;
		result['c'] = true;
		result['C'] = true;
		result['g'] = true;
		result['G'] = true;
		result['t'] = true;
		result['T'] = true;
		return result;
	}

	// https://naml.us/post/inverse-of-a-hash-function/
	uint64_t hash(uint64_t key) {
		key = (~key) + (key << 21); // key = (key << 21) - key - 1;
		key = key ^ (key >> 24);
		key = (key + (key << 3)) + (key << 8); // key * 265
		key = key ^ (key >> 14);
		key = (key + (key << 2)) + (key << 4); // key * 21
		key = key ^ (key >> 28);
		key = key + (key << 31);
		return key;
	}

	std::vector<bool> validChar = getValidChars();

}
