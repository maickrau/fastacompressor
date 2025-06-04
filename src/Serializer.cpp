#include "Serializer.h"

namespace Serializer
{
	void writeUint8_t(uint8_t value, std::ostream& stream)
	{
		stream.write((char*)&value, 1);
	}
	uint8_t readUint8_t(std::istream& stream)
	{
		uint8_t c = 0;
		stream.read((char*)&c, 1);
		return c;
	}
	void writeUint16_t(uint16_t value, std::ostream& stream)
	{
		// endianness
		for (size_t i = 0; i < 2; i++)
		{
			writeUint8_t(value % 256, stream);
			value /= 256;
		}
	}
	uint16_t readUint16_t(std::istream& stream)
	{
		uint16_t result = 0;
		for (size_t i = 0; i < 2; i++)
		{
			result >>= 8;
			result += readUint8_t(stream) << 8;
		}
		return result;
	}
	void writeUint64_t(uint64_t value, std::ostream& stream)
	{
		// endianness
		for (size_t i = 0; i < 8; i++)
		{
			writeUint8_t(value % 256, stream);
			value /= 256;
		}
	}
	uint64_t readUint64_t(std::istream& stream)
	{
		uint64_t result = 0;
		for (size_t i = 0; i < 8; i++)
		{
			result >>= 8;
			result += ((uint64_t)readUint8_t(stream)) << 56ull;
		}
		return result;
	}
	void writeString(const std::string& str, std::ostream& stream)
	{
		writeUint64_t(str.size(), stream);
		for (size_t i = 0; i < str.size(); i++)
		{
			writeUint8_t(str[i], stream);
		}
	}
	std::string readString(std::istream& stream)
	{
		size_t size = readUint64_t(stream);
		std::string result;
		result.resize(size);
		for (size_t i = 0; i < result.size(); i++)
		{
			result[i] = readUint8_t(stream);
		}
		return result;
	}
}
