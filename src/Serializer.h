#include <cstdint>
#include <vector>
#include <string>
#include <fstream>

namespace Serializer
{
	void writeUint8_t(uint8_t value, std::ostream& stream);
	uint8_t readUint8_t(std::istream& stream);
	void writeUint16_t(uint16_t value, std::ostream& stream);
	uint16_t readUint16_t(std::istream& stream);
	void writeUint64_t(uint64_t value, std::ostream& stream);
	uint64_t readUint64_t(std::istream& stream);
	template <typename T>
	void write(const T value, std::ostream& stream)
	{
		// weird workaround, force SFINAE to use specialized template functions. if specialized function missing, causes compilation failure
		// https://stackoverflow.com/questions/76375725/why-cant-i-use-static-assert-for-unimplemented-and-unused-template-function
		static_assert(sizeof(T) == -1);
	}
	template <typename T>
	T read(std::istream& stream)
	{
		// weird workaround, force SFINAE to use specialized template functions. if specialized function missing, causes compilation failure
		// https://stackoverflow.com/questions/76375725/why-cant-i-use-static-assert-for-unimplemented-and-unused-template-function
		static_assert(sizeof(T) == -1);
	}
	template<>
	inline void write(const uint8_t value, std::ostream& stream)
	{
		writeUint8_t(value, stream);
	}
	template<>
	inline void write(const uint16_t value, std::ostream& stream)
	{
		writeUint16_t(value, stream);
	}
	template<>
	inline void write(const uint64_t value, std::ostream& stream)
	{
		writeUint64_t(value, stream);
	}
	template<>
	inline uint8_t read(std::istream& stream)
	{
		return readUint8_t(stream);
	}
	template<>
	inline uint16_t read(std::istream& stream)
	{
		return readUint16_t(stream);
	}
	template<>
	inline uint64_t read(std::istream& stream)
	{
		return readUint64_t(stream);
	}
	void writeString(const std::string& str, std::ostream& stream);
	std::string readString(std::istream& stream);
	template <typename T>
	void writeVector(const std::vector<T>& vector, std::ostream& stream)
	{
		writeUint64_t(vector.size(), stream);
		for (size_t i = 0; i < vector.size(); i++)
		{
			write<T>(vector[i], stream);
		}
	}
	template <typename T>
	std::vector<T> readVector(std::istream& stream)
	{
		size_t size = readUint64_t(stream);
		std::vector<T> result;
		result.reserve(size);
		for (size_t i = 0; i < size; i++)
		{
			result.emplace_back(read<T>(stream));
		}
		return result;
	}
}
