#include <string>
#include <cassert>
#include "CompressedStringIndex.h"

int main(int argc, char** argv)
{
	// these parameters work decently for ONT ultralong r10.4.1 reads
	FastaCompressor::CompressedStringIndex index { 5, 100 };
	size_t numThreads = 8;

	// add sequences: if ID given as a parameter, the sequence ID in the index will be the specific ID regardless of the order the reads are added
	size_t readID = index.addString(2, "read2", "AGCGGGCAGTCAGTCAGTAGTGTCAGTCAGTAGTC");
	assert(readID == 2);
	readID = index.addString(0, "read0", "ACCCGAGTCTCATACTCAACGTCAG");
	assert(readID == 0);
	readID = index.addString(1, "read1", "ACGAGTCTTTCATACTCAACGTCAGACGAGTCTCATACTCAACGTCAG");
	assert(readID == 1);
	// add sequences: if ID is NOT given as a parameter, the sequence ID in the index will be arbitrary
	readID = index.addString("readN", "ATCCCTCATCTCATCAGTACTCAGTTCATC");
	// all sequences must be added to the index before sequences can be read from the index
	// only ACTG characters are allowed in the sequences. Lowercase and uppercase are considered the same character
	// any characters are allowed in names
	// adding sequences is thread safe and multiple threads can add sequences at the same time without locking

	// after all reads have been added, freeze the index
	index.removeConstructionVariables(numThreads);
	// now the index is compressed and immutable

	// read names and sequences can be recovered
	assert(index.getSequence(2) == "AGCGGGCAGTCAGTCAGTAGTGTCAGTCAGTAGTC");
	assert(index.getName(2) == "read2");
	assert(index.getName(readID) == "readN");
	// accessing sequences and names is thread safe and multiple threads can read at the same time

	// substrings can be extracted slightly faster with a specialized function instead of extracting the entire read
	assert(index.getSubstring(2, 5, 10) == index.getSequence(2).substr(5, 10));
}
