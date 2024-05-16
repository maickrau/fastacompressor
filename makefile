GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -g -Izstr/src -Iparallel-hashmap/parallel_hashmap/ -Wno-unused-parameter `pkg-config --cflags zlib`

ODIR=obj
BINDIR=bin
SRCDIR=src

LIBS=`pkg-config --libs zlib`

_DEPS = fastqloader.h MinimizerIterator.h CommonUtils.h CompressedIndex.h CompressedStringIndex.h VariableWidthIntVector.h StringContainer.h RankBitvector.h StringHashIndex.h HierarchyIndex.h
DEPS = $(patsubst %, $(SRCDIR)/%, $(_DEPS))

_OBJ = fastqloader.o MinimizerIterator.o CommonUtils.o CompressedIndex.o CompressedStringIndex.o VariableWidthIntVector.o StringContainer.o RankBitvector.o StringHashIndex.o HierarchyIndex.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++

VERSION := Branch $(shell git rev-parse --abbrev-ref HEAD) commit $(shell git rev-parse HEAD) $(shell git show -s --format=%ci)

$(shell mkdir -p bin)
$(shell mkdir -p obj)

$(BINDIR)/fastacompress: $(OBJ) $(ODIR)/main.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/fastacompress_multithread: $(OBJ) $(ODIR)/main_multithread.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

all: $(BINDIR)/fastacompress $(BINDIR)/fastacompress_multithread

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
