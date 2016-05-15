DEP_DIR:=./deps
SRC_DIR:=src
BIN_DIR:=bin
OBJ_DIR:=obj
LIB_DIR:=lib
INC_DIR:=include
CPP_DIR:=cpp

EXE:=hhga

CXX:=g++
CXXFLAGS:=-O3 -msse4.1 -fopenmp -std=c++11 -ggdb

CWD:=$(shell pwd)

LD_INCLUDE_FLAGS:=-I$(CWD)/$(INC_DIR) -I. -I$(CWD)/$(SRC_DIR) -I$(CWD)/$(CPP_DIR) -I$(CWD)/$(INC_DIR)/bamtools
LD_LIB_FLAGS:= -ggdb -L$(CWD)/$(LIB_DIR) -lvg -lvcflib -lhts -lbamtools -lpthread -lz -lm -lbz2

RAPTOR_INCLUDE:=/usr/include/
ifeq ($(shell uname -s),Darwin)
    # We may need libraries from Macports
    # TODO: where does Homebrew keep libraries?
    ifeq ($(shell if [ -d /opt/local/lib ];then echo 1;else echo 0;fi), 1)
       LD_LIB_FLAGS += -L/opt/local/lib
       RAPTOR_INCLUDE=/opt/local/include/
    endif
    ifeq ($(shell if [ -d /usr/local/lib ];then echo 1;else echo 0;fi), 1)
       LD_LIB_FLAGS += -L/usr/local/lib
       RAPTOR_INCLUDE=/usr/local/include/
    endif
    ROCKSDB_PORTABLE=PORTABLE=1 # needed to build rocksdb without weird assembler options
else
    # Not on OS X, we can have librt
    LD_LIB_FLAGS += -lrt
endif

OBJ:=$(OBJ_DIR)/hhga.o

SDSL_DIR:=deps/sdsl-lite
FASTAHACK_DIR:=deps/fastahack
HTSLIB_DIR:=deps/htslib
VCFLIB_DIR:=deps/vcflib
BAMTOOLS_DIR=deps/bamtools
VGLIB_DIR:=deps/vg

STATIC_FLAGS=-static -static-libstdc++ -static-libgcc

.PHONY: clean get-deps test set-path static

$(BIN_DIR)/hhga: $(LIB_DIR)/libhhga.a $(OBJ_DIR)/main.o deps
	. ./source_me.sh && $(CXX) $(CXXFLAGS) -o $(BIN_DIR)/hhga $(OBJ_DIR)/main.o $(LD_INCLUDE_FLAGS) -lhhga $(LD_LIB_FLAGS)

static: $(LIB_DIR)/libhhga.a $(OBJ_DIR)/main.o $(OBJ) deps
	$(CXX) $(CXXFLAGS) -o $(BIN_DIR)/hhga $(OBJ_DIR)/main.o $(OBJ) $(STATIC_FLAGS) $(LD_INCLUDE_FLAGS) -lhhga $(LD_LIB_FLAGS)	

$(LIB_DIR)/libhhga.a: $(OBJ)
	ar rs $@ $^

#get-deps:
#	sudo apt-get install -qq -y protobuf-compiler libprotoc-dev libjansson-dev libbz2-dev libncurses5-dev automake libtool jq samtools curl unzip redland-utils librdf-dev cmake pkg-config wget bc

test: $(BIN_DIR)/hhga
	. ./source_me.sh && cd test && $(MAKE)

deps: $(LIB_DIR)/libvg.a $(LIB_DIR)/libvcflib.a $(LIB_DIR)/libhts.a $(LIB_DIR)/libbamtools.a $(OBJ_DIR)/Fasta.o

$(OBJ_DIR)/Fasta.o: .pre-build
	+cd $(FASTAHACK_DIR) && make && mv Fasta.o $(CWD)/$(OBJ_DIR) && cp Fasta.h $(CWD)/$(INC_DIR)

$(LIB_DIR)/libhts.a: .pre-build
	+cd $(HTSLIB_DIR) && $(MAKE) lib-static && mv libhts.a $(CWD)/$(LIB_DIR) && cp *.h $(CWD)/$(INC_DIR) && cp -r htslib $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libvcflib.a: .pre-build
	+. ./source_me.sh && cd $(VCFLIB_DIR) && $(MAKE) libvcflib.a && cp lib/* $(CWD)/$(LIB_DIR)/ && cp include/* $(CWD)/$(INC_DIR)/ && cp src/*.h* $(CWD)/$(INC_DIR)/

$(LIB_DIR)/libvg.a: .pre-build
	+. ./source_me.sh && cd $(VGLIB_DIR) && $(MAKE) && cp lib/libvg.a $(CWD)/$(LIB_DIR)/ && cp include/* $(CWD)/$(INC_DIR)/ && cp src/*.h* $(CWD)/$(INC_DIR)/

# builds bamtools static lib, and copies into root
$(LIB_DIR)/libbamtools.a:
	+cd $(BAMTOOLS_DIR) && mkdir -p build && cd build && cmake .. && make && cp ../lib/libbamtools.a $(CWD)/$(LIB_DIR) && cp -r ../src $(CWD)/$(INC_DIR)/bamtools


###################################
## HHGA source code compilation begins here
####################################

$(OBJ_DIR)/hhga.o: $(SRC_DIR)/hhga.cpp $(SRC_DIR)/hhga.hpp deps
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

$(OBJ_DIR)/main.o: $(SRC_DIR)/main.cpp $(SRC_DIR)/hhga.hpp deps
	+. ./source_me.sh && $(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_INCLUDE_FLAGS) $(LD_LIB_FLAGS)

.pre-build:
	if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi
	if [ ! -d $(LIB_DIR) ]; then mkdir -p $(LIB_DIR); fi
	if [ ! -d $(OBJ_DIR) ]; then mkdir -p $(OBJ_DIR); fi
	if [ ! -d $(INC_DIR) ]; then mkdir -p $(INC_DIR); fi
	if [ ! -d $(CPP_DIR) ]; then mkdir -p $(CPP_DIR); fi
	touch .pre-build

# for rebuilding just vg
clean-hhga:
	$(RM) -r $(BIN_DIR)/hhga
	$(RM) -r $(OBJ_DIR)/*
	$(RM) -r $(CPP_DIR)/*

clean:
	$(RM) -r $(BIN_DIR)
	$(RM) -r $(LIB_DIR)
	$(RM) -r $(OBJ_DIR)
	$(RM) -r $(INC_DIR)
	$(RM) -r $(CPP_DIR)
	$(RM) -r share/
	$(RM) -f .pre-build
	cd $(DEP_DIR) && cd vcflib && $(MAKE) clean
	cd $(DEP_DIR) && cd vg && $(MAKE) clean
	cd $(DEP_DIR) && cd htslib && $(MAKE) clean
	cd $(DEP_DIR) && cd fastahack && $(MAKE) clean
	cd $(DEP_DIR) && rm -rf bamtools/build
