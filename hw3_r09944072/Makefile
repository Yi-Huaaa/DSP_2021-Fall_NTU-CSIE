SRIPATH ?= /root/srilm-1.5.10
MACHINE_TYPE ?= i686-m64
INC_PATH ?= inc
SRC_PATH ?= src
LM=bigram.lm

CXX=g++
CXXFLAGS=-O3 -I$(SRIPATH)/include -w --std=c++11

vpath lib%.a $(SRIPATH)/lib/$(MACHINE_TYPE)
vpath %.c $(SRC_PATH)
vpath %.cpp $(SRC_PATH)

TARGET = mydisambig
SRC = mydisambig.cpp
OBJ = $(SRC:.cpp=.o)
FROM ?= Big5-ZhuYin.map
TO ?= ZhuYin-Big5.map
.PHONY: all clean map

all: $(TARGET)

$(TARGET): $(OBJ) -loolm -ldstruct -lmisc
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

map:
	@python3 src/mapping.py $(FROM) $(TO)


run:
	./mydisambig test_data/1_seg.txt ZhuYin-Big5.map bigram.lm result2/1.txt
	./mydisambig test_data/2_seg.txt ZhuYin-Big5.map bigram.lm result2/2.txt
	./mydisambig test_data/3_seg.txt ZhuYin-Big5.map bigram.lm result2/3.txt
	./mydisambig test_data/4_seg.txt ZhuYin-Big5.map bigram.lm result2/4.txt
	./mydisambig test_data/5_seg.txt ZhuYin-Big5.map bigram.lm result2/5.txt
	./mydisambig test_data/6_seg.txt ZhuYin-Big5.map bigram.lm result2/6.txt
	./mydisambig test_data/7_seg.txt ZhuYin-Big5.map bigram.lm result2/7.txt
	./mydisambig test_data/8_seg.txt ZhuYin-Big5.map bigram.lm result2/8.txt
	./mydisambig test_data/9_seg.txt ZhuYin-Big5.map bigram.lm result2/9.txt
	./mydisambig test_data/10_seg.txt ZhuYin-Big5.map bigram.lm result2/10.txt


clean:
	$(RM) $(OBJ) $(TARGET)