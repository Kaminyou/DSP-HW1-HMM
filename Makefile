.PHONY: all clean run
CC=g++
CFLAGS=-std=c++11
TARGET=train test
TRAIN_ITER=100

all: $(TARGET)

train: src/train.cpp
	$(CC) -o $@ $^ $(CFLAGS) -Iinc

test: src/test.cpp
	$(CC) -o $@ $^ $(CFLAGS) -Iinc

clean:
	rm -f $(TARGET)

