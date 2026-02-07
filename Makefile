# Makefile for kallisto C implementation

CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99
LDFLAGS = -lm -lz -lpthread

TARGET = kallisto
BINARY_CONVERT = binary_convert
SOURCES = kallisto.c kthread.c
CONVERT_SOURCES = Binary_convert.c
OBJECTS = $(SOURCES:.c=.o)
CONVERT_OBJECTS = $(CONVERT_SOURCES:.c=.o)

.PHONY: all clean

all: $(TARGET) $(BINARY_CONVERT)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BINARY_CONVERT): $(CONVERT_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(TARGET) $(BINARY_CONVERT) $(OBJECTS) $(CONVERT_OBJECTS)

debug: CFLAGS += -g -DDEBUG
debug: clean all
