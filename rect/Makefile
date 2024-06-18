CC=g++
CFLAGS=-c -std=c++11 -O3 # -ggdb3 -MMD -Wall -Wextra -Wpedantic
LDFLAGS=-static
SOURCES=edge.cpp rectangle.cpp vertex.cpp wall.cpp rectangulation.cpp main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=rect
OTHERS=edge.hpp rectangle.hpp rectangulation.hpp vertex.hpp wall.hpp Makefile
ZIP=zip
ZIPFILE=rect.zip

all: $(EXECUTABLE) $(ZIPFILE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

%.o : %.cpp
	$(CC) $(CFLAGS) $< -o $@

$(ZIPFILE): $(SOURCES) $(OTHERS)
	$(ZIP) $(ZIPFILE) $(SOURCES) $(OTHERS)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE) $(ZIPFILE)

.PHONY: all clean
