
CXX = clang++
CXX_FLAGS=-Iinclude/ -I/usr/include/ -Wall
TARGET = htx_assessment

SRC_FILES = $(wildcard src/*.cpp)

$(TARGET): $(SRC_FILES)
	$(CXX) $(CXX_FLAGS) -o $(TARGET) $(SRC_FILES)

.PHONY: clean
clean:
	rm -f ./htx_assessment

all: $(TARGET)

run: $(TARGET)
	./$(TARGET)