
CXX = clang++
CXX_FLAGS=-Iinclude/
TARGET = htx_assessment

SRC_FILES = $(wildcard src/*.cpp)

$(TARGET): $(SRC_FILES)
	$(CXX) $(CXX_FLAGS) -o $(TARGET) $(SRC_FILES)

all: $(TARGET)

run: $(TARGET)
	./$(TARGET)