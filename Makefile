all: datastat.x

CXX = g++
CXXFLAGS = -O3 -Wall

datastat.x: datastat.C
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f *.o 

