all: rc clean

%.o: %.cpp
	g++ -Wall -c -o $@ $< -W -O 


rc: rc-Main.o rc.o rc-Burst.o rc-Par.o rc-Est.o nrutil.c
	g++ -o $@ -L'/usr/local/lib' $^ -W -lgsl -lgslcblas -lm -lz

clean:
	rm -rf *.o
