# build library
ipTest: main.o
	rm -f ipTest
	g++ -o3 -o ipTest  main.o -L../Lib -lip
#	g++ -static -o ipTest  main.o -L../../lib -lip -lma -lcompare -lfingerSelector  -lwizard73
#	g++ -o ipTest  main.o -L../../lib -lip 

main.o: main.cpp
	g++ $(CFLAGS) -I../include -c main.cpp

# remove old files
clean:
	rm -f *.o
	rm -f ipTest