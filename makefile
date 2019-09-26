all: gsgd

gsgd: .obj/main.o .obj/Graph.o .obj/Utility.o
	g++ .obj/main.o .obj/Graph.o .obj/Utility.o -o gsgd -Wno-deprecated 

.obj/main.o: main.cpp
	g++ -c -O3 -o .obj/main.o main.cpp -Wno-deprecated 

.obj/Graph.o: Graph.cpp
	g++ -c -O3 -o .obj/Graph.o Graph.cpp -Wno-deprecated 

.obj/Utility.o: Utility.cpp
	g++ -c -O3 -o .obj/Utility.o Utility.cpp -Wno-deprecated 

clean:
	rm -rf *o .obj/
	mkdir .obj
