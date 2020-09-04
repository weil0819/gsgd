all: gsgd

gsgd: .obj/main.o .obj/Graph.o .obj/Utility.o
	g++ .obj/main.o .obj/Graph.o .obj/Utility.o -o gsgd -std=c++11 -Wno-deprecated -Wno-unused-result 

.obj/main.o: main.cpp
	g++ -c -O3 -o .obj/main.o main.cpp -std=c++11 -Wno-deprecated -Wno-unused-result 

.obj/Graph.o: Graph.cpp
	g++ -c -O3 -o .obj/Graph.o Graph.cpp -std=c++11 -Wno-deprecated -Wno-unused-result 

.obj/Utility.o: Utility.cpp
	g++ -c -O3 -o .obj/Utility.o Utility.cpp -std=c++11 -Wno-deprecated -Wno-unused-result 

clean:
	rm -rf *o .obj/
	mkdir .obj
