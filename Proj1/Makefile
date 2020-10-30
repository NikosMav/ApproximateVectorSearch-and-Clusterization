.PHONY: all

all: lsh cluster cube

lsh: lib/Aux.o lib/Hashtable.o lib/Image.o lib/LSH.o lib/Bucket.o lib/HyperCube.o lib/distFunc.o source/LSHmain.cpp 
	g++ -o lsh source/LSHmain.cpp lib/Aux.o lib/Image.o lib/LSH.o lib/Hashtable.o lib/Bucket.o lib/HyperCube.o lib/distFunc.o -ggdb3 -std=c++11 

cluster: lib/Aux.o lib/Hashtable.o lib/Image.o lib/LSH.o lib/Bucket.o lib/HyperCube.o lib/distFunc.o source/ClusteringMain.cpp 
	g++ -o cluster source/ClusteringMain.cpp lib/Aux.o lib/Image.o lib/LSH.o lib/Hashtable.o lib/Bucket.o lib/HyperCube.o lib/distFunc.o -ggdb3 -std=c++11 

cube: lib/Aux.o lib/Hashtable.o lib/Image.o lib/LSH.o lib/Bucket.o lib/HyperCube.o lib/distFunc.o source/Cubemain.cpp 
	g++ -o cube source/Cubemain.cpp lib/Aux.o lib/Image.o lib/LSH.o lib/Hashtable.o lib/Bucket.o lib/HyperCube.o lib/distFunc.o -ggdb3 -std=c++11 

lib/Aux.o: source/Aux.cpp headers/Aux.h headers/Hashtable.h
	g++ -c source/Aux.cpp -ggdb3 -std=c++11
	mv Aux.o lib/Aux.o

lib/Hashtable.o: source/Hashtable.cpp headers/Hashtable.h
	g++ -c source/Hashtable.cpp -ggdb3 -std=c++11
	mv Hashtable.o lib/Hashtable.o

lib/Image.o: source/Image.cpp headers/Image.h headers/Hashtable.h
	g++ -c source/Image.cpp -ggdb3 -std=c++11
	mv Image.o lib/Image.o

lib/LSH.o: source/LSH.cpp headers/LSH.h headers/Hashtable.h
	g++ -c source/LSH.cpp -ggdb3 -std=c++11
	mv LSH.o lib/LSH.o

lib/Bucket.o: source/Bucket.cpp headers/Bucket.h headers/Hashtable.h
	g++ -c source/Bucket.cpp -ggdb3 -std=c++11
	mv Bucket.o lib/Bucket.o


lib/HyperCube.o: source/HyperCube.cpp headers/HyperCube.h headers/Hashtable.h
	g++ -c source/HyperCube.cpp -ggdb3 -std=c++11
	mv HyperCube.o lib/HyperCube.o

lib/distFunc.o: source/distFunc.cpp headers/distFunc.h
	g++ -c source/distFunc.cpp -ggdb3 -std=c++11
	mv distFunc.o lib/distFunc.o

clean: 
	rm lsh cube cluster lib/*