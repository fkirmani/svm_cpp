CC=g++
CFLAGS=-c 
#CFLAGS=-std=c++11 -stdlib=libc++

all: model

model: main.o svm.o
	$(CC) -lm main.o -o model libsvm.so.2

main.o: main.cpp calculate_features.h svm.h
	$(CC) $(CFLAGS) main.cpp -Isvm.h

svm.o: svm.cpp svm.h
	$(CC) $(CFLAGS) svm.cpp

clean:
	rm *o model
