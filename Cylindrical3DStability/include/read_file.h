#ifndef READ_FILE_H
#define READ_FILE_H

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <Eigen/Dense>
using namespace std;
#define MAXLINE 1024

class some_data {

public:
	some_data() {
	}
	;
	~some_data() {
	}
	;

	double h = 0.01;
	int refinelevel = 5;
	int Nsteps = 10;
	double nu = 0.5;
	char FolderName[MAXLINE];
	char FileName[MAXLINE];

};

class read_file {
public:
	read_file() {
	}
	;
	~read_file() {
	}
	;

	void readInputFile(char *filename, some_data &dat);

private:
	void getNextDataLine(FILE *const filePtr, char *nextLinePtr,
			int const maxSize, int *const endOfFileFlag);
};

#endif
