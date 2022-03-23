#ifndef READ_FILE_H
#define READ_FILE_H

#include <fstream>
#include <iostream>
#include <stdio.h>
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
	int NstepsPre = 10;
	int Nsteps = 10;
	int NFmode = 4;
	double nu = 0.5;
	double bendingdefmag = 0.0;
	double inplanedefmag = 0.0;
	char FolderName[MAXLINE];
	char FileName[MAXLINE];
	double bending00 = 0.0;
	double bending01 = 0.0;
	double bending11 = 0.0;
	double inplane00 = 0.0;
	double inplane01 = 0.0;
	double inplane11 = 0.0;

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
