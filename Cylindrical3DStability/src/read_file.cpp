#ifndef READ_FILE_CC
#define READ_FILE_CC
#include "read_file.h"

#include <iostream>
#include <string.h>

void read_file::readInputFile(char *filename, some_data &dat) {
	FILE *fid;
	int endOfFileFlag;
	char nextLine[MAXLINE];

	int valuesWritten;
	bool fileReadErrorFlag = false;

	fid = std::fopen(filename, "r");
	if (fid == NULL) {
		std::cout << "Unable to open file \"" << filename << "\"" << std::endl;
		fileReadErrorFlag = true;
	} else {
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%lg", &(dat.h));
		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// read in an Nn
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%u", &(dat.refinelevel));

		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}

		// read in Nstep
		getNextDataLine(fid, nextLine, MAXLINE, &endOfFileFlag);
		valuesWritten = sscanf(nextLine, "%u", &(dat.Nsteps));
		if (valuesWritten != 1) {
			fileReadErrorFlag = true;
			goto fileClose;
		}




		fileClose: {
			fclose(fid);
		}
	}

	if (fileReadErrorFlag) {
		// default parameter values
		std::cout << "Error reading input file, Exiting.\n" << std::endl;
		exit(1);
	} else
		std::cout << "Input file successfully read" << std::endl;

}

void read_file::getNextDataLine(FILE *const filePtr, char *nextLinePtr,
		int const maxSize, int *const endOfFileFlag) {
	*endOfFileFlag = 0;
	do {
		if (fgets(nextLinePtr, maxSize, filePtr) == NULL) {
			*endOfFileFlag = 1;
			break;
		}
		while ((nextLinePtr[0] == ' ' || nextLinePtr[0] == '\t')
				|| (nextLinePtr[0] == '\n' || nextLinePtr[0] == '\r')) {
			nextLinePtr = (nextLinePtr + 1);
		}
	} while ((strncmp("#", nextLinePtr, 1) == 0) || (strlen(nextLinePtr) == 0));
}

#endif
