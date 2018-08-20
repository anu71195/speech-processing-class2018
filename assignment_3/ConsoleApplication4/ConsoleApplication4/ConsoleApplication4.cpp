// ConsoleApplication3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
int _tmain(int argc, _TCHAR* argv[])
{
	ifstream inFile;
	string end_line_text,item;
	int  count = 0, bitspersample = 16, channels = 1, samplerate = 16000, normalized = 0;
	inFile.open("data/ambient_sound.txt");
	if (inFile.fail())
	{
		cerr << "Error opening file" << endl;
		exit(1);
	}
	//while (!inFile.eof())
	//{
	//	inFile >> item;
	//	cout << item<<endl;
	//	count++;
	//}
	//cout << "number of items are " << count;
	//cin >> end_line_text;

	return 0;
}
