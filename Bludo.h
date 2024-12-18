#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class Bludo
{
public:
	Bludo(istream& fin);
	void Write();
	vector<double> get_param();
	char* get_name();
	~Bludo();

private:
	char* Name = new char[50];
	vector<double> param;
};