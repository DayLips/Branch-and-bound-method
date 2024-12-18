#pragma once
#include <iostream>
#include <vector>
#include "Bludo.h"

using namespace std;

class Menu
{
public:
	Menu(string filename);
	void Write();
	vector<Bludo> get_menu();
	~Menu();

private:
	vector<Bludo> menu;
};