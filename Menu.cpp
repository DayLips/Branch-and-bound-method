#include "Menu.h"

Menu::Menu(string filename)
{
	ifstream fin(filename);
	while (fin) menu.push_back(Bludo(fin));
	menu.pop_back();
	fin.close();
}

void Menu::Write()
{
	for (int i = 0; i < menu.size(); i++) menu[i].Write();
}

vector<Bludo> Menu::get_menu()
{
	return menu;
}

Menu::~Menu()
{
}