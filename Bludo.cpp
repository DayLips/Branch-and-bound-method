#include "Bludo.h"

Bludo::Bludo(istream& fin)
{
    string line;
    getline(fin, line);
    istringstream fi(line);
    double R;
    fi >> this->Name;
    while (fi >> R) {
        param.push_back(R);
    }
}

void Bludo::Write()
{
    cout << Name << " ";
    for (int i = 0; i < param.size(); i++)
    {
        cout << param[i] << " ";
    }
    cout << endl;
}

vector<double> Bludo::get_param()
{
    return param;
}

char* Bludo::get_name()
{
    return Name;
}

Bludo::~Bludo()
{
}