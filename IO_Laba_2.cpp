#include <iostream>
#include <vector>
#include <iomanip>
#include "Bludo.h"
#include "Menu.h"
//#include <glpk.h>

using namespace std;

const double M = 1000000.0;                        // Штраф

bool Proverka_optimum(vector<double> t);

struct M_metod
{
    vector<vector<double>> matr;
    vector<vector<int>> basis;
    bool succes;
    M_metod() : succes(false) {}
};

M_metod method(Menu data, int max_min, vector<double> limit)
{
    M_metod res;
    // Создание симплекс таблицы
    int n = 2 * data.get_menu()[0].get_param().size() - 1 + data.get_menu().size();
    int m = data.get_menu().size() + (n - 1) + data.get_menu()[0].get_param().size();
    vector<vector<double>> matr(n, vector<double>(m));          // Симплекс таблица
    vector<vector<int>> basis(n - 1, vector<int>(2));           // Вектор для отслеживания базисных векторов и искуств.перем.
    int limit_products = 1;
    // Заполнение последнего столбца в симплекс таблице до ограничений на кол-во продуктов
    for (int i = 0; i < (2 * data.get_menu()[0].get_param().size() - 2); i += 2) 
    {
        matr[i + 1][m - 1] = limit[i + 1];
        matr[i][m - 1] = limit[i];
    }
    // Заполнение последнего столбца в симплекс таблице а именно ограничения на кол-во продуктов
    for (int i = (2 * data.get_menu()[0].get_param().size() - 2); i < (n - 1); i++) matr[i][m - 1] = limit_products;
    // Заполняем вектор basis
    for (int i = 0, j = data.get_menu().size(), k = n - 1; i < basis.size(); i++) 
    {
        if ((i % 2 == 0) && (j < (n - 1)))
        {
            basis[i][0] = j;
            basis[i][1] = 0;
            j += 2;
        }
        else
        {
            basis[i][0] = k;
            basis[i][1] = 1;
            k++;
        }
    }

    for (int j = 0; j < (m - 1); j++)
    {
        for (int i = 0; i < n; i++)
        {
            if (i == (n - 1))                   // Заполнение строки z
            {
                if (j < data.get_menu().size())
                {
                    matr[i][j] = -max_min * data.get_menu()[j].get_param()[data.get_menu()[0].get_param().size() - 1];
                }
            }
            if (i < (n - 1))                    // Заполнение векторов блюд, искуственных перемнных и дополнительных
            {
                if (j < data.get_menu().size())
                {
                    if (((i % 2) == 0) && (i <= data.get_menu()[0].get_param().size()))
                    {
                        matr[i][j] = data.get_menu()[j].get_param()[i / 2];
                        matr[i + 1][j] = data.get_menu()[j].get_param()[i / 2];
                    }
                }
                else
                {
                    if (basis[i][0] == j)
                    {
                        matr[i][j] = 1;
                    }
                    else if (((j - data.get_menu().size()) == i) && (j < (n - 1)))
                    {
                        matr[i - 1][j] = -1;
                    }
                }
            }

        }
    }
    for (int j = 0, i = 2*data.get_menu()[0].get_param().size() - 2; j < data.get_menu().size(); j++)
    {
        matr[i][j] = 1;
        i++;
    }
    for (int i = 0; i < basis.size(); i++)
    {
        if (basis[i][1] == 0)
        {
            matr[n - 1][basis[i][0]] = -M;
        }
    }
    for (int i = 0; i < n - 1; i++)     // Приводим симплекс таблицу к применению метода(зануляем коэффициенты у искуст. перем.)
    {
        for (int j = 0; j < m; j++)
        {
            if (basis[i][1] == 0)
            {
                matr[n - 1][j] += M * matr[i][j];
            }

        }
    }
    bool flag = false;
    while (!Proverka_optimum(matr[n - 1]))
    {
        flag = false;
        vector<int> desision_el = { 0,0 };                              // Находим разрешающий элемент
        double drob = -1;
        for (int j = 0; j < (m - 1); j++)
        {
            if ((matr[n - 1][j] > 0) && (matr[n - 1][j] > drob))
            {
                drob = matr[n - 1][j];
                desision_el[1] = j;
            }
        }
        drob = matr[0][m - 1] / matr[0][desision_el[1]];
        for (int i = 1; i < (n - 1); i++)
        {
            if (matr[i][desision_el[1]] > 0)
            {
                flag = true;
                if ((matr[i][m - 1] / matr[i][desision_el[1]]) < drob)
                {
                    drob = matr[i][m - 1] / matr[i][desision_el[1]];
                    desision_el[0] = i;
                }
            }
        }
        if (!flag)
        {
            if (matr[0][desision_el[1]] <= 0)
            {
                return res;
                break;
            }
        }
        if (basis[desision_el[0]][1] == 1)                                  // Удаляем столбец с искуст. перем.
        {
            for (auto& row : matr)
            {
                if (row.size() > (basis[desision_el[0]][0]))
                {
                    row.erase(row.begin() + basis[desision_el[0]][0]);
                }
            }
            m -= 1;

            for (int i = 0; i < basis.size(); i++)
            {
                if (basis[i][0] > basis[desision_el[0]][0])
                {
                    basis[i][0] -= 1;
                }
            }
            basis[desision_el[0]][0] = desision_el[1];
            basis[desision_el[0]][1] = 0;
        }
        else
        {
            basis[desision_el[0]][0] = desision_el[1];
            basis[desision_el[0]][1] = 0;
        }

        double desision = matr[desision_el[0]][desision_el[1]];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                if ((i != desision_el[0]) && (j != desision_el[1]))
                {
                    matr[i][j] = (matr[i][j] * desision - matr[i][desision_el[1]] * matr[desision_el[0]][j]) / desision;
                }

            }
        }
        for (int i = 0; i < n; i++)
        {
            if (i != desision_el[0]) matr[i][desision_el[1]] = 0;
        }
        for (int j = 0; j < m; j++)
        {
            if (j != desision_el[1]) matr[desision_el[0]][j] = matr[desision_el[0]][j] / desision;
        }
        matr[desision_el[0]][desision_el[1]] = 1;
    }

    res.succes = true;
    res.matr = matr;
    res.basis = basis;
    return res;
}

void Print_res(M_metod res, Menu data);



int main()
{
    setlocale(LC_ALL, "RUS");
    Menu data("Menu.txt");
    int max_min = 0;
    double delta = 0;
    vector<double> limit(2 * data.get_menu()[0].get_param().size() - 2);
    cout << "Если хотите найти максимум или минимум, выберите -1 или 1 соотвественно: ";
    cin >> max_min;
    cout << "Введите максимальное отклонение ограничений в %: ";
    cin >> delta;
    for (int i = 0; i < limit.size(); i += 2)
    {
        cout << "Введите значение ограничения:";
        cin >> limit[i];
        limit[i + 1] = round((1.0 - delta / 100.0) * limit[i]);
        limit[i] = round((1.0 + delta / 100.0) * limit[i]);
    }
    cout << endl;
    M_metod table = method(data, max_min, limit);
    Print_res(table, data);
}

//int main() {
//    int a = 10;
//    int b = 99;
//
//}

bool Proverka_optimum(vector<double> t)
{
    bool flag = true;
    for (int i = 0; i < (t.size() - 1); i++)
    {
        if ((t[i] / 1.0) > 0)
        {
            flag = false;
            break;
        }
    }
    return flag;
}

void Print_res(M_metod res, Menu data)
{
    int n = res.matr.size();
    int m = res.matr[1].size();
    if (!res.succes) cout << "Нет решения!!!" << endl;
    else
    {
        double summ = 0;
        for (int i = 0; i < res.basis.size(); i++)
        {
            if (res.basis[i][0] < data.get_menu().size())
            {
                double x = res.matr[i][m - 1];
                if (i % 2 == 0)
                {
                    cout << fixed << setprecision(2) << "Количество " << data.get_menu()[res.basis[i][0]].get_name() << " = " << x << ", парметры блюда умноженное на кол-во =";
                    for (int j = 0; j < data.get_menu()[res.basis[i][0]].get_param().size(); j++)
                    {
                        cout << " " << fixed << setprecision(2) << x * data.get_menu()[res.basis[i][0]].get_param()[j] << ",";
                        if (j == data.get_menu()[res.basis[i][0]].get_param().size() - 1) summ += x * data.get_menu()[res.basis[i][0]].get_param()[j];
                    }
                }
                else
                {
                    cout << fixed << setw(6) << setprecision(2) << "Количество " << data.get_menu()[res.basis[i][0]].get_name() << " = " << x << ", парметры блюда умноженное на кол-во =";
                    for (int j = 0; j < data.get_menu()[res.basis[i][0]].get_param().size(); j++)
                    {
                        cout << " " << fixed << setprecision(2) << x * data.get_menu()[res.basis[i][0]].get_param()[j] << ",";
                        if (j == data.get_menu()[res.basis[i][0]].get_param().size() - 1) summ += x * data.get_menu()[res.basis[i][0]].get_param()[j];
                    }
                }
                cout << endl;

            }
        }
        cout << endl << "Итоговая цена = " << summ << " рублей " << endl;
    }
}