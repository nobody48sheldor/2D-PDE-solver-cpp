#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include<string>
// #include "/usr/include/python3.9/Python.h"
// #include "/usr/include/python3.9/matplotlibcpp.h"
// #include "matplotlibcpp.h"

using namespace std;

auto linspace(double l_0, double l_1, int n)
{
    vector<double> array;
    array.resize(n);

    for (int i = 0; i < n; i++)
    {
        array[i] = l_0 + (i*(l_1 - l_0)/n);
    }
    return array;
}

auto tau_init(vector<double> x, vector<double> y, double l, double a, int n)
{
    vector<vector<double>> yarray;
    yarray.resize(n);
    vector<double> xarray;
    xarray.resize(n);

    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            double value = 10*exp(-(pow((x[i]-(l/4)), 2) + (pow((y[j]-(l/4)), 2)))/a) / (sqrt(2)*3.14159*a);
            xarray[i] = value;
        }
        yarray[j] = xarray;
    }
    return yarray;
}

int output_data(vector<vector<vector<double>>> data,int nt,int n)
{
    ofstream output_file;
    output_file.open("data.txt");
    for (int timeT=0; timeT < 120; timeT++)
    {
        for (int ypos=0; ypos < n; ypos++)
        {
            //cout << ypos << endl;
            for (int xpos=0; xpos < n; xpos++)
            {
                output_file << data[(int)round(timeT*nt/120)][ypos][xpos] << endl;
                //cout << data[(int)round(timeT*nt/120)][ypos][xpos] << endl;
            }
        }
        output_file << endl;
    }
    output_file.close();
    return 0;
}

int write_data(vector<vector<double>> data, int nt, int n, int k)
{
    ofstream output_file;
    output_file.open("data.txt");
    for (int ypos=0; ypos < n; ypos++)
    {
        //cout << ypos << endl;
        for (int xpos=0; xpos < n; xpos++)
        {
            output_file << data[ypos][xpos] << '/';
            //cout << data[ypos][xpos] << ";";
            //cout << data[(int)round(timeT*nt/120)][ypos][xpos] << endl;
        }
    }
    output_file.close();
    return 0;
}

auto tau_evolution(vector<vector<double>> tau_0, double d, int treshold, double dx, double dy, double dt, int n, int nt)
{
    vector<vector<vector<double>>> tau;
    tau.resize(nt);
    vector<vector<double>> yarray;
    yarray.resize(n);
    vector<double> xarray;
    xarray.resize(n);

    tau[0] = tau_0;

    vector<vector<double>> empty;
    empty.resize(1);
    vector<double> empty_;
    empty_.resize(1);

    empty_[0] = 0;
    empty[0] = empty_;

    for(int k = 1; k < nt; k++)
    {
        if (k > treshold)
        {
            tau[k - treshold + 1] = empty;
        }
        if (100*k%nt == 0)
        {
            cout << 100*k/(nt) << " %" << endl;
        }
        yarray[0] = tau[k-1][0];
        for(int j = 1; j < n-1; j++)
        {
            xarray[0] = tau[k-1][j][0];
            for(int i = 1; i < n-1; i++)
            {
                xarray[i] = tau[k-1][j][i] + ((d * (((tau[k-1][j][i+1] - 2*tau[k-1][j][i] + tau[k-1][j][i-1]) / (dx*dx)) + ((tau[k-1][j+1][i] - 2*tau[k-1][j][i] + tau[k-1][j-1][i]) / (dy*dy)) ) ) *dt);
            }
            xarray[n-1] = xarray[n-2];
            yarray[j] = xarray;
        }
        yarray[n-1] = xarray;
        tau[k] = yarray;
        if (k%((int)round(nt/240)) == 0)
        {
            write_data(yarray, nt, n, (int)round(k/240));
            system("python3 plot-a-frame.py");
        }
    }
    return tau;
}

int main()
{
    int n = 200;
    int nt = 80000;
    int treshold = 40;
    double tmax = 60;
    double l = 5;
    double a = 2;
    double d = 0.1;
    vector<double> x = linspace(-l, l, n);
    vector<double> y = linspace(-l, l, n);
    vector<double> t = linspace(0, tmax, nt);

    double dx = x[1] - x[0];
    double dy = y[1] - y[0];
    double dt = t[1] - t[0];
    cout << dx << "  " << dy << "  " << dt << endl;

    vector<vector<double>> tau_0 = tau_init(x, y, l, a, n);
    vector<vector<vector<double>>> tau = tau_evolution(tau_0, d, treshold, dx, dy, dt, n, nt);
    cout << "done" << endl;
    cout << "now writing data into a file" << endl;
    return 0;
}
