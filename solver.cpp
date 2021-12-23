#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include<string>

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

auto tau_init_gaussian(vector<double> x, vector<double> y, double l, double a, int n)
{
    vector<vector<double>> yarray;
    yarray.resize(n);
    vector<double> xarray;
    xarray.resize(n);

    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            double value = 10*exp(- a*(pow((x[i]-(l/4)), 2) + (pow((y[j]-(l/4)), 2)))) / (sqrt(2)*3.14159*a);
            xarray[i] = value;
        }
        yarray[j] = xarray;
    }
    return yarray;
}

auto tau_init_sin_gaussian(vector<double> x, vector<double> y, double l, double a, int n)
{
    vector<vector<double>> yarray;
    yarray.resize(n);
    vector<double> xarray;
    xarray.resize(n);

    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            double value = 1.2*(sin(x[i]*3.14159*a/l) * sin(y[j]*3.14159*a/l)) * exp(- 0.1*(pow((x[i]-(l/4)), 2) + (pow((y[j]-(l/4)), 2)))/a) / (sqrt(2)*3.14159*a);
            xarray[i] = value;
        }
        yarray[j] = xarray;
    }
    return yarray;
}

auto tau_init_sin(vector<double> x, vector<double> y, double l, double a, int n)
{
    vector<vector<double>> yarray;
    yarray.resize(n);
    vector<double> xarray;
    xarray.resize(n);

    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            double value = 0.6*(sin(x[i]*3.14159*a/l) * sin(y[j]*3.14159*a/l));
            xarray[i] = value;
        }
        yarray[j] = xarray;
    }
    return yarray;
}

auto tau_init_q(vector<double> x, vector<double> y, double l, double a, int n)
{
    vector<vector<double>> yarray;
    yarray.resize(n);
    vector<double> xarray;
    xarray.resize(n);

    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            double value = 0.6*(1/pow(l, 4))*((x[i]*x[i]) - (l*l))*((y[j]*y[j]) - (l*l));
            xarray[i] = value;
        }
        yarray[j] = xarray;
    }
    return yarray;
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

auto conduction(vector<vector<double>> tau_0, double d, int treshold, double dx, double dy, double dt, int n, int nt)
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
            //xarray[n-1] = xarray[n-2];
            xarray[n-1] = tau[k-1][j][n-1];
            yarray[j] = xarray;
        }
        //yarray[n-1] = xarray;
        yarray[n-1] = tau[k-1][n-1];
        tau[k] = yarray;
        if (k%((int)round(nt/240)) == 0)
        {
            write_data(yarray, nt, n, (int)round(k/240));
            system("python3 plot-a-frame.py");
        }
    }
    return tau;
}

auto wave_eq(vector<vector<double>> tau_0, double d, int treshold, double dx, double dy, double dt, int n, int nt)
{
    vector<vector<vector<double>>> tau;
    tau.resize(nt);
    vector<vector<double>> yarray;
    yarray.resize(n);
    vector<double> xarray;
    xarray.resize(n);

    tau[0] = tau_0;
    tau[1] = tau_0;

    vector<vector<double>> empty;
    empty.resize(1);
    vector<double> empty_;
    empty_.resize(1);

    empty_[0] = 0;
    empty[0] = empty_;

    for(int k = 2; k < nt; k++)
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
                xarray[i] = (2*tau[k-1][j][i] - tau[k-2][j][i]) + ((d * (((tau[k-1][j][i+1] - 2*tau[k-1][j][i] + tau[k-1][j][i-1]) / (dx*dx)) + ((tau[k-1][j+1][i] - 2*tau[k-1][j][i] + tau[k-1][j-1][i]) / (dy*dy)) ) ) *dt*dt);
            }
            //xarray[n-1] = xarray[n-2];
            xarray[n-1] = tau[k-1][j][n-1];
            yarray[j] = xarray;
        }
        //yarray[n-1] = xarray;
        yarray[n-1] = tau[k-1][n-1];
        tau[k] = yarray;
        if (k%((int)round(nt/240)) == 0)
        {
            write_data(yarray, nt, n, (int)round(k/240));
            system("python3 plot-a-frame.py");
        }
    }
    return tau;
}

auto wave_eq_a(vector<vector<double>> tau_0, double d, double q, int treshold, double dx, double dy, double dt, int n, int nt)
{
    vector<vector<vector<double>>> tau;
    tau.resize(nt);
    vector<vector<double>> yarray;
    yarray.resize(n);
    vector<double> xarray;
    xarray.resize(n);

    tau[0] = tau_0;
    tau[1] = tau_0;

    vector<vector<double>> empty;
    empty.resize(1);
    vector<double> empty_;
    empty_.resize(1);

    empty_[0] = 0;
    empty[0] = empty_;

    for(int k = 2; k < nt; k++)
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
                xarray[i] = (2*tau[k-1][j][i] - tau[k-2][j][i]) + ((d * (((tau[k-1][j][i+1] - 2*tau[k-1][j][i] + tau[k-1][j][i-1]) / (dx*dx)) + ((tau[k-1][j+1][i] - 2*tau[k-1][j][i] + tau[k-1][j-1][i]) / (dy*dy)) ) - (q * ((tau[k-1][j][i] - tau[k-2][j][i])/dt)) ) *dt*dt);
            }
            //xarray[n-1] = xarray[n-2];
            xarray[n-1] = tau[k-1][j][n-1];
            yarray[j] = xarray;
        }
        //yarray[n-1] = xarray;
        yarray[n-1] = tau[k-1][n-1];
        tau[k] = yarray;
        if (k%((int)round(nt/240)) == 0)
        {
            write_data(yarray, nt, n, (int)round(k/240));
            system("python3 plot-a-frame.py");
        }
    }
    return tau;
}

auto wave_eq_aq(vector<vector<double>> tau_0, double d, double q, int treshold, double dx, double dy, double dt, int n, int nt)
{
    vector<vector<vector<double>>> tau;
    tau.resize(nt);
    vector<vector<double>> yarray;
    yarray.resize(n);
    vector<double> xarray;
    xarray.resize(n);

    tau[0] = tau_0;
    tau[1] = tau_0;

    vector<vector<double>> empty;
    empty.resize(1);
    vector<double> empty_;
    empty_.resize(1);

    empty_[0] = 0;
    empty[0] = empty_;

    for(int k = 2; k < nt; k++)
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
                xarray[i] = (2*tau[k-1][j][i] - tau[k-2][j][i]) + ((d * (((tau[k-1][j][i+1] - 2*tau[k-1][j][i] + tau[k-1][j][i-1]) / (dx*dx)) + ((tau[k-1][j+1][i] - 2*tau[k-1][j][i] + tau[k-1][j-1][i]) / (dy*dy)) ) - (q * pow(((tau[k-1][j][i] - tau[k-2][j][i])/dt), 2)) ) *dt*dt);
            }
            //xarray[n-1] = xarray[n-2];
            xarray[n-1] = tau[k-1][j][n-1];
            yarray[j] = xarray;
        }
        //yarray[n-1] = xarray;
        yarray[n-1] = tau[k-1][n-1];
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
    string DE;
    string IC;
    string EM;
    bool asked = false;
    vector<vector<double>> tau_0;
    int n = 200;
    int nt = 80000;
    int treshold = 40;
    double tmax = 60;
    double l = 5;
    double a = 2;
    double d = 0.2;
    double k = 0.3;
    vector<double> x = linspace(-l, l, n);
    vector<double> y = linspace(-l, l, n);
    vector<double> t = linspace(0, tmax, nt);

    double dx = x[1] - x[0];
    double dy = y[1] - y[0];
    double dt = t[1] - t[0];
    cout << dx << "  " << dy << "  " << dt << endl;

    cout << "" << endl;
    cout << "What kind of differential equation do you want to solve ?" << endl;
    cout << "" << endl;
    cout << "1 - Conduction [Heat equation]" << endl;
    cout << "2 - Wave [Wave equation]" << endl;
    cout << "3 - Wave with linear energy lost [Wave equation]" << endl;
    cout << "4 - Wave with quadratic energy lost [Wave equation]" << endl;
    cout << "" << endl;
    cin >> DE;

    cout << "" << endl;
    cout << "What kind of initial condition do you want to use ?" << endl;
    cout << "" << endl;
    cout << "1 - gaussian" << endl;
    cout << "2 - sinusoidal with gaussian" << endl;
    cout << "3 - sinusoidal" << endl;
    cout << "4 - quadratic" << endl;
    cout << "" << endl;
    cin >> IC;
    cout << "" << endl;

    cout << "" << endl;
    cout << "Do you want to use the edit mode ? [y/N]" << endl;
    cout << "" << endl;
    cin >> EM;

    if (EM == "y")
    {
        cout << "" << endl;
        cout << "a (size coef) [double/float] = " << endl;
        cout << "" << endl;
        cin >> a;
        cout << "" << endl;
        cout << "d (time coef) [double/float] = " << endl;
        cout << "" << endl;
        cin >> d;
        cout << "" << endl;
        if (DE == "3")
        {
            cout << "" << endl;
            cout << "k (energy lost coef) = " << endl;
            cout << "" << endl;
            cin >> k;
        }
        if (DE == "4")
        {
            cout << "" << endl;
            cout << "k (energy lost coef) = " << endl;
            cout << "" << endl;
            cin >> k;
        }
        asked = true;
    }
    if (EM == "Y")
    {
        cout << "" << endl;
        cout << "a (size coef) [double/float] = " << endl;
        cout << "" << endl;
        cin >> a;
        cout << "" << endl;
        cout << "d (time coef) [double/float] = " << endl;
        cout << "" << endl;
        cin >> d;
        cout << "" << endl;
        if (DE == "3")
        {
            cout << "" << endl;
            cout << "k (energy lost coef) = " << endl;
            cout << "" << endl;
            cin >> k;
        }
        asked = true;
    }
    if (asked == false)
    {
        cout << "" << endl;
        cout << "Edit mode off" << endl;
        cout << "" << endl;
    }

    if (IC == "1")
    {
        tau_0 = tau_init_gaussian(x, y, l, a, n);
    }
    if (IC == "2")
    {
        tau_0 = tau_init_sin_gaussian(x, y, l, a, n);
    }
    if (IC == "3")
    {
        tau_0 = tau_init_sin(x, y, l, a, n);
    }
    if (IC == "4")
    {
        tau_0 = tau_init_q(x, y, l, a, n);
    }

    if (DE == "1")
    {
        vector<vector<vector<double>>> tau = conduction(tau_0, d, treshold, dx, dy, dt, n, nt);
        cout << "running the conduction equation" << endl;
    }
    if (DE == "2")
    {
        vector<vector<vector<double>>> tau = wave_eq(tau_0, d, treshold, dx, dy, dt, n, nt);
        cout << "running the wave equation" << endl;
    }
    if (DE == "3")
    {
        vector<vector<vector<double>>> tau = wave_eq_a(tau_0, d, k, treshold, dx, dy, dt, n, nt);
        cout << "running the wave equation with linear energy lost" << endl;
    }
    if (DE == "4")
    {
        vector<vector<vector<double>>> tau = wave_eq_aq(tau_0, d, k, treshold, dx, dy, dt, n, nt);
        cout << "running the wave equation with quadratic energy lost" << endl;
    }


    cout << "done" << endl;
    cout << "now writing data into a file" << endl;
    return 0;
}
