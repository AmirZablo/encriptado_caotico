#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>

using namespace std;

//Global parameters
double h = 0.0001;
double sigma=10;
double b=8.0/3.0;
double r=60;
double epsilon=1830300;
double u0=1.5;
double v0=1.5;
double w0=1.5;

//Compute k1_u
double k1_u(double u, double v,double w)
{
    double k1 = sigma*(v-u);

    return k1;
}

//Compute k1_v
double k1_v(double u, double v,double w)
{
    double k1=r*u-v-20*u*w;

    return k1;
}

//Compute k1_w
double k1_w(double u, double v,double w)
{
    double k1=5*u*v-b*w;

    return k1;
}

//Compute k2_u
double k2_u(double u, double v,double w)
{
    double k2 = sigma*((v+(h*k1_v(u,v,w))/2)-(u+(h*k1_u(u,v,w))/2));

    return k2;
}

//Compute k2_v
double k2_v(double u, double v,double w)
{
    double k2=r*(u+(h*k1_u(u,v,w))/2)-(v+(h*k1_v(u,v,w))/2)-20*(u+(h*k1_u(u,v,w))/2)*(w+(h*k1_w(u,v,w))/2);

    return k2;
}

//Compute k2_w
double k2_w(double u, double v,double w)
{
    double k2=5*(u+(h*k1_u(u,v,w))/2)*(v+(h*k1_v(u,v,w))/2)-b*(w+(h*k1_w(u,v,w))/2);

    return k2;
}

//Compute k3_u
double k3_u(double u, double v,double w)
{
    double k3 =sigma*((v+(h*k2_v(u,v,w))/2)-(u+(h*k2_u(u,v,w))/2));

    return k3;
}

//Compute k3_v
double k3_v(double u, double v,double w)
{
    double k3=r*(u+(h*k2_u(u,v,w))/2)-(v+(h*k2_v(u,v,w))/2)-20*(u+(h*k2_u(u,v,w))/2)*(w+(h*k2_w(u,v,w))/2);

    return k3;
}

//Compute k3_w
double k3_w(double u, double v,double w)
{
    double k3=5*(u+(h*k2_u(u,v,w))/2)*(v+(h*k2_v(u,v,w))/2)-b*(w+(h*k2_w(u,v,w))/2);

    return k3;
}

//Compute k4_u
double k4_u(double u, double v,double w)
{
    double k4 =sigma*((v+(h*k3_v(u,v,w)))-(u+(h*k3_u(u,v,w))));

    return k4;
}

//Compute k4_v
double k4_v(double u, double v,double w)
{
    double k4=r*(u+(h*k3_u(u,v,w)))-(v+(h*k3_v(u,v,w)))-20*(u+(h*k3_u(u,v,w)))*(w+(h*k3_w(u,v,w)));

    return k4;
}

//Compute k4_w
double k4_w(double u, double v,double w)
{
    double k4=5*(u+(h*k3_u(u,v,w)))*(v+(h*k3_v(u,v,w)))-b*(w+(h*k3_w(u,v,w)));

    return k4;
}

//Compute u_(j+1)
double RK4_u(double u, double v,double w)
{
    double RK4_u = u + (h/6)*(k1_u(u,v,w) + 2*k2_u(u,v,w) + 2*k3_u(u,v,w) + k4_u(u,v,w));
    return RK4_u;
}

//Compute v_(j+1)
double RK4_v(double u, double v,double w)
{
    double RK4_v = v + (h/6)*(k1_v(u,v,w) + 2*k2_v(u,v,w) + 2*k3_v(u,v,w) + k4_v(u,v,w));
    return RK4_v;
}

//Compute w_(j+1)
double RK4_w(double u, double v,double w)
{
    double RK4_w = w + (h/6)*(k1_w(u,v,w) + 2*k2_w(u,v,w) + 2*k3_w(u,v,w) + k4_w(u,v,w));
    return RK4_w;
}


int main()
{
    //Vectors
    vector<double> u; //u(t)
    vector<double> v; //v(t)
    vector<double> w; //w(t)

    vector<double> m; //m(t)
    vector<double> s; //s(t)

    //Save initial conditions
    u.push_back(u0);
    v.push_back(v0);
    w.push_back(w0);

    //Load and save m(t)
    ifstream song("africa.dat");
    while(!song.eof())
    {
        double n;
        song>>n;
        m.push_back(n);
    }

    //Fill u(t), v(t) and w(t) with their corresponding RK4
    for(double t=h; t<=130; t=t+h)
    {
        double newU = RK4_u(u.back(),v.back(),w.back());
        double newV = RK4_v(u.back(),v.back(),w.back());
        double newW = RK4_w(u.back(),v.back(),w.back());

        u.push_back(newU);
        v.push_back(newV);
        w.push_back(newW);
    }

    //Mask m(t)
    for(int i=0;i<m.size();i++)
    {
        s.push_back(u.at(i)+m.at(i)/epsilon);
    }

    ofstream oFile;
    oFile.open("Trayectoria.dat");

    ofstream oFile2;
    oFile2.open("s(t).dat");

    ofstream oFile3;
    oFile3.open("m(t).dat");

    oFile <<"t u v w\n";
    oFile3<<"t m\n";

    //Save trajectory, s(t) and m(t)
    for(int i=0; i<m.size(); i++)
    {
        oFile <<i*h<<" "<<setprecision(7)<<u.at(i)<<" "<<setprecision(7)<<v.at(i)<<" "<<setprecision(7)<<w.at(i)<<'\n';
        oFile2<<i*h<<" "<<setprecision(7)<<s.at(i)<<'\n';
        oFile3<<i*h<<" "<<setprecision(7)<<m.at(i)<<'\n';
    }

    oFile.close();
    oFile2.close();
    oFile3.close();

    return 0;
}
