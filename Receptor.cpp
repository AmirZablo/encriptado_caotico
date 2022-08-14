#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iomanip>

using namespace std;

//Global parameters
double h = 0.0001;
double sigma=10;
double b=8.0/3.0;
double r=60;
double u0=4;
double v0=4;
double w0=4;
double epsilon=1830300;

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
double k2_u(double u, double v,double w, double k1u, double k1v, double k1w)
{
    double k2 = sigma*((v+(h*k1v)/2)-(u+(h*k1u)/2));

    return k2;
}

//Compute k2_v
double k2_v(double u, double v,double w, double k1u, double k1v, double k1w)
{
    double k2=r*(u+(h*k1u)/2)-(v+(h*k1v)/2)-20*(u+(h*k1u)/2)*(w+(h*k1w)/2);

    return k2;
}

//Compute k2_w
double k2_w(double u, double v,double w, double k1u, double k1v, double k1w)
{
    double k2=5*(u+(h*k1_u(u,v,w))/2)*(v+(h*k1_v(u,v,w))/2)-b*(w+(h*k1_w(u,v,w))/2);

    return k2;
}

//Compute k3_u
double k3_u(double u, double v,double w, double k2u, double k2v, double k2w)
{
    double k3 =sigma*((v+(h*k2v)/2)-(u+(h*k2u)/2));

    return k3;
}

//Compute k3_v
double k3_v(double u, double v,double w, double k2u, double k2v, double k2w)
{
    double k3=r*(u+(h*k2u)/2)-(v+(h*k2v)/2)-20*(u+(h*k2u)/2)*(w+(h*k2w)/2);

    return k3;
}

//Compute k3_w
double k3_w(double u, double v,double w, double k2u, double k2v, double k2w)
{
    double k3=5*(u+(h*k2u)/2)*(v+(h*k2v)/2)-b*(w+(h*k2w)/2);

    return k3;
}

//Compute k4_u
double k4_u(double u, double v,double w, double k3u, double k3v, double k3w)
{
    double k4 =sigma*((v+(h*k3v))-(u+(h*k3u)));

    return k4;
}

//Compute k4_v
double k4_v(double u, double v,double w, double k3u, double k3v, double k3w)
{
    double k4=r*(u+(h*k3u))-(v+(h*k3v))-20*(u+(h*k3u))*(w+(h*k3w));

    return k4;
}

//Compute k4_w
double k4_w(double u, double v,double w, double k3u, double k3v, double k3w)
{
    double k4=5*(u+(h*k3u))*(v+(h*k3v))-b*(w+(h*k3w));

    return k4;
}

void RK4(double &newU, double &newV, double &newW, double u, double v, double w, double s)
{
    double k1u=k1_u(u,v,w);
    double k1v=k1_v(s,v,w);
    double k1w=k1_w(s,v,w);

    double k2u=k2_u(u,v,w,k1u,k1v,k1w);
    double k2v=k2_v(s,v,w,k1u,k1v,k1w);
    double k2w=k2_w(s,v,w,k1u,k1v,k1w);

    double k3u=k3_u(u,v,w,k2u,k2v,k2w);
    double k3v=k3_v(s,v,w,k2u,k2v,k2w);
    double k3w=k3_w(s,v,w,k2u,k2v,k2w);

    double k4u=k4_u(u,v,w,k3u,k3v,k3w);
    double k4v=k4_v(s,v,w,k3u,k3v,k3w);
    double k4w=k4_w(s,v,w,k3u,k3v,k3w);

    newU = u + (h/6)*(k1u + 2*k2u + 2*k3u + k4u);
    newV = v + (h/6)*(k1v + 2*k2v + 2*k3v + k4v);
    newW = w + (h/6)*(k1w + 2*k2w + 2*k3w + k4w);
}


int main()
{
    //Vectors
    vector<double> u; //ur(t)
    vector<double> v; //vr(t)
    vector<double> w; //wr(t)

    vector<double> s; //s(t)
    vector<double> m; //m(t)

    /*//Emisor's trajectory
    vector<double> u_t;
    vector<double> v_t;
    vector<double> w_t;

    //Load emisor's trajectory
    ifstream forzado("forzado.dat");
    while(!forzado.eof())
    {
        double ut;
        double vt;
        double wt;
        double descarte;
        forzado>>descarte;
        forzado>>ut;
        forzado>>vt;
        forzado>>wt;
        u_t.push_back(ut);
        v_t.push_back(vt);
        w_t.push_back(wt);
    }*/

    //Load s(t)
    ifstream encriptado("s(t).dat");
    while(!encriptado.eof())
    {
        double s_t;
        double descarte;
        encriptado>>descarte;
        encriptado>>s_t;
        s.push_back(s_t);
    }

    //Save initial conditions
    u.push_back(u0);
    v.push_back(v0);
    w.push_back(w0);

    /*//Compute ur(t), vr(t) y wr(t) with their respectives RK4 (forcing w)
    for(int cont=0;cont<u_t.size();cont++)
    {
        double newU = RK4_u(u.back(),v.back(),w_t.at(cont));
        double newV = RK4_v(u.back(),v.back(),w_t.at(cont));
        double newW = RK4_w(u.back(),v.back(),w_t.at(cont));

        u.push_back(newU);
        v.push_back(newV);
        w.push_back(newW);
    }*/

    //Compute ur(t), vr(t) y wr(t) with their respectives RK4
    for(int cont=0; cont<s.size(); cont++) //cont<s.zize()
    {
        double newU=0.0;
        double newV=0.0;
        double newW=0.0;

        RK4(newU,newV,newW,u.back(),v.back(),w.back(),s.at(cont));

        u.push_back(newU);
        v.push_back(newV);
        w.push_back(newW);
    }

    //Unmask m
    for(int i=0;i<s.size();i++)
    {
        m.push_back(s.at(i)-u.at(i));
    }

    /*ofstream oFile;
    oFile.open("Trayectoria.dat");
    ofstream oFile2;
    oFile2.open("Distancia.dat");

    oFile <<"t u v w\n";
    oFile2<<"t dist\n";

    //Save treajectory and distance to emisor's trajectory
    for(int i=0; i<u_t.size(); i++)
    {
        oFile <<i*h<<" "<<setprecision(7)<<u.at(i)<<" "<<setprecision(7)<<v.at(i)<<" "<<setprecision(7)<<w.at(i)<<'\n';
        oFile2<<i*h<<" "<<setprecision(7)<<sqrt(pow(u_t.at(i)-u.at(i),2)+pow(v_t.at(i)-v.at(i),2)+pow(w_t.at(i)-w.at(i),2))<<'\n';
    }

    oFile.close();*/

    ofstream oFile3;
    oFile3.open("Desencriptado.dat");

    //Save decrypted message
    for(int i=0;i<m.size();i++)
    {
        oFile3<<setprecision(7)<<i*h<<" "<<epsilon*m.at(i)<<'\n';
    }

    return 0;
}
