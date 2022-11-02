#ifndef VANDERPOL_H
#define VANDERPOL_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <fstream>

using namespace std;

class Vanderpol{

    public:
    Vanderpol( string ); //constructor

    double funcion( const double, const double, const double);
    double funcion_y( const double, const double, const double);
    double funcion_yp( const double, const double, const double);

    void disparo(const double, const double, const double, const double, const double, const double, const double, const double);

    void setA(const double ) ;
    void setB(const double) ;
    void setAlpha(const double) ;
    void setBeta(const double) ;
    void setMu(const double) ;

    double getA( ) ;
    double getB( ) ;
    double getAlpha( ) ;
    double getBeta( ) ;
    double getMu( );

    private:
    double A;
    double B;
    double alpha;
    double beta;
    double mu;
};

#endif


