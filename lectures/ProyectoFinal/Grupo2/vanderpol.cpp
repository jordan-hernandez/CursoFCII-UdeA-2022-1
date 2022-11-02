#include "vanderpol.h"
#include <iostream>
#include <fstream>

using namespace std;

Vanderpol::Vanderpol( string name )
{
cout <<"Solucion de ecuaciones no lineales por el metodo de disparo\n"<<endl;

}

double Vanderpol::funcion( const double x, const double y, const double yp){
    //ec de Van der pol
    mu = getMu();
    return mu*yp*(y*y-1)-y;
    
    //ec de Van der pol forzada; F = A sin(wx)
    //mu = getMu();
    //return mu*yp*(y*y-1)-y+10*sin(15*x);
    
    //ejemplo
    //return (32+2*x*x*x-y*yp)/8; 
    
    //Ecuacion de Duffing ypp+ayp+cy^3 = dcos(wx)
    // a = 0.02, b = 1, c = 5, d = 8 , w = 0.5
    //return 8*cos(0.5*x)-0.02*yp-5*y*y*y-y;
    
    
}

double Vanderpol::funcion_y( const double x, const double y, const double yp){
    //ec de Van der pol
    mu = getMu();
    return 2*mu*y*yp-1;
    
    //ejemplo
    //return -1*yp/8;
    
    //ec de Duffing
    //return 1+3*5*x*x;
}

double Vanderpol::funcion_yp( const double x, const double y, const double yp){
    //ec de Van der pol
    mu = getMu();
    return mu*y*y;
    
    //ejemplo
    //return -1*y/8;
    
    //ec de Duffing
    //return 0.02;
}


void Vanderpol::disparo(const double a1, const double b1 ,const double alpha1 ,const double beta1, const double N1, const double tol1, const double m1, const double mu1){

    setA(a1);
    setB(b1);
    setAlpha(alpha1);
    setBeta(beta1);
    setMu(mu1);
    
    int M = m1; //numero maximo de iteraciones
    int N = N1; //numero de subintervalos
    double TOL = tol1; //tolerancia
    
    // inicializacion de vectores solucion
    double W1[N+1];
    double W2[N+1];
    double x;
    
    // inicializacion de vectores de trabajo (Runge-Kutta)
    
    double U[3];
    double k1[3];
    double k2[3];
    double k3[3];
    double k4[3];
    double k1p[3];
    double k2p[3];
    double k3p[3];
    double k4p[3];
    
    //STEP 1
    double h = (getB()-getA())/N;
    int k = 1;
    double TK = (getBeta()-getAlpha())/(getB()-getA());
    
    //STEP 2
    
    while (k<=M){
    
    	//STEP 3
    
    	W1[0] = getAlpha();
    	W2[0] = TK;
    	U[1] = 0;
    	U[2] = 1;
    
    	//STEP 4 
    
    	for(int i = 1; i<=N ; ++i){
    		
    		// Pasos del 5 al 6 es el metodo de Runge-Kutta
    		//STEP 5
    		
    		x = getA() + (i-1)*h ;
    		
    		//STEP 6
    		
    		k1[1] = h* W2[i-1];
       		k1[2] = h* funcion(x,W1[i-1],W2[i-1]); 
		k2[1] = h* (W2[i-1]+0.5*k1[2]);
		k2[2] = h* funcion(x+0.5*h,W1[i-1]+0.5*k1[1],W2[i-1]+0.5*k1[2]);
		k3[1] = h* (W2[i-1]+0.5*k2[2]);
		k3[2] = h* funcion(x+0.5*h,W1[i-1]+0.5*k2[1],W2[i-1]+0.5*k2[2]);
		k4[1] = h* (W2[i-1]+k3[2]);
		k4[2] = h* funcion(x+h,W1[i-1]+k3[1],W2[i-1]+k3[2]);
		
		W1[i] = W1[i-1] + (k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
		W2[i] = W2[i-1] + (k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;
		
		k1p[1] = h* U[2];
		k1p[2] = h* (funcion_y(x,W1[i-1],W2[i-1])*U[1]+funcion_yp(x,W1[i-1],W2[i-1])*U[2]);
		k2p[1] = h* (U[2]+0.5*k1p[2]);
		k2p[2] = h* (funcion_y(x+0.5*h,W1[i-1],W2[i-1])*(U[1]+0.5*k1p[1])+funcion_yp(x+0.5*h,W1[i-1],W2[i-1])*(U[2]+0.5*k1p[2]));
		k3p[1] = h* (U[2]+0.5*k2p[2]);
		k3p[2] = h* (funcion_y(x+0.5*h,W1[i-1],W2[i-1])*(U[1]+0.5*k2p[1])+funcion_yp(x+0.5*h,W1[i-1],W2[i-1])*(U[2]+0.5*k2p[2]));
		k4p[1] = h* (U[2]+k3p[2]);
		k4p[2] = h* (funcion_y(x+h,W1[i-1],W2[i-1])*(U[1]+k3p[1])+funcion_yp(x+h,W1[i-1],W2[i-1])*(U[2]+k3p[2]));
		
		U[1] = U[1] + (k1p[1]+2*k2p[1]+2*k3p[1]+k4p[1])/6;
		U[2] = U[2] + (k1p[2]+2*k2p[2]+2*k3p[2]+k4p[2])/6; 
   	}
   	
   	//STEP 7
   	
   	if (abs(W1[N]-getBeta())<=TOL){
   		
   		//STEP 8
   		
		string nombreArchivo = "graficaVanderPol" + std::to_string(getMu())+ ".txt";
		ofstream archivo;
		archivo.open(nombreArchivo.c_str(), fstream::out);
		
   		
		for (int i = 0 ; i<=N ; i++ ){
   			x = getA() + i*h;
   			//OUTPUT
   		
			//archivo 
			archivo <<x<<","<<W1[i]<<","<<W2[i]<<endl; //guarda txt

			archivo << endl;
   		}
		archivo.close();

   		
   		//STEP 9
   		
   		// STOP 
   		break;
   	}
   	
   	//STEP 10
   	
   	TK = TK - (W1[N]-getBeta())/U[1];
   	k = k+1;
    
    }
    
    //STEP 11
    
    // se excede el numero maximo de iteraciones
    // no convergio
    
}

void Vanderpol::setA(const double a){
A = a;
}

void Vanderpol::setB(const double a){
B = a;
}

void Vanderpol::setAlpha(const double a){
alpha = a;
}

void Vanderpol::setBeta(const double a){
beta = a;
}

void Vanderpol::setMu(const double a){
mu = a;
}


double Vanderpol::getA(){
return A;
}

double Vanderpol::getB(){
return B;
}

double Vanderpol::getAlpha(){
return alpha;
}

double Vanderpol::getBeta(){
return beta;
}

double Vanderpol::getMu(){
return mu;
}
