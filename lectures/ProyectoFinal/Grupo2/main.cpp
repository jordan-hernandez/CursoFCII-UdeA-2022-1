#include "vanderpol.h"
using namespace std;

int main(){
    
    Vanderpol ejemplo("Oscilador de Van der pol");
    
    // void disparo(a, b, alpha, beta, N, TOL, M, mu)
    
    //cout<<"Ejemplo del libro\n"<<endl;
    //ejemplo.disparo(1, 3, 17, 43/3, 20, 10e-5, 10, 1);
    
    //ec de Duffing
    //ejemplo.disparo(1, 2, 0, 0, 50, 10e-2, 10000, 0);
    
    // mu=0
    ejemplo.disparo(1, 2, 0, 1, 50, 10e-2, 10000, 0);
    
    // mu=0.2
    ejemplo.disparo(1, 2, 0, 1, 50, 10e-2, 10000, 0.2);
    
    // mu=0.5
    ejemplo.disparo(1, 2, 0, 1, 50, 10e-2, 10000, 0.5);
    
    // mu=1
    ejemplo.disparo(1, 2, 0, 1, 50, 10e-2, 10000, 1);
    
    // mu=2
    ejemplo.disparo(1, 2, 0, 1, 50, 10e-2, 10000, 2);
    
    // mu=4
    ejemplo.disparo(1, 2, 0, 1, 50, 10e-2, 10000, 4);
    
    // mu=6
    ejemplo.disparo(1, 2, 0, 1, 50, 10e-2, 10000, 6);
    
    // mu=8
    ejemplo.disparo(1, 2, 0, 1, 50, 10e-2, 10000, 8);
    
    // mu=10
    ejemplo.disparo(1, 2, 0, 1, 50, 10e-2, 10000, 10);
    
    cout<<"Ejecucion finalizada"<<endl;
    

    
//system("pause");
return 0;
        
}
