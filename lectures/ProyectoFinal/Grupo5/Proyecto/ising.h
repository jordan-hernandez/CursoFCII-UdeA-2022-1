#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>

#define UP 0
#define RIGHT 1
#define LEFT 2
#define DOWN 3

#define L 16
#define SIZE L*L

#define DATA 9

#define MAG 0
#define MAG2 1
#define MAG4 2
#define MAGERR 3
#define SUSERR 4
#define ENE 5
#define ENE2 6
#define ENE4 7
#define ENERR 8
#define CHERR 9



using namespace std;

class Ising
{
bool  spins[SIZE]; //Stores the spins
    
int i,j; //Counters

int N; //Number of averages done into the system

double h[5]; //Values of the exp(-J/kT)

double energy; //Value of the energy of the system

double m[DATA]; //This array contains several moments of the magnetization and energy

int neighs[SIZE][4]; //To store nearest neighbours

mt19937 gen; //Mersenne Twister RNG

uniform_int_distribution<int> brandom; //Get any random integer
uniform_int_distribution<int> ran_pos; //Get any random integer
uniform_real_distribution<double> ran_u; //Our uniform variable generator

double tstar; //Control parameter
double deltat, deltat_crit; //Change in control parameter by iteration
double tmax, tmin, tcrit_up, tcrit_down; //Max and min temperature, and interval where we apply Wolff

 /*  
 tmax = 5.0;
 tmin = 0.1;
 tcrit_up = 2.4;
 tcrit_down = 2.2;

 deltat = 0.1;
 deltat_crit = 0.01;
*/
ofstream output; //Output of the stream

public:
Ising():gen(958431198),brandom(0, 1),ran_pos(0, SIZE-1),ran_u(0.0, 1.0){}

void initialize(bool spins[SIZE], mt19937& gen, uniform_int_distribution<int>& brandom);
void get_neighbors(int neighs[SIZE][4]);


void do_step(bool spins[SIZE],  int neighs[SIZE][4], double tstar, int N, double h[5], double& energy, mt19937& gen, uniform_real_distribution<double>& ran_u, uniform_int_distribution<int>& ran_pos, double m[DATA]);
void do_step_wolff(bool spins[SIZE],  int neighs[SIZE][4], double tstar, int N, double& energy, mt19937& gen, uniform_real_distribution<double>& ran_u, uniform_int_distribution<int>& ran_pos, double m[DATA]);
double magnetization(bool spins[SIZE]);

void flip_spin(bool spins[SIZE], int neighs[SIZE][4], double h[5], double& energy, mt19937& gen, uniform_real_distribution<double>& ran_u, uniform_int_distribution<int>& ran_pos);
void add_to_cluster(bool spins[SIZE], int neighs[SIZE][4], int pos, double& energy, double p, mt19937& gen, uniform_real_distribution<double>& ran_u);
double get_energy(bool spins[SIZE], int neighs[SIZE][4]);

void write(bool spins[SIZE]);

void w_output(ofstream& file, double tstar, int N, double m[DATA]);

 

private:
bool spin;



/*
void initialize(bool , mt19937& , uniform_int_distribution<int>& );
void get_neighbors(int );


void do_step(bool ,  int , double , int , double , double& , mt19937& , uniform_real_distribution<double>& , uniform_int_distribution<int>& , double );
void do_step_wolff(bool ,  int , double , int , double& , mt19937& , uniform_real_distribution<double>& , uniform_int_distribution<int>& , double );
double magnetization(bool );

void flip_spin(bool , int , double , double& , mt19937& , uniform_real_distribution<double>& , uniform_int_distribution<int>& );
void add_to_cluster(bool , int , int , double& , double , mt19937& , uniform_real_distribution<double>& );
double get_energy(bool , int );

void write(bool );

void w_output(ofstream& , double , int , double );
*/

 
};
