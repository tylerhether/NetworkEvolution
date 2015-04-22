//
//  main.cpp
//  main
//
//  Created by Tyler Hether on 4/17/15.
//
//

#include <iostream> // input output
// #include <stdafx> // only used for Windows compiling
#include <cmath> // math operations and (psuedo)random number generation
#include <new>
#include <memory>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <math.h>
#include <complex>
#include <algorithm>
#include <fstream>
#include <cstdlib> // allows us to halt with exit() function and use things like atoi
#include <random>


using namespace std; // This let's you use the shorthand the std library's functions

//Clasess and Strcutures
struct locus                                // Initialize the locus structure
{
    string regulatory;                         // regulatory alleles can only be the integers 0,1,2
    double coding;                          // coding alleles are continuous
};

class Populations                            // Initialize the population class
{
//private:                                  // Shouldn't we have some private members?
public:
    // Define the structure of the population.
    locus**** pop;
    
    // Initialize variables that will make numPops populations
    int numPops,numInd,numChromo,numLoci;
    
    // Start with a NULL population
    Populations();
    
    // Initialize population
    void initilizePop(string reg_pattern, double theta, double gamma, char mod, double a1, double a2, double allelic_Stdev);
    
    // Destroy population. Can't we just used a destructor here?
    void deletePop();
    
    // Print population to screen (for debugging)
    void printPop();
    
};

// Forward Declarations of Functions:       // Can be moved to a header file if gets too long.
void input(Populations *popPtr, int POPS, int INDS);             // These are function prototypes
void printLocus(locus Locus);               // These are function prototypes
void Pheno_to_Geno(string reg_pattern, double x, double y, double theta, double gamma, char *mod, double &a1, double &a2);
double make_genos(double geno_value, double allelic_stdev);


// Main Function to run:
int main(int argc, char *argv[])
{

    // Pull in command line arguments
    if (argc != 11) {
        // Inform the user of how to use the program if not entered in correctly
        std::cout << "Usage is <num_pops> <num_individuals> <initial_x> <initial_y> <initial_reg_pattern> <theta> <gamma> <model> <num_generations> <allelic_stdev>\n";
        std::cout << argc-1 << " given" << endl;
        exit(0);
    }
    
    int numPops(atoi(argv[1]));                 // Number of Populations
    int numInds(atoi(argv[2]));                 // Number of Individuals per population
    double x1(atoi(argv[3]));                   // Initial mean trait value for trait # 1
    double x2(atoi(argv[4]));                   // Initial mean trait value for trait # 2
    char* reg_pattern = argv[5];                // r11, r12, r21, r22 where rij is the reg allele for locus i and allele j
    double theta(atoi(argv[6]));                // Theta value for network regulation
    double gamma(atoi(argv[7]));                // gamma (decay) rate for network regulation
    char *mod = argv[8];                        // model used: 'A' or 'B'
    int num_generations(atoi(argv[9]));         // Number of generations to simulate
    double allelic_Stdev(atoi(argv[10]));       // Initial standard dev. of allelic values to seed standing genetic variation
    
    
    
    
    
//    srand((unsigned int)time(NULL));        // Seeding Random
    Populations Pop;                     // Initialize Populations
    
    // Running:
    input(&Pop, numPops, numInds);

    // Initialize genotypic values. These will be updated in the next function.
    double a1, a2;
    
    // Find the genotypic values that make the starting two-trait phenotype
    Pheno_to_Geno(reg_pattern, x1, x2, theta, gamma, mod, a1, a2);

    // Initialize the populations:
    Pop.initilizePop(reg_pattern, theta, gamma, *mod, a1, a2, allelic_Stdev); // TO DO: remove hard-coded variances in random normal generation
    

    // Recursion:
    for(int g=0; g<num_generations; g++){
        Pop.printPop();             // Troubleshooting: print populations
        
        // Need selection
        
        // Need mutation
        
        // Need mating
        
        // Need migration
        
        // Need output summary stats to file
        
        
    }

   
    // Cleaning up:
    Pop.deletePop();
    
    /* Troubleshooting: run this block of code to determine the version of C++ that is used:
     * if( __cplusplus == 201103L ) std::cout << "C++11\n" ;
     * else if( __cplusplus == 199711L ) std::cout << "C++98\n" ;
     * else std::cout << "pre-standard C++\n" ;*/
    
    return 0;
}




//Class Functions
Populations::Populations()
{
    pop=NULL;
}

void Populations::initilizePop(string reg_pattern, double theta, double gamma, char mod, double a1, double a2, double allelic_Stdev)
{
    pop=new locus***[numPops];
    for(int p=0; p<numPops; p++){
        pop[p]=new locus**[numInd];
        for(int i=0; i<numInd; i++){
            pop[p][i]=new locus*[numChromo];
            int counter(0);
            for(int c=0;c<numChromo;c++){
                pop[p][i][c]=new locus[numLoci];

                for(int l=0; l<numLoci;l++){
                    // Need to make genotype for the individual
                    
                    if(c==0){
                        pop[p][i][c][l].coding=a1+make_genos(a1/2, allelic_Stdev);
                        pop[p][i][c][l].regulatory=static_cast<char>(reg_pattern[counter]);
                    }
                    if(c==1){
                        pop[p][i][c][l].coding=a2+make_genos(a2/2, allelic_Stdev);
                        pop[p][i][c][l].regulatory=static_cast<char>(reg_pattern[counter]);
                    }
                    counter++;
                    if(counter==4){
                        counter=0;
                    }
                    
//                    pop[p][i][c][l].coding=gamma;           // Seed coding allele expression
//                    pop[p][i][c][l].regulatory=rand() % 3; // Seed regulatory allele
                }
            }
        }
    }
    
}

void Populations::deletePop()
{
    for(int p=0; p<numPops; p++){
        for(int i=0;i<numInd;i++){
            for(int c=0;c<numChromo;c++){
                delete[] pop[p][i][c];
            }
            delete[] pop[p][i];
        }
        delete[] pop[p];
    }
    delete[] pop;
}
void Populations::printPop()
{
    for(int p=0; p<numPops; p++){
        for(int i=0; i<numInd; i++){
            for(int c=0;c<numChromo;c++){
                for(int l=0; l<numLoci;l++){
                    cout<< "Pop=" << p << "\t" << "Ind=" << i << "\t" << "chrom=" << c << "\t"<< "locus=" << l << "\t(" << pop[p][i][c][l].regulatory<<" , "<<pop[p][i][c][l].coding<<")" << endl;
                }
//                cout<<endl;
            }
            cout<<endl;
        }
    }
}


//General Functions
void input(Populations *popPtr, int POPS, int INDS)
{
    
    (*popPtr).numPops = POPS;
    cout<<POPS<<endl;
    cout<<"# of populations are "<<(*popPtr).numPops<<endl;
    (*popPtr).numInd = INDS;
    cout<<"# of individuals are "<<(*popPtr).numInd<<endl;
    (*popPtr).numChromo = 2;
    cout<<"# of chromosomes are "<<(*popPtr).numChromo<<endl;
    (*popPtr).numLoci =2;
    cout<<"# of loci are "<<(*popPtr).numLoci<<endl;
    
}

void printLocus(locus Locus)
{
    cout<<"Locus"<<endl;
    cout<<"regulatory =" <<Locus.regulatory<<endl;
    cout<<"coding = " <<Locus.coding<<endl;
    
}

void Pheno_to_Geno(string reg_pattern, double x, double y, double theta, double gamma, char *mod, double &a1, double &a2)
{
    switch(*mod)
    {
            case 'A':
            if(reg_pattern=="0000"){
                a1 = (gamma*x*(theta+y))/theta;
                a2 = (gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="0010"){
                a1 = (gamma*x*(theta+y))/theta;
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="0020"){
                a1 = (gamma*x*(theta+y))/theta;
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="0001"){
                a1 = (gamma*x*(theta+y))/theta;
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="0011"){
                a1 = (gamma*x*(theta+y))/theta;
                a2 = gamma*y;
            }
            if(reg_pattern=="0021"){
                a1 = (gamma*x*(theta+y))/theta;
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="0002"){
                a1 = (gamma*x*(theta+y))/theta;
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="0012"){
                a1 = (gamma*x*(theta+y))/theta;
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="0022"){
                a1 = (gamma*x*(theta+y))/theta;
                a2 = (gamma*(theta+x)*y)/x;
            }
            if(reg_pattern=="1000"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = (gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="1010"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="1020"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="1001"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="1011"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = gamma*y;
            }
            if(reg_pattern=="1021"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="1002"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="1012"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="1022"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = (gamma*(theta+x)*y)/x;
            }
            if(reg_pattern=="2000"){
                a1 = 2*gamma*x;
                a2 = (gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="2010"){
                a1 = 2*gamma*x;
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="2020"){
                a1 = 2*gamma*x;
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="2001"){
                a1 = 2*gamma*x;
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="2011"){
                a1 = 2*gamma*x;
                a2 = gamma*y;
            }
            if(reg_pattern=="2021"){
                a1 = 2*gamma*x;
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="2002"){
                a1 = 2*gamma*x;
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="2012"){
                a1 = 2*gamma*x;
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="2022"){
                a1 = 2*gamma*x;
                a2 = (gamma*(theta+x)*y)/x;
            }
            if(reg_pattern=="0100"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = (gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="0110"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="0120"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="0101"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="0111"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = gamma*y;
            }
            if(reg_pattern=="0121"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="0102"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="0112"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="0122"){
                a1 = (2*gamma*x*(theta+y))/(2*theta+y);
                a2 = (gamma*(theta+x)*y)/x;
            }
            if(reg_pattern=="1100"){
                a1 = gamma*x;
                a2 = (gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="1110"){
                a1 = gamma*x;
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="1120"){
                a1 = gamma*x;
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="1101"){
                a1 = gamma*x;
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="1111"){
                a1 = gamma*x;
                a2 = gamma*y;
            }
            if(reg_pattern=="1121"){
                a1 = gamma*x;
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="1102"){
                a1 = gamma*x;
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="1112"){
                a1 = gamma*x;
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="1122"){
                a1 = gamma*x;
                a2 = (gamma*(theta+x)*y)/x;
            }
            if(reg_pattern=="2100"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = (gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="2110"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="2120"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="2101"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="2111"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = gamma*y;
            }
            if(reg_pattern=="2121"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="2102"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="2112"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="2122"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = (gamma*(theta+x)*y)/x;
            }
            if(reg_pattern=="0200"){
                a1 = 2*gamma*x;
                a2 = (gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="0210"){
                a1 = 2*gamma*x;
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="0220"){
                a1 = 2*gamma*x;
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="0201"){
                a1 = 2*gamma*x;
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="0211"){
                a1 = 2*gamma*x;
                a2 = gamma*y;
            }
            if(reg_pattern=="0221"){
                a1 = 2*gamma*x;
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="0202"){
                a1 = 2*gamma*x;
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="0212"){
                a1 = 2*gamma*x;
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="0222"){
                a1 = 2*gamma*x;
                a2 = (gamma*(theta+x)*y)/x;
            }
            if(reg_pattern=="1200"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = (gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="1210"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="1220"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="1201"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="1211"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = gamma*y;
            }
            if(reg_pattern=="1221"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="1202"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="1212"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="1222"){
                a1 = (2*gamma*x*(theta+y))/(theta+2*y);
                a2 = (gamma*(theta+x)*y)/x;
            }
            if(reg_pattern=="2200"){
                a1 = (gamma*x*(theta+y))/y;
                a2 = (gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="2210"){
                a1 = (gamma*x*(theta+y))/y;
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="2220"){
                a1 = (gamma*x*(theta+y))/y;
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="2201"){
                a1 = (gamma*x*(theta+y))/y;
                a2 = (2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="2211"){
                a1 = (gamma*x*(theta+y))/y;
                a2 = gamma*y;
            }
            if(reg_pattern=="2221"){
                a1 = (gamma*x*(theta+y))/y;
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="2202"){
                a1 = (gamma*x*(theta+y))/y;
                a2 = 2*gamma*y;
            }
            if(reg_pattern=="2212"){
                a1 = (gamma*x*(theta+y))/y;
                a2 = (2*gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="2222"){
                a1 = (gamma*x*(theta+y))/y;
                a2 = (gamma*(theta+x)*y)/x;
            }
            break;
            
            case 'B':
            if(reg_pattern=="0000"){
                a1=(gamma*x*(theta+y))/theta;
                a2=(gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="0010"){
                a1=(gamma*x*(theta+y))/theta;
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="0020"){
                a1=(gamma*x*(theta+y))/theta;
                a2=gamma*y;
            }
            if(reg_pattern=="0001"){
                a1=(gamma*x*(theta+y))/theta;
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="0011"){
                a1=(gamma*x*(theta+y))/theta;
                a2=gamma*y;
            }
            if(reg_pattern=="0021"){
                a1=(gamma*x*(theta+y))/theta;
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="0002"){
                a1=(gamma*x*(theta+y))/theta;
                a2=gamma*y;
            }
            if(reg_pattern=="0012"){
                a1=(gamma*x*(theta+y))/theta;
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="0022"){
                a1=(gamma*x*(theta+y))/theta;
                a2=(gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="1000"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=(gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="1010"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="1020"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=gamma*y;
            }
            if(reg_pattern=="1001"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="1011"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=gamma*y;
            }
            if(reg_pattern=="1021"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="1002"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=gamma*y;
            }
            if(reg_pattern=="1012"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="1022"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=(gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="2000"){
                a1=gamma*x;
                a2=(gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="2010"){
                a1=gamma*x;
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="2020"){
                a1=gamma*x;
                a2=gamma*y;
            }
            if(reg_pattern=="2001"){
                a1=gamma*x;
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="2011"){
                a1=gamma*x;
                a2=gamma*y;
            }
            if(reg_pattern=="2021"){
                a1=gamma*x;
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="2002"){
                a1=gamma*x;
                a2=gamma*y;
            }
            if(reg_pattern=="2012"){
                a1=gamma*x;
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="2022"){
                a1=gamma*x;
                a2=(gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="0100"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=(gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="0110"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="0120"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=gamma*y;
            }
            if(reg_pattern=="0101"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="0111"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=gamma*y;
            }
            if(reg_pattern=="0121"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="0102"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=gamma*y;
            }
            if(reg_pattern=="0112"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="0122"){
                a1=(2*gamma*x*(theta+y))/(2*theta+y);
                a2=(gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="1100"){
                a1=gamma*x;
                a2=(gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="1110"){
                a1=gamma*x;
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="1120"){
                a1=gamma*x;
                a2=gamma*y;
            }
            if(reg_pattern=="1101"){
                a1=gamma*x;
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="1111"){
                a1=gamma*x;
                a2=gamma*y;
            }
            if(reg_pattern=="1121"){
                a1=gamma*x;
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="1102"){
                a1=gamma*x;
                a2=gamma*y;
            }
            if(reg_pattern=="1112"){
                a1=gamma*x;
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="1122"){
                a1=gamma*x;
                a2=(gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="2100"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=(gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="2110"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="2120"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=gamma*y;
            }
            if(reg_pattern=="2101"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="2111"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=gamma*y;
            }
            if(reg_pattern=="2121"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="2102"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=gamma*y;
            }
            if(reg_pattern=="2112"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="2122"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=(gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="0200"){
                a1=gamma*x;
                a2=(gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="0210"){
                a1=gamma*x;
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="0220"){
                a1=gamma*x;
                a2=gamma*y;
            }
            if(reg_pattern=="0201"){
                a1=gamma*x;
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="0211"){
                a1=gamma*x;
                a2=gamma*y;
            }
            if(reg_pattern=="0221"){
                a1=gamma*x;
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="0202"){
                a1=gamma*x;
                a2=gamma*y;
            }
            if(reg_pattern=="0212"){
                a1=gamma*x;
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="0222"){
                a1=gamma*x;
                a2=(gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="1200"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=(gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="1210"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="1220"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=gamma*y;
            }
            if(reg_pattern=="1201"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="1211"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=gamma*y;
            }
            if(reg_pattern=="1221"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="1202"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=gamma*y;
            }
            if(reg_pattern=="1212"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="1222"){
                a1=(2*gamma*x*(theta+y))/(2*theta+3*y);
                a2=(gamma*(theta+x)*y)/(theta+2*x);
            }
            if(reg_pattern=="2200"){
                a1=(gamma*x*(theta+y))/(theta+2*y);
                a2=(gamma*(theta+x)*y)/theta;
            }
            if(reg_pattern=="2210"){
                a1=(gamma*x*(theta+y))/(theta+2*y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="2220"){
                a1=(gamma*x*(theta+y))/(theta+2*y);
                a2=gamma*y;
            }
            if(reg_pattern=="2201"){
                a1=(gamma*x*(theta+y))/(theta+2*y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+x);
            }
            if(reg_pattern=="2211"){
                a1=(gamma*x*(theta+y))/(theta+2*y);
                a2=gamma*y;
            }
            if(reg_pattern=="2221"){
                a1=(gamma*x*(theta+y))/(theta+2*y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="2202"){
                a1=(gamma*x*(theta+y))/(theta+2*y);
                a2=gamma*y;
            }
            if(reg_pattern=="2212"){
                a1=(gamma*x*(theta+y))/(theta+2*y);
                a2=(2*gamma*(theta+x)*y)/(2*theta+3*x);
            }
            if(reg_pattern=="2222"){
                a1=(gamma*x*(theta+y))/(theta+2*y);
                a2=(gamma*(theta+x)*y)/(theta+2*x);
            }

            break;
    }
    
    
    cout << "reg = " << reg_pattern << " x = " << x << " y = " << y << " theta = " << theta << " gamma= " << gamma << " model= " << mod << " a1= " << a1 << " a2= " << a2 << endl;
    x = 20.3;
}

double make_genos(double geno_value, double allelic_stdev)
{
    using namespace std;
    random_device rd;
    mt19937 e2(rd());
    normal_distribution<double> dist(0, 20);
    
//    cout << "Here's the result: " <<  dist(e2) << endl;
    return dist(e2);
}























