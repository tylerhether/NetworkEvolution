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
//using namespace arma;


//Clasess and Strcutures
struct locus                                    // Initialize the locus structure
{
    string regulatory;                          // regulatory alleles can only be the integers 0,1,2
    double coding;                              // coding alleles are continuous
    //    double xx;                                  // Trait 1 value
    //    double yy;                                  // Trait 2 value
};
struct phenotypes
{
    double xx;                                  // Trait 1's value
    double yy;                                  // Trait 2's value
    double ww;                                  // Individual's fitness
};

struct optima
{
    double x_opt;                               // Trait 1's optimum value
    double y_opt;                               // Trait 2's optimum value
    double om11;                                // variance in stabilizing selection
    double om12;                                // covariance in stabilizing selection
};

class Populations                            // Initialize the population class
{
    //private:                                  // Shouldn't we have some private members?
public:
    // Define the genotypic structure of the population.
    locus**** pop;
    
    // Initialize variables that will make numPops populations
    int numPops,numInd,numLoci,numAlleles;
    
    // Start with a NULL population
    Populations();
    
    // Initialize population
    void initilizePop(string reg_pattern, double theta, double gamma, char mod, double a1, double a2, double allelic_Stdev);
    
    // Define the structure of the Phenotypes
    phenotypes** xys;
    
    // Initialize Phenotypes
    void initilizeXYs();
    
    // Get Phenotypes
    void getPheno(double theta, double gamma, char mod);
    void geno_to_pheno(double a11, double a12, double a21, double a22, string r00, string r01, string r10, string r11, double theta, double gamma, char mod, double &XX, double &YY);
    
    // Get Fitness
    void getFitness();
    void pheno_to_fitness(double xx, double yy, double xopt, double yopt, double om11, double om12, double &W);
    
    // Intialize Optima
    optima* opts;
    void initializeOPTIMA(double xxx, double yyy, double omm11, double omm12);
    
    
    // For cleaning up.
    void deletePops_XYs();
    
    // Destroy Populations.
    ~Populations() // destructor
    {
        // We need to deallocate our buffer
        deletePops_XYs();
    }
    
    // Print population to screen (for debugging)
    void printPop();
    
    // Calculate the fitness given the trait values, selection, and triat optimum
    
};

// Forward Declarations of Functions:       // Can be moved to a header file if gets too long.
void input(Populations *popPtr, int POPS, int INDS);             // These are function prototypes
void printLocus(locus Locus);               // These are function prototypes
void Pheno_to_Geno(string reg_pattern, double x, double y, double theta, double gamma, char *mod, double &a1, double &a2);
double make_genos(double geno_value, double allelic_stdev);
double getFitness(double xx, double yy, double xopt=300, double yopt=300, double om11=1000, double om12=500);



// Main Function to run:
int main(int argc, char *argv[])
{
    
    // Pull in command line arguments
//    if (argc != 15) {
//        // Inform the user of how to use the program if not entered in correctly
//        std::cout << "Usage is <num_pops> <num_individuals> <initial_x> <initial_y> <initial_reg_pattern> <theta> <gamma> <model> <num_generations> <allelic_stdev> <xopts> <yopts> <om11> <om12>\n";
//        std::cout << argc-1 << " given" << endl;
//        std::cout << "Be sure that the number of optima correspond to the number of populations" << endl;
//        exit(0);
//    }
//    

    
    int numPops(0);
    int numInds(1);
    double x1(1);
    double x2(1);
    char* reg_pattern;
    double theta(300);
    double gamma(1);
    char *mod;
    int num_generations(1);
    double allelic_Stdev(1);
    double XOPT(300);                // Trait 1 optimum
    double YOPT(300);                // Trait 2 optimum
    double OM11(10000);                // Variance in stabilizing selection (for all pops)
    double OM12(0);
    
//
    cout << "Number of arguments provided = " << (argc-1)/2 << endl;
    for(int i=1; i<argc; i++){
        if(string(argv[i]) == "-p"){
            cout << "# number of pops = " << argv[i+1] << endl;
            numPops = atoi(argv[i+1]);
        }
        if(string(argv[i]) == "-i"){
            cout << "# number of individuals = " << argv[i+1] << endl;
            numInds = atoi(argv[i+1]);
        }
        if(string(argv[i]) == "-x"){
            cout << "Initial mean trait x value = " << argv[i+1] << endl;
            x1 = atoi(argv[i+1]);
        }
        if(string(argv[i]) == "-y"){
            cout << "Initial mean trait y value = " << argv[i+1] << endl;
            x2 = atoi(argv[i+1]);
        }
        if(string(argv[i]) == "-n"){
            cout << "Initial network architecture = " << argv[i+1] << endl;
            reg_pattern = argv[i+1];
        }
        if(string(argv[i]) == "-t"){
            cout << "theta = " << argv[i+1] << endl;
            theta = atoi(argv[i+1]);
        }
        if(string(argv[i]) == "-g"){
            cout << "gamma = " << argv[i+1] << endl;
            theta = atoi(argv[i+1]);
        }
        if(string(argv[i]) == "-m"){
            cout << "Regulatory model = " << argv[i+1] << endl;
            mod = argv[i+1];
        }
        if(string(argv[i]) == "-G"){
            cout << "# of Generations = " << argv[i+1] << endl;
            num_generations = atoi(argv[i+1]);
        }
        if(string(argv[i]) == "-d"){
            cout << "Initial standard deviation of allelic values = " << argv[i+1] << endl;
            allelic_Stdev = atoi(argv[i+1]);
        }
        if(string(argv[i]) == "-X"){
            cout << "Triat 1's optimum = " << argv[i+1] << endl;
            XOPT = atoi(argv[i+1]);
        }
        if(string(argv[i]) == "-Y"){
            cout << "Trait 2's optimum = " << argv[i+1] << endl;
            YOPT = atoi(argv[i+1]);
        }
        if(string(argv[i]) == "-s"){
            cout << "Variance in stabilizing selection = " << argv[i+1] << endl;
            OM11 = atoi(argv[i+1]);
        }
        if(string(argv[i]) == "-S"){
            cout << "Covariance in stabilizing selection = " << argv[i+1] << endl;
            OM11 = atoi(argv[i+1]);
        }
        
    }
    
    
//    int numPops(atoi(argv[1]));                 // Number of Populations
//    int numInds(atoi(argv[3]));                 // Number of Individuals per population
//    double x1(atoi(argv[4]));                   // Initial mean trait value for trait # 1
//    double x2(atoi(argv[5]));                   // Initial mean trait value for trait # 2
//    char* reg_pattern = argv[6];                // r11, r12, r21, r22 where rij is the reg allele for locus i and allele j
//    double theta(atoi(argv[7]));                // Theta value for network regulation
//    double gamma(atoi(argv[8]));                // gamma (decay) rate for network regulation
//    char *mod = argv[9];                        // model used: 'A' or 'B'
//    int num_generations(atoi(argv[10]));         // Number of generations to simulate
//    double allelic_Stdev(atoi(argv[11]));       // Initial standard dev. of allelic values to seed standing genetic variation
//    double XOPT(atoi(argv[12]));                // Trait 1 optimum
//    double YOPT(atoi(argv[13]));                // Trait 2 optimum
//    double OM11(atoi(argv[14]));                // Variance in stabilizing selection (for all pops)
//    double OM12(atoi(argv[15]));                // Covariance in stabilizing selection (for all pops)
    
    
//    int numPops(atoi(argv[1]));                 // Number of Populations
//    int numInds(atoi(argv[2]));                 // Number of Individuals per population
//    double x1(atoi(argv[3]));                   // Initial mean trait value for trait # 1
//    double x2(atoi(argv[4]));                   // Initial mean trait value for trait # 2
//    char* reg_pattern = argv[5];                // r11, r12, r21, r22 where rij is the reg allele for locus i and allele j
//    double theta(atoi(argv[6]));                // Theta value for network regulation
//    double gamma(atoi(argv[7]));                // gamma (decay) rate for network regulation
//    char *mod = argv[8];                        // model used: 'A' or 'B'
//    int num_generations(atoi(argv[9]));         // Number of generations to simulate
//    double allelic_Stdev(atoi(argv[10]));       // Initial standard dev. of allelic values to seed standing genetic variation
//    double XOPT(atoi(argv[11]));                // Trait 1 optimum
//    double YOPT(atoi(argv[12]));                // Trait 2 optimum
//    double OM11(atoi(argv[13]));                // Variance in stabilizing selection (for all pops)
//    double OM12(atoi(argv[14]));                // Covariance in stabilizing selection (for all pops)
//    
    
    
    
    
    
    //
    // Running:
    Populations Pop;
    input(&Pop, numPops, numInds);
    
    // Initialize genotypic values. These will be updated in the next function.
    double a1, a2;
    
    // Find the genotypic values that make the starting two-trait phenotype
    Pheno_to_Geno(reg_pattern, x1, x2, theta, gamma, mod, a1, a2);
    
    // Initialize the populations by their...
    // ...Genotypes:
    Pop.initilizePop(reg_pattern, theta, gamma, *mod, a1, a2, allelic_Stdev); // TO DO: remove hard-coded variances in random normal generation
    
    // ...Phenotypes:
    Pop.initilizeXYs();
    Pop.getPheno(theta, gamma, *mod);
    
    Pop.initializeOPTIMA(XOPT, YOPT, OM11, OM12);
    //
    
    // Recursion:
    for(int g=0; g<num_generations; g++){
        
        // Fitness:
        Pop.getFitness();
        //            num_generations = num_generations + 1 -1;
                Pop.printPop();             // Troubleshooting: print populations
        // Need selection
        
        // Need mutation
        
        // Need mating
        
        // Need migration
        
        // Need to estimate hybrid fitness
        
        // Need output summary stats to file
        
        
    }
    
    
    // Cleaning up:
//    Pop.deletePops_XYs();
   
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
            pop[p][i]=new locus*[numLoci];
            int counter(0);
            for(int c=0;c<numLoci;c++){
                pop[p][i][c]=new locus[numAlleles];
                
                for(int l=0; l<numAlleles;l++){
                    // Need to make genotype for the individual
                    
                    if(c==0){
                        pop[p][i][c][l].coding=a1/numAlleles+make_genos(a1/numAlleles, allelic_Stdev);
                        pop[p][i][c][l].regulatory=static_cast<char>(reg_pattern[counter]);
                    }
                    if(c==1){
                        pop[p][i][c][l].coding=a2/numAlleles+make_genos(a2/numAlleles, allelic_Stdev);
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

void Populations::initilizeXYs()
{
    xys=new phenotypes*[numPops];
    for(int p=0; p<numPops; p++){
        xys[p]=new phenotypes[numInd];
        for(int i=0; i<numInd; i++){
            xys[p][i].xx = 0.0;
            xys[p][i].yy = 0.0;
            //            cout << xys[p][i].xx << endl;
        }
        
    }
}

void Populations::deletePops_XYs()
{
    for(int p=0; p<numPops; p++)
    {
        for(int i=0; i<numInd; i++)
        {
            for(int c=0; c<numLoci; c++)
            {
                delete[] pop[p][i][c];
            }
            delete pop[p][i];

        }
        delete pop[p];
        delete xys[p];
    }
    delete pop;
    delete xys;
    delete opts;
}


void Populations::initializeOPTIMA(double xxx, double yyy, double omm11, double omm12)
{
    opts=new optima[numPops];
    for(int p=0; p<numPops; p++){
        
        opts[p].x_opt = xxx;
        opts[p].y_opt = yyy;
        opts[p].om11 = omm11;
        opts[p].om12 = omm12;
    }
}

void Populations::getPheno(double theta, double gamma, char mod)
{
    
    for(int p=0; p<numPops; p++){
        for(int i=0; i<numInd; i++){
            
            /* There are two structs with the Populations class: locus and phenotypes.
             locus holds the geneotypes and phenotypes hold the phenotypes.
             This function calls the geno_to_pheno function which actually retrieves the phenotypes.
             
             */
            
            geno_to_pheno(pop[p][i][0][0].coding,
                          pop[p][i][0][1].coding,
                          pop[p][i][1][0].coding,
                          pop[p][i][1][1].coding,
                          pop[p][i][0][0].regulatory,
                          pop[p][i][0][1].regulatory,
                          pop[p][i][1][0].regulatory,
                          pop[p][i][1][1].regulatory,
                          theta,
                          gamma,
                          mod,
                          xys[p][i].xx,
                          xys[p][i].yy);
            
        }
        
    }
}

void Populations::geno_to_pheno(double a11, double a12, double a21, double a22, string r00, string r01, string r10, string r11, double theta, double gamma, char mod, double &XX, double &YY)
{
    //    cout << "model is: " << mod << endl;
    switch(mod)
    {
        case 'A':
            if(r00 == "0" && r01 == "0" && r10=="0" && r11=="0"){
                XX = (a11+a12-a21-a22-gamma*theta+sqrt(4*(a11+a12)*gamma*theta+(-a11-a12+a21+a22+gamma*theta)*exp(2)))/(2*gamma);
                YY = (-a11-a12+a21+a22-gamma*theta+sqrt(4*(a11+a12)*gamma*theta+(-a11-a12+a21+a22+gamma*theta)*exp(2)))/(2*gamma);
            }
            if(r00 == "0" && r01 == "0" && r10=="1" && r11=="0"){
                XX = ((a11+a12-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a11+a12)*(a21+gamma*theta)+(-a11-a12+a21+a22+gamma*theta)*exp(2))))/(2*(a21+gamma*theta));
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a11+a12)*(a21+gamma*theta)+(-a11-a12+a21+a22+gamma*theta)*exp(2))))/(2*gamma*theta);
            }
            if(r00 == "0" && r01 == "0" && r10=="2" && r11=="0"){
                XX = ((a11+a12-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((a11+a12-a22-gamma*theta)*exp(2)+4*(a11+a12)*(a21+gamma*theta))))/(2*(a21+gamma*theta));
                YY = (-(a11+a12-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((a11+a12-a22-gamma*theta)*exp(2)+4*(a11+a12)*(a21+gamma*theta))))/(2*gamma*theta);
            }
            if(r00 == "0" && r01 == "0" && r10=="0" && r11=="1"){
                XX = ((a11+a12-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a11+a12)*(a22+gamma*theta)+(-a11-a12+a21+a22+gamma*theta)*exp(2))))/(2*(a22+gamma*theta));
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a11+a12)*(a22+gamma*theta)+(-a11-a12+a21+a22+gamma*theta)*exp(2))))/(2*gamma*theta);
            }
            if(r00 == "0" && r01 == "0" && r10=="1" && r11=="1"){
                XX = ((a11+a12)*theta)/(a21+a22+gamma*theta);
                YY = (a21+a22)/gamma;
            }
            if(r00 == "0" && r01 == "0" && r10=="2" && r11=="1"){
                XX = (theta*(a11+a12-a22-gamma*theta)+sqrt(theta*exp(2)*(a11*exp(2)+a12*exp(2)+(a22+gamma*theta)*exp(2)+2*a12*(2*a21+a22+gamma*theta)+2*a11*(a12+2*a21+a22+gamma*theta))))/(2*(a21+a22+gamma*theta));
                YY = (-theta*(a11+a12-a22+gamma*theta)+sqrt(theta*exp(2)*(a11*exp(2)+a12*exp(2)+(a22+gamma*theta)*exp(2)+2*a12*(2*a21+a22+gamma*theta)+2*a11*(a12+2*a21+a22+gamma*theta))))/(2*gamma*theta);
            }
            if(r00 == "0" && r01 == "0" && r10=="0" && r11=="2"){
                XX = ((a11+a12-a21)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((a11+a12-a21-gamma*theta)*exp(2)+4*(a11+a12)*(a22+gamma*theta))))/(2*(a22+gamma*theta));
                YY = (-(a11+a12-a21)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((a11+a12-a21-gamma*theta)*exp(2)+4*(a11+a12)*(a22+gamma*theta))))/(2*gamma*theta);
            }
            if(r00 == "0" && r01 == "0" && r10=="1" && r11=="2"){
                XX = (theta*(a11+a12-a21-gamma*theta)+sqrt(theta*exp(2)*(a11*exp(2)+a12*exp(2)+(a21+gamma*theta)*exp(2)+2*a12*(a21+2*a22+gamma*theta)+2*a11*(a12+a21+2*a22+gamma*theta))))/(2*(a21+a22+gamma*theta));
                YY = (-theta*(a11+a12-a21+gamma*theta)+sqrt(theta*exp(2)*(a11*exp(2)+a12*exp(2)+(a21+gamma*theta)*exp(2)+2*a12*(a21+2*a22+gamma*theta)+2*a11*(a12+a21+2*a22+gamma*theta))))/(2*gamma*theta);
            }
            if(r00 == "0" && r01 == "0" && r10=="2" && r11=="2"){
                XX = ((a11+a12)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((a11+a12-gamma*theta)*exp(2)+4*(a11+a12)*(a21+a22+gamma*theta))))/(2*(a21+a22+gamma*theta));
                YY = (-theta*(a11+a12+gamma*theta)+sqrt(theta*exp(2)*((a11+a12)*(a11+a12+4*(a21+a22))+2*(a11+a12)*gamma*theta+gamma*exp(2)*theta*exp(2))))/(2*gamma*theta);
            }
            if(r00 == "1" && r01 == "0" && r10=="0" && r11=="0"){
                XX = ((a11+a12-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a21+a22)*(a11+gamma*theta)+(a11+a12-a21-a22+gamma*theta)*exp(2))))/(2*gamma*theta);
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a21+a22)*(a11+gamma*theta)+(a11+a12-a21-a22+gamma*theta)*exp(2))))/(2*(a11+gamma*theta));
            }
            if(r00 == "1" && r01 == "0" && r10=="1" && r11=="0"){
                XX = (1/(2*gamma*(a21+gamma*theta)))*(a11*(a21+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(a11*(a21+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+gamma*theta)))*(a11*(a21-gamma*theta)+gamma*theta*(-a12+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(a11*(a21+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)));
            }
            if(r00 == "1" && r01 == "0" && r10=="2" && r11=="0"){
                XX = (1/(2*gamma*(a21+gamma*theta)))*(a11*(a21+gamma*theta)-gamma*theta*(-a12+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a11*a22+(a11+a12)*gamma*theta)+(a11*(a21+gamma*theta)-gamma*theta*(-a12+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+gamma*theta)))*(a11*(a21-gamma*theta)-gamma*theta*(a12-a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a11*a22+(a11+a12)*gamma*theta)+(a11*(a21+gamma*theta)-gamma*theta*(-a12+a22+gamma*theta))*exp(2)));
            }
            if(r00 == "1" && r01 == "0" && r10=="0" && r11=="1"){
                XX = (1/(2*gamma*(a22+gamma*theta)))*(a11*(a22+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(a11*(a22+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+gamma*theta)))*(a11*(a22-gamma*theta)+gamma*theta*(-a12+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(a11*(a22+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)));
            }
            if(r00 == "1" && r01 == "0" && r10=="1" && r11=="1"){
                XX = (a11*(a21+a22)+(a11+a12)*gamma*theta)/(gamma*(a21+a22+gamma*theta));
                YY = (a21+a22)/gamma;
            }
            if(r00 == "1" && r01 == "0" && r10=="2" && r11=="1"){
                XX = (1/(2*gamma*(a21+a22+gamma*theta)))*(-gamma*theta*(-a12+a22+gamma*theta)+a11*(a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a11*a22+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a22+gamma*theta)-a11*(a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+gamma*theta)))*(a11*(a21+a22-gamma*theta)-gamma*theta*(a12-a22+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a11*a22+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a22+gamma*theta)-a11*(a21+a22+gamma*theta))*exp(2)));
            }
            if(r00 == "1" && r01 == "0" && r10=="0" && r11=="2"){
                XX = (1/(2*gamma*(a22+gamma*theta)))*(-gamma*theta*(-a12+a21+gamma*theta)+a11*(a22+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a11*a21+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+gamma*theta)-a11*(a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+gamma*theta)))*(a11*(a22-gamma*theta)-gamma*theta*(a12-a21+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a11*a21+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+gamma*theta)-a11*(a22+gamma*theta))*exp(2)));
            }
            if(r00 == "1" && r01 == "0" && r10=="1" && r11=="2"){
                XX = (1/(2*gamma*(a21+a22+gamma*theta)))*(-gamma*theta*(-a12+a21+gamma*theta)+a11*(a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a11*a21+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+gamma*theta)-a11*(a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+gamma*theta)))*(a11*(a21+a22-gamma*theta)-gamma*theta*(a12-a21+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a11*a21+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+gamma*theta)-a11*(a21+a22+gamma*theta))*exp(2)));
            }
            if(r00 == "1" && r01 == "0" && r10=="2" && r11=="2"){
                XX = (gamma*theta*(a12-gamma*theta)+a11*(a21+a22+gamma*theta)+sqrt(4*(a11+a12)*gamma*exp(2)*theta*exp(2)*(a21+a22+gamma*theta)+(gamma*theta*(a12-gamma*theta)+a11*(a21+a22+gamma*theta))*exp(2)))/(2*gamma*(a21+a22+gamma*theta));
                YY = (a11*(a21+a22-gamma*theta)-gamma*theta*(a12+gamma*theta)+sqrt(4*(a11+a12)*gamma*exp(2)*theta*exp(2)*(a21+a22+gamma*theta)+(gamma*theta*(a12-gamma*theta)+a11*(a21+a22+gamma*theta))*exp(2)))/(2*gamma*(a11+gamma*theta));
            }
            if(r00 == "2" && r01 == "0" && r10=="0" && r11=="0"){
                XX = (a12*theta-theta*(a21+a22+gamma*theta)+sqrt(theta*exp(2)*(a12*exp(2)+4*a11*(a21+a22)-2*a12*(a21+a22-gamma*theta)+(a21+a22+gamma*theta)*exp(2))))/(2*gamma*theta);
                YY = (-a12*theta+a21*theta+a22*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((-a12+a21+a22-gamma*theta)*exp(2)+4*(a21+a22)*(a11+gamma*theta))))/(2*(a11+gamma*theta));
            }
            if(r00 == "2" && r01 == "0" && r10=="1" && r11=="0"){
                XX = (a11*a21-gamma*theta*(-a12+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a11*(a21+a22)+a12*gamma*theta)+(a11*a21-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)))/(2*gamma*(a21+gamma*theta));
                YY = (a11*a21+gamma*theta*(-a12+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a11*(a21+a22)+a12*gamma*theta)+(a11*a21-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)))/(2*gamma*(a11+gamma*theta));
            }
            if(r00 == "2" && r01 == "0" && r10=="2" && r11=="0"){
                XX = (a11*a21-gamma*theta*(-a12+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a11*a22+a12*gamma*theta)+(a11*a21-gamma*theta*(-a12+a22+gamma*theta))*exp(2)))/(2*gamma*(a21+gamma*theta));
                YY = (a11*a21-gamma*theta*(a12-a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a11*a22+a12*gamma*theta)+(a11*a21-gamma*theta*(-a12+a22+gamma*theta))*exp(2)))/(2*gamma*(a11+gamma*theta));
            }
            if(r00 == "2" && r01 == "0" && r10=="0" && r11=="1"){
                XX = (a11*a22-gamma*theta*(-a12+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a11*(a21+a22)+a12*gamma*theta)+(a11*a22-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)))/(2*gamma*(a22+gamma*theta));
                YY = (a11*a22+gamma*theta*(-a12+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a11*(a21+a22)+a12*gamma*theta)+(a11*a22-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)))/(2*gamma*(a11+gamma*theta));
            }
            if(r00 == "2" && r01 == "0" && r10=="1" && r11=="1"){
                XX = (a11*(a21+a22)+a12*gamma*theta)/(gamma*(a21+a22+gamma*theta));
                YY = (a21+a22)/gamma;
            }
            if(r00 == "2" && r01 == "0" && r10=="2" && r11=="1"){
                XX = (a11*(a21+a22)-gamma*theta*(-a12+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a11*a22+a12*gamma*theta)+(a11*(a21+a22)-gamma*theta*(-a12+a22+gamma*theta))*exp(2)))/(2*gamma*(a21+a22+gamma*theta));
                YY = (a11*(a21+a22)-gamma*theta*(a12-a22+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a11*a22+a12*gamma*theta)+(a11*(a21+a22)-gamma*theta*(-a12+a22+gamma*theta))*exp(2)))/(2*gamma*(a11+gamma*theta));
            }
            if(r00 == "2" && r01 == "0" && r10=="0" && r11=="2"){
                XX = (a11*a22-gamma*theta*(-a12+a21+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a11*a21+a12*gamma*theta)+(a11*a22-gamma*theta*(-a12+a21+gamma*theta))*exp(2)))/(2*gamma*(a22+gamma*theta));
                YY = (a11*a22-gamma*theta*(a12-a21+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a11*a21+a12*gamma*theta)+(a11*a22-gamma*theta*(-a12+a21+gamma*theta))*exp(2)))/(2*gamma*(a11+gamma*theta));
            }
            if(r00 == "2" && r01 == "0" && r10=="1" && r11=="2"){
                XX = (a11*(a21+a22)-gamma*theta*(-a12+a21+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a11*a21+a12*gamma*theta)+(a11*(a21+a22)-gamma*theta*(-a12+a21+gamma*theta))*exp(2)))/(2*gamma*(a21+a22+gamma*theta));
                YY = (a11*(a21+a22)-gamma*theta*(a12-a21+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a11*a21+a12*gamma*theta)+(a11*(a21+a22)-gamma*theta*(-a12+a21+gamma*theta))*exp(2)))/(2*gamma*(a11+gamma*theta));
            }
            if(r00 == "2" && r01 == "0" && r10=="2" && r11=="2"){
                XX = (a11*(a21+a22)+gamma*theta*(a12-gamma*theta)+sqrt(4*a12*gamma*exp(2)*theta*exp(2)*(a21+a22+gamma*theta)+(a11*(a21+a22)+gamma*theta*(a12-gamma*theta))*exp(2)))/(2*gamma*(a21+a22+gamma*theta));
                YY = (a11*(a21+a22)-gamma*theta*(a12+gamma*theta)+sqrt(4*a12*gamma*exp(2)*theta*exp(2)*(a21+a22+gamma*theta)+(a11*(a21+a22)+gamma*theta*(a12-gamma*theta))*exp(2)))/(2*gamma*(a11+gamma*theta));
            }
            if(r00 == "0" && r01 == "1" && r10=="0" && r11=="0"){
                XX = (a11*theta-theta*(-a12+a21+a22+gamma*theta)+sqrt(theta*exp(2)*(a11*exp(2)+2*a11*(a12-a21-a22+gamma*theta)+(a12+a21+a22+gamma*theta)*exp(2))))/(2*gamma*theta);
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a21+a22)*(a12+gamma*theta)+(a11+a12-a21-a22+gamma*theta)*exp(2))))/(2*(a12+gamma*theta));
            }
            if(r00 == "0" && r01 == "1" && r10=="1" && r11=="0"){
                XX = (1/(2*gamma*(a21+gamma*theta)))*(a12*(a21+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(a12*(a21+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a12+gamma*theta)))*(a12*(a21-gamma*theta)+gamma*theta*(-a11+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(a12*(a21+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)));
            }
            if(r00 == "0" && r01 == "1" && r10=="2" && r11=="0"){
                XX = (1/(2*gamma*(a21+gamma*theta)))*(a12*(a21+gamma*theta)-gamma*theta*(-a11+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a12*a22+(a11+a12)*gamma*theta)+(a12*(a21+gamma*theta)-gamma*theta*(-a11+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a12+gamma*theta)))*(a12*(a21-gamma*theta)-gamma*theta*(a11-a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a12*a22+(a11+a12)*gamma*theta)+(a12*(a21+gamma*theta)-gamma*theta*(-a11+a22+gamma*theta))*exp(2)));
            }
            if(r00 == "0" && r01 == "1" && r10=="0" && r11=="1"){
                XX = (1/(2*gamma*(a22+gamma*theta)))*(a12*(a22+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(a12*(a22+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a12+gamma*theta)))*(a12*(a22-gamma*theta)+gamma*theta*(-a11+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(a12*(a22+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)));
            }
            if(r00 == "0" && r01 == "1" && r10=="1" && r11=="1"){
                XX = (a12*(a21+a22)+(a11+a12)*gamma*theta)/(gamma*(a21+a22+gamma*theta));
                YY = (a21+a22)/gamma;
            }
            if(r00 == "0" && r01 == "1" && r10=="2" && r11=="1"){
                XX = (1/(2*gamma*(a21+a22+gamma*theta)))*(-gamma*theta*(-a11+a22+gamma*theta)+a12*(a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a12*a22+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a22+gamma*theta)-a12*(a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a12+gamma*theta)))*(a12*(a21+a22-gamma*theta)-gamma*theta*(a11-a22+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a12*a22+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a22+gamma*theta)-a12*(a21+a22+gamma*theta))*exp(2)));
            }
            if(r00 == "0" && r01 == "1" && r10=="0" && r11=="2"){
                XX = (1/(2*gamma*(a22+gamma*theta)))*(-gamma*theta*(-a11+a21+gamma*theta)+a12*(a22+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a12*a21+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+gamma*theta)-a12*(a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a12+gamma*theta)))*(a12*(a22-gamma*theta)-gamma*theta*(a11-a21+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a12*a21+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+gamma*theta)-a12*(a22+gamma*theta))*exp(2)));
            }
            if(r00 == "0" && r01 == "1" && r10=="1" && r11=="2"){
                XX = (1/(2*gamma*(a21+a22+gamma*theta)))*(-gamma*theta*(-a11+a21+gamma*theta)+a12*(a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a12*a21+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+gamma*theta)-a12*(a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a12+gamma*theta)))*(a12*(a21+a22-gamma*theta)-gamma*theta*(a11-a21+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a12*a21+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+gamma*theta)-a12*(a21+a22+gamma*theta))*exp(2)));
            }
            if(r00 == "0" && r01 == "1" && r10=="2" && r11=="2"){
                XX = (gamma*theta*(a11-gamma*theta)+a12*(a21+a22+gamma*theta)+sqrt(4*(a11+a12)*gamma*exp(2)*theta*exp(2)*(a21+a22+gamma*theta)+(gamma*theta*(a11-gamma*theta)+a12*(a21+a22+gamma*theta))*exp(2)))/(2*gamma*(a21+a22+gamma*theta));
                YY = (a12*(a21+a22-gamma*theta)-gamma*theta*(a11+gamma*theta)+sqrt(4*(a11+a12)*gamma*exp(2)*theta*exp(2)*(a21+a22+gamma*theta)+(gamma*theta*(a11-gamma*theta)+a12*(a21+a22+gamma*theta))*exp(2)))/(2*gamma*(a12+gamma*theta));
            }
            if(r00 == "1" && r01 == "1" && r10=="0" && r11=="0"){
                XX = (a11+a12)/gamma;
                YY = ((a21+a22)*theta)/(a11+a12+gamma*theta);
            }
            if(r00 == "1" && r01 == "1" && r10=="1" && r11=="0"){
                XX = (a11+a12)/gamma;
                YY = ((a11+a12)*a21+(a21+a22)*gamma*theta)/(gamma*(a11+a12+gamma*theta));
            }
            if(r00 == "1" && r01 == "1" && r10=="2" && r11=="0"){
                XX = (a11+a12)/gamma;
                YY = ((a11+a12)*a21+a22*gamma*theta)/(gamma*(a11+a12+gamma*theta));
            }
            if(r00 == "1" && r01 == "1" && r10=="0" && r11=="1"){
                XX = (a11+a12)/gamma;
                YY = ((a11+a12)*a22+(a21+a22)*gamma*theta)/(gamma*(a11+a12+gamma*theta));
            }
            if(r00 == "1" && r01 == "1" && r10=="1" && r11=="1"){
                XX = (a11+a12)/gamma;
                YY = (a21+a22)/gamma;
            }
            if(r00 == "1" && r01 == "1" && r10=="2" && r11=="1"){
                XX = (a11+a12)/gamma;
                YY = ((a11+a12)*(a21+a22)+a22*gamma*theta)/(gamma*(a11+a12+gamma*theta));
            }
            if(r00 == "1" && r01 == "1" && r10=="0" && r11=="2"){
                XX = (a11+a12)/gamma;
                YY = ((a11+a12)*a22+a21*gamma*theta)/(gamma*(a11+a12+gamma*theta));
            }
            if(r00 == "1" && r01 == "1" && r10=="1" && r11=="2"){
                XX = (a11+a12)/gamma;
                YY = ((a11+a12)*(a21+a22)+a21*gamma*theta)/(gamma*(a11+a12+gamma*theta));
            }
            if(r00 == "1" && r01 == "1" && r10=="2" && r11=="2"){
                XX = (a11+a12)/gamma;
                YY = ((a11+a12)*(a21+a22))/(gamma*(a11+a12+gamma*theta));
            }
            if(r00 == "2" && r01 == "1" && r10=="0" && r11=="0"){
                XX = (a12*theta-theta*(a21+a22+gamma*theta)+sqrt(theta*exp(2)*(4*a11*(a21+a22)+(a12+a21+a22+gamma*theta)*exp(2))))/(2*gamma*theta);
                YY = (theta*(-a12+a21+a22-gamma*theta)+sqrt(theta*exp(2)*(4*a11*(a21+a22)+(a12+a21+a22+gamma*theta)*exp(2))))/(2*(a11+a12+gamma*theta));
            }
            if(r00 == "2" && r01 == "1" && r10=="1" && r11=="0"){
                XX = (1/(2*gamma*(a21+gamma*theta)))*(a11*a21+a12*(a21+gamma*theta)-gamma*theta*(a21+a22+gamma*theta)+sqrt(a11*exp(2)*a21*exp(2)+(a12*(a21+gamma*theta)+gamma*theta*(a21+a22+gamma*theta))*exp(2)+2*a11*(a12*a21*(a21+gamma*theta)+gamma*theta*(a21*(a21+a22)+(a21+2*a22)*gamma*theta))));
                YY = (1/(2*gamma*(a11+a12+gamma*theta)))*(a11*a21+a12*a21-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+a12+gamma*theta)*(a12*a21+(a21+a22)*gamma*theta)+(a11*a21+a12*(a21-gamma*theta)+gamma*theta*(a21+a22-gamma*theta))*exp(2)));
            }
            if(r00 == "2" && r01 == "1" && r10=="2" && r11=="0"){
                XX = (1/(2*gamma*(a21+gamma*theta)))*(a11*a21+a12*a21+a12*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+a12+gamma*theta)*(a12*a21+a22*gamma*theta)+(a11*a21+a12*(a21-gamma*theta)+gamma*theta*(a22-gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+a12+gamma*theta)))*(a11*a21+a12*a21-a12*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+a12+gamma*theta)*(a12*a21+a22*gamma*theta)+(a11*a21+a12*(a21-gamma*theta)+gamma*theta*(a22-gamma*theta))*exp(2)));
            }
            if(r00 == "2" && r01 == "1" && r10=="0" && r11=="1"){
                XX = (1/(2*gamma*(a22+gamma*theta)))*(a11*a22+a12*(a22+gamma*theta)-gamma*theta*(a21+a22+gamma*theta)+sqrt(a11*exp(2)*a22*exp(2)+(a12*(a22+gamma*theta)+gamma*theta*(a21+a22+gamma*theta))*exp(2)+2*a11*(a12*a22*(a22+gamma*theta)+gamma*theta*(a22*(a21+a22)+(2*a21+a22)*gamma*theta))));
                YY = (1/(2*gamma*(a11+a12+gamma*theta)))*(a11*a22+a12*a22-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+a12+gamma*theta)*(a12*a22+(a21+a22)*gamma*theta)+(a11*a22+a12*(a22-gamma*theta)+gamma*theta*(a21+a22-gamma*theta))*exp(2)));
            }
            if(r00 == "2" && r01 == "1" && r10=="1" && r11=="1"){
                XX = ((a11+a12)*(a21+a22)+a12*gamma*theta)/(gamma*(a21+a22+gamma*theta));
                YY = (a21+a22)/gamma;
            }
            if(r00 == "2" && r01 == "1" && r10=="2" && r11=="1"){
                XX = (1/(2*gamma*(a21+a22+gamma*theta)))*(a11*a21+a11*a22-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+a12*(a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*((a11+a12)*a22+a12*gamma*theta)+(a11*(a21+a22)-gamma*theta*(a22+gamma*theta)+a12*(a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+a12+gamma*theta)))*(a11*a21+a11*a22+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+a12*(a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*((a11+a12)*a22+a12*gamma*theta)+(a11*(a21+a22)-gamma*theta*(a22+gamma*theta)+a12*(a21+a22+gamma*theta))*exp(2)));
            }
            if(r00 == "2" && r01 == "1" && r10=="0" && r11=="2"){
                XX = (1/(2*gamma*(a22+gamma*theta)))*(a11*a22+a12*a22+a12*gamma*theta-a21*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+a12+gamma*theta)*(a12*a22+a21*gamma*theta)+(a11*a22+gamma*theta*(a21-gamma*theta)+a12*(a22-gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+a12+gamma*theta)))*(a11*a22+a12*a22-a12*gamma*theta+a21*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+a12+gamma*theta)*(a12*a22+a21*gamma*theta)+(a11*a22+gamma*theta*(a21-gamma*theta)+a12*(a22-gamma*theta))*exp(2)));
            }
            if(r00 == "2" && r01 == "1" && r10=="1" && r11=="2"){
                XX = (1/(2*gamma*(a21+a22+gamma*theta)))*(a11*a21+a11*a22-a21*gamma*theta-gamma*exp(2)*theta*exp(2)+a12*(a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*((a11+a12)*a21+a12*gamma*theta)+(a11*(a21+a22)-gamma*theta*(a21+gamma*theta)+a12*(a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+a12+gamma*theta)))*(a11*a21+a11*a22+a21*gamma*theta-gamma*exp(2)*theta*exp(2)+a12*(a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*((a11+a12)*a21+a12*gamma*theta)+(a11*(a21+a22)-gamma*theta*(a21+gamma*theta)+a12*(a21+a22+gamma*theta))*exp(2)));
            }
            if(r00 == "2" && r01 == "1" && r10=="2" && r11=="2"){
                XX = (a11*a21+a11*a22-gamma*exp(2)*theta*exp(2)+a12*(a21+a22+gamma*theta)+sqrt(4*a12*gamma*exp(2)*theta*exp(2)*(a21+a22+gamma*theta)+((a11+a12)*(a21+a22)+a12*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)))/(2*gamma*(a21+a22+gamma*theta));
                YY = (a11*a21+a11*a22-gamma*exp(2)*theta*exp(2)+a12*(a21+a22-gamma*theta)+sqrt(4*a12*gamma*exp(2)*theta*exp(2)*(a21+a22+gamma*theta)+((a11+a12)*(a21+a22)+a12*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)))/(2*gamma*(a11+a12+gamma*theta));
            }
            if(r00 == "0" && r01 == "2" && r10=="0" && r11=="0"){
                XX = ((a11-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((-a11+a21+a22-gamma*theta)*exp(2)+4*(a21+a22)*(a12+gamma*theta))))/(2*gamma*theta);
                YY = (-a11*theta+a21*theta+a22*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((-a11+a21+a22-gamma*theta)*exp(2)+4*(a21+a22)*(a12+gamma*theta))))/(2*(a12+gamma*theta));
            }
            if(r00 == "0" && r01 == "2" && r10=="1" && r11=="0"){
                XX = (a12*a21-gamma*theta*(-a11+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a12*(a21+a22)+a11*gamma*theta)+(a12*a21-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)))/(2*gamma*(a21+gamma*theta));
                YY = (a12*a21+gamma*theta*(-a11+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a12*(a21+a22)+a11*gamma*theta)+(a12*a21-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)))/(2*gamma*(a12+gamma*theta));
            }
            if(r00 == "0" && r01 == "2" && r10=="2" && r11=="0"){
                XX = (a12*a21-gamma*theta*(-a11+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a12*a22+a11*gamma*theta)+(a12*a21-gamma*theta*(-a11+a22+gamma*theta))*exp(2)))/(2*gamma*(a21+gamma*theta));
                YY = (a12*a21-gamma*theta*(a11-a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a12*a22+a11*gamma*theta)+(a12*a21-gamma*theta*(-a11+a22+gamma*theta))*exp(2)))/(2*gamma*(a12+gamma*theta));
            }
            if(r00 == "0" && r01 == "2" && r10=="0" && r11=="1"){
                XX = (a12*a22-gamma*theta*(-a11+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a12*(a21+a22)+a11*gamma*theta)+(a12*a22-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)))/(2*gamma*(a22+gamma*theta));
                YY = (a12*a22+gamma*theta*(-a11+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a12*(a21+a22)+a11*gamma*theta)+(a12*a22-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)))/(2*gamma*(a12+gamma*theta));
            }
            if(r00 == "0" && r01 == "2" && r10=="1" && r11=="1"){
                XX = (a12*(a21+a22)+a11*gamma*theta)/(gamma*(a21+a22+gamma*theta));
                YY = (a21+a22)/gamma;
            }
            if(r00 == "0" && r01 == "2" && r10=="2" && r11=="1"){
                XX = (a12*(a21+a22)-gamma*theta*(-a11+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a12*a22+a11*gamma*theta)+(a12*(a21+a22)-gamma*theta*(-a11+a22+gamma*theta))*exp(2)))/(2*gamma*(a21+a22+gamma*theta));
                YY = (a12*(a21+a22)-gamma*theta*(a11-a22+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a12*a22+a11*gamma*theta)+(a12*(a21+a22)-gamma*theta*(-a11+a22+gamma*theta))*exp(2)))/(2*gamma*(a12+gamma*theta));
            }
            if(r00 == "0" && r01 == "2" && r10=="0" && r11=="2"){
                XX = (a12*a22-gamma*theta*(-a11+a21+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a12*a21+a11*gamma*theta)+(a12*a22-gamma*theta*(-a11+a21+gamma*theta))*exp(2)))/(2*gamma*(a22+gamma*theta));
                YY = (a12*a22-gamma*theta*(a11-a21+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a12*a21+a11*gamma*theta)+(a12*a22-gamma*theta*(-a11+a21+gamma*theta))*exp(2)))/(2*gamma*(a12+gamma*theta));
            }
            if(r00 == "0" && r01 == "2" && r10=="1" && r11=="2"){
                XX = (a12*(a21+a22)-gamma*theta*(-a11+a21+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a12*a21+a11*gamma*theta)+(a12*(a21+a22)-gamma*theta*(-a11+a21+gamma*theta))*exp(2)))/(2*gamma*(a21+a22+gamma*theta));
                YY = (a12*(a21+a22)-gamma*theta*(a11-a21+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*(a12*a21+a11*gamma*theta)+(a12*(a21+a22)-gamma*theta*(-a11+a21+gamma*theta))*exp(2)))/(2*gamma*(a12+gamma*theta));
            }
            if(r00 == "0" && r01 == "2" && r10=="2" && r11=="2"){
                XX = (a12*(a21+a22)+gamma*theta*(a11-gamma*theta)+sqrt(4*a11*gamma*exp(2)*theta*exp(2)*(a21+a22+gamma*theta)+(a12*(a21+a22)+gamma*theta*(a11-gamma*theta))*exp(2)))/(2*gamma*(a21+a22+gamma*theta));
                YY = (a12*(a21+a22)-gamma*theta*(a11+gamma*theta)+sqrt(4*a11*gamma*exp(2)*theta*exp(2)*(a21+a22+gamma*theta)+(a12*(a21+a22)+gamma*theta*(a11-gamma*theta))*exp(2)))/(2*gamma*(a12+gamma*theta));
            }
            if(r00 == "1" && r01 == "2" && r10=="0" && r11=="0"){
                XX = (-theta*(-a11+a21+a22+gamma*theta)+sqrt(theta*exp(2)*(a11*exp(2)+4*a12*(a21+a22)+2*a11*(a21+a22+gamma*theta)+(a21+a22+gamma*theta)*exp(2))))/(2*gamma*theta);
                YY = (theta*(-a11+a21+a22-gamma*theta)+sqrt(theta*exp(2)*(a11*exp(2)+4*a12*(a21+a22)+2*a11*(a21+a22+gamma*theta)+(a21+a22+gamma*theta)*exp(2))))/(2*(a11+a12+gamma*theta));
            }
            if(r00 == "1" && r01 == "2" && r10=="1" && r11=="0"){
                XX = (1/(2*gamma*(a21+gamma*theta)))*(a11*a21+a12*a21+a11*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+a12+gamma*theta)*(a11*a21+(a21+a22)*gamma*theta)+(a12*a21+a11*(a21-gamma*theta)+gamma*theta*(a21+a22-gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+a12+gamma*theta)))*(a11*a21+a12*a21-a11*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+a12+gamma*theta)*(a11*a21+(a21+a22)*gamma*theta)+(a12*a21+a11*(a21-gamma*theta)+gamma*theta*(a21+a22-gamma*theta))*exp(2)));
            }
            if(r00 == "1" && r01 == "2" && r10=="2" && r11=="0"){
                XX = (1/(2*gamma*(a21+gamma*theta)))*(a11*a21+a12*a21+a11*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+a12+gamma*theta)*(a11*a21+a22*gamma*theta)+(a12*a21+a11*(a21-gamma*theta)+gamma*theta*(a22-gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+a12+gamma*theta)))*(a11*a21+a12*a21-a11*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+a12+gamma*theta)*(a11*a21+a22*gamma*theta)+(a12*a21+a11*(a21-gamma*theta)+gamma*theta*(a22-gamma*theta))*exp(2)));
            }
            if(r00 == "1" && r01 == "2" && r10=="0" && r11=="1"){
                XX = (1/(2*gamma*(a22+gamma*theta)))*(a11*a22+a12*a22+a11*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+a12+gamma*theta)*(a11*a22+(a21+a22)*gamma*theta)+(a12*a22+a11*(a22-gamma*theta)+gamma*theta*(a21+a22-gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+a12+gamma*theta)))*(a11*a22+a12*a22-a11*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+a12+gamma*theta)*(a11*a22+(a21+a22)*gamma*theta)+(a12*a22+a11*(a22-gamma*theta)+gamma*theta*(a21+a22-gamma*theta))*exp(2)));
            }
            if(r00 == "1" && r01 == "2" && r10=="1" && r11=="1"){
                XX = ((a11+a12)*(a21+a22)+a11*gamma*theta)/(gamma*(a21+a22+gamma*theta));
                YY = (a21+a22)/gamma;
            }
            if(r00 == "1" && r01 == "2" && r10=="2" && r11=="1"){
                XX = (1/(2*gamma*(a21+a22+gamma*theta)))*(a12*a21+a12*a22-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+a11*(a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*((a11+a12)*a22+a11*gamma*theta)+(a12*(a21+a22)-gamma*theta*(a22+gamma*theta)+a11*(a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+a12+gamma*theta)))*(a12*a21+a12*a22+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+a11*(a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*((a11+a12)*a22+a11*gamma*theta)+(a12*(a21+a22)-gamma*theta*(a22+gamma*theta)+a11*(a21+a22+gamma*theta))*exp(2)));
            }
            if(r00 == "1" && r01 == "2" && r10=="0" && r11=="2"){
                XX = (1/(2*gamma*(a22+gamma*theta)))*(a11*a22+a12*a22+a11*gamma*theta-a21*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+a12+gamma*theta)*(a11*a22+a21*gamma*theta)+(a12*a22+gamma*theta*(a21-gamma*theta)+a11*(a22-gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+a12+gamma*theta)))*(a11*a22+a12*a22-a11*gamma*theta+a21*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+a12+gamma*theta)*(a11*a22+a21*gamma*theta)+(a12*a22+gamma*theta*(a21-gamma*theta)+a11*(a22-gamma*theta))*exp(2)));
            }
            if(r00 == "1" && r01 == "2" && r10=="1" && r11=="2"){
                XX = (1/(2*gamma*(a21+a22+gamma*theta)))*(a12*a21+a12*a22-a21*gamma*theta-gamma*exp(2)*theta*exp(2)+a11*(a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*((a11+a12)*a21+a11*gamma*theta)+(a12*(a21+a22)-gamma*theta*(a21+gamma*theta)+a11*(a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+a12+gamma*theta)))*(a12*a21+a12*a22+a21*gamma*theta-gamma*exp(2)*theta*exp(2)+a11*(a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a21+a22+gamma*theta)*((a11+a12)*a21+a11*gamma*theta)+(a12*(a21+a22)-gamma*theta*(a21+gamma*theta)+a11*(a21+a22+gamma*theta))*exp(2)));
            }
            if(r00 == "1" && r01 == "2" && r10=="2" && r11=="2"){
                XX = (a12*a21+a12*a22-gamma*exp(2)*theta*exp(2)+a11*(a21+a22+gamma*theta)+sqrt(4*a11*gamma*exp(2)*theta*exp(2)*(a21+a22+gamma*theta)+((a11+a12)*(a21+a22)+a11*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)))/(2*gamma*(a21+a22+gamma*theta));
                YY = (a12*a21+a12*a22-gamma*exp(2)*theta*exp(2)+a11*(a21+a22-gamma*theta)+sqrt(4*a11*gamma*exp(2)*theta*exp(2)*(a21+a22+gamma*theta)+((a11+a12)*(a21+a22)+a11*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)))/(2*gamma*(a11+a12+gamma*theta));
            }
            if(r00 == "2" && r01 == "2" && r10=="0" && r11=="0"){
                XX = (-theta*(a21+a22+gamma*theta)+sqrt(theta*exp(2)*(4*a11*(a21+a22)+4*a12*(a21+a22)+(a21+a22+gamma*theta)*exp(2))))/(2*gamma*theta);
                YY = (theta*(a21+a22-gamma*theta)+sqrt(theta*exp(2)*(4*a11*(a21+a22)+4*a12*(a21+a22)+(a21+a22+gamma*theta)*exp(2))))/(2*(a11+a12+gamma*theta));
            }
            if(r00 == "2" && r01 == "2" && r10=="1" && r11=="0"){
                XX = (a11*a21+a12*a21-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*(a21+a22)*gamma*exp(2)*theta*exp(2)*(a11+a12+gamma*theta)+(a11*a21+a12*a21+gamma*theta*(a21+a22-gamma*theta))*exp(2)))/(2*gamma*(a21+gamma*theta));
                YY = (a11*a21+a12*a21+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*(a21+a22)*gamma*exp(2)*theta*exp(2)*(a11+a12+gamma*theta)+(a11*a21+a12*a21+gamma*theta*(a21+a22-gamma*theta))*exp(2)))/(2*gamma*(a11+a12+gamma*theta));
            }
            if(r00 == "2" && r01 == "2" && r10=="2" && r11=="0"){
                XX = (a11*a21+a12*a21-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*a22*gamma*exp(2)*theta*exp(2)*(a11+a12+gamma*theta)+((a11+a12)*a21+gamma*theta*(a22-gamma*theta))*exp(2)))/(2*gamma*(a21+gamma*theta));
                YY = (a11*a21+a12*a21+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*a22*gamma*exp(2)*theta*exp(2)*(a11+a12+gamma*theta)+((a11+a12)*a21+gamma*theta*(a22-gamma*theta))*exp(2)))/(2*gamma*(a11+a12+gamma*theta));
            }
            if(r00 == "2" && r01 == "2" && r10=="0" && r11=="1"){
                XX = (a11*a22+a12*a22-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*(a21+a22)*gamma*exp(2)*theta*exp(2)*(a11+a12+gamma*theta)+(a11*a22+a12*a22+gamma*theta*(a21+a22-gamma*theta))*exp(2)))/(2*gamma*(a22+gamma*theta));
                YY = (a11*a22+a12*a22+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*(a21+a22)*gamma*exp(2)*theta*exp(2)*(a11+a12+gamma*theta)+(a11*a22+a12*a22+gamma*theta*(a21+a22-gamma*theta))*exp(2)))/(2*gamma*(a11+a12+gamma*theta));
            }
            if(r00 == "2" && r01 == "2" && r10=="1" && r11=="1"){
                XX = ((a11+a12)*(a21+a22))/(gamma*(a21+a22+gamma*theta));
                YY = (a21+a22)/gamma;
            }
            if(r00 == "2" && r01 == "2" && r10=="2" && r11=="1"){
                XX = (1/(2*gamma*(a21+a22+gamma*theta)))*(a12*a21+a12*a22+a11*(a21+a22)-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*(a11+a12)*a22*gamma*theta*(a21+a22+gamma*theta)+(a11*(a21+a22)+a12*(a21+a22)-gamma*theta*(a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+a12+gamma*theta)))*(a12*a21+a12*a22+a11*(a21+a22)+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*(a11+a12)*a22*gamma*theta*(a21+a22+gamma*theta)+(a11*(a21+a22)+a12*(a21+a22)-gamma*theta*(a22+gamma*theta))*exp(2)));
            }
            if(r00 == "2" && r01 == "2" && r10=="0" && r11=="2"){
                XX = (a11*a22+a12*a22-a21*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*a21*gamma*exp(2)*theta*exp(2)*(a11+a12+gamma*theta)+((a11+a12)*a22+gamma*theta*(a21-gamma*theta))*exp(2)))/(2*gamma*(a22+gamma*theta));
                YY = (a11*a22+a12*a22+a21*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*a21*gamma*exp(2)*theta*exp(2)*(a11+a12+gamma*theta)+((a11+a12)*a22+gamma*theta*(a21-gamma*theta))*exp(2)))/(2*gamma*(a11+a12+gamma*theta));
            }
            if(r00 == "2" && r01 == "2" && r10=="1" && r11=="2"){
                XX = (1/(2*gamma*(a21+a22+gamma*theta)))*(a12*a21+a12*a22+a11*(a21+a22)-a21*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*(a11+a12)*a21*gamma*theta*(a21+a22+gamma*theta)+(a11*(a21+a22)+a12*(a21+a22)-gamma*theta*(a21+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+a12+gamma*theta)))*(a12*a21+a12*a22+a11*(a21+a22)+a21*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*(a11+a12)*a21*gamma*theta*(a21+a22+gamma*theta)+(a11*(a21+a22)+a12*(a21+a22)-gamma*theta*(a21+gamma*theta))*exp(2)));
            }
            if(r00 == "2" && r01 == "2" && r10=="2" && r11=="2"){
                XX = ((a11+a12)*(a21+a22)-gamma*exp(2)*theta*exp(2))/(gamma*(a21+a22+gamma*theta));
                YY = ((a11+a12)*(a21+a22)-gamma*exp(2)*theta*exp(2))/(gamma*(a11+a12+gamma*theta));
            }
            break;
            
        case 'B':
            if(r00 == "0" && r01 == "0" && r10 == "0" && r11 == "0"){
                XX = (a11+a12-a21-a22-gamma*theta+sqrt(4*(a11+a12)*gamma*theta+(-a11-a12+a21+a22+gamma*theta)*exp(2)))/(2*gamma);
                YY = (-a11-a12+a21+a22-gamma*theta+sqrt(4*(a11+a12)*gamma*theta+(-a11-a12+a21+a22+gamma*theta)*exp(2)))/(2*gamma);
            }
            
            
            if(r00 == "0" && r01 == "0" && r10 == "1" && r11 == "0"){
                XX = ((a11+a12-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a11+a12)*(a21+gamma*theta)+(-a11-a12+a21+a22+gamma*theta)*exp(2))))/(2*(a21+gamma*theta));
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a11+a12)*(a21+gamma*theta)+(-a11-a12+a21+a22+gamma*theta)*exp(2))))/(2*gamma*theta);
            }
            
            
            if(r00 == "0" && r01 == "0" && r10 == "2" && r11 == "0"){
                XX = ((a11+a12-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a11+a12)*(2*a21+gamma*theta)+(-a11-a12+a21+a22+gamma*theta)*exp(2))))/(4*a21+2*gamma*theta);
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a11+a12)*(2*a21+gamma*theta)+(-a11-a12+a21+a22+gamma*theta)*exp(2))))/(2*gamma*theta);
            }
            
            
            if(r00 == "0" && r01 == "0" && r10 == "0" && r11 == "1"){
                XX = ((a11+a12-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a11+a12)*(a22+gamma*theta)+(-a11-a12+a21+a22+gamma*theta)*exp(2))))/(2*(a22+gamma*theta));
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a11+a12)*(a22+gamma*theta)+(-a11-a12+a21+a22+gamma*theta)*exp(2))))/(2*gamma*theta);
            }
            
            
            if(r00 == "0" && r01 == "0" && r10 == "1" && r11 == "1"){
                XX = ((a11+a12)*theta)/(a21+a22+gamma*theta);
                YY = (a21+a22)/gamma;
            }
            
            
            if(r00 == "0" && r01 == "0" && r10 == "2" && r11 == "1"){
                XX = ((a11+a12-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((-a11-a12+a21+a22+gamma*theta)*exp(2)+4*(a11+a12)*(2*a21+a22+gamma*theta))))/(2*(2*a21+a22+gamma*theta));
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((-a11-a12+a21+a22+gamma*theta)*exp(2)+4*(a11+a12)*(2*a21+a22+gamma*theta))))/(2*gamma*theta);
            }
            
            
            if(r00 == "0" && r01 == "0" && r10 == "0" && r11 == "2"){
                XX = ((a11+a12-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((-a11-a12+a21+a22+gamma*theta)*exp(2)+4*(a11+a12)*(2*a22+gamma*theta))))/(4*a22+2*gamma*theta);
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((-a11-a12+a21+a22+gamma*theta)*exp(2)+4*(a11+a12)*(2*a22+gamma*theta))))/(2*gamma*theta);
            }
            
            
            if(r00 == "0" && r01 == "0" && r10 == "1" && r11 == "2"){
                XX = ((a11+a12-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((-a11-a12+a21+a22+gamma*theta)*exp(2)+4*(a11+a12)*(a21+2*a22+gamma*theta))))/(2*(a21+2*a22+gamma*theta));
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((-a11-a12+a21+a22+gamma*theta)*exp(2)+4*(a11+a12)*(a21+2*a22+gamma*theta))))/(2*gamma*theta);
            }
            
            
            if(r00 == "0" && r01 == "0" && r10 == "2" && r11 == "2"){
                XX = ((a11+a12-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((-a11-a12+a21+a22+gamma*theta)*exp(2)+4*(a11+a12)*(2*(a21+a22)+gamma*theta))))/(4*(a21+a22)+2*gamma*theta);
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*((-a11-a12+a21+a22+gamma*theta)*exp(2)+4*(a11+a12)*(2*(a21+a22)+gamma*theta))))/(2*gamma*theta);
            }
            
            
            if(r00 == "1" && r01 == "0" && r10 == "0" && r11 == "0"){
                XX = ((a11+a12-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a21+a22)*(a11+gamma*theta)+(a11+a12-a21-a22+gamma*theta)*exp(2))))/(2*gamma*theta);
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a21+a22)*(a11+gamma*theta)+(a11+a12-a21-a22+gamma*theta)*exp(2))))/(2*(a11+gamma*theta));
            }
            
            
            if(r00 == "1" && r01 == "0" && r10 == "1" && r11 == "0"){
                XX = (1/(2*gamma*(a21+gamma*theta)))*(a11*(a21+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(a11*(a21+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+gamma*theta)))*(a11*(a21-gamma*theta)+gamma*theta*(-a12+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(a11*(a21+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "1" && r01 == "0" && r10 == "2" && r11 == "0"){
                XX = (1/(2*gamma*(2*a21+gamma*theta)))*(a11*(2*a21+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(2*a21+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(a11*(2*a21+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+gamma*theta)))*(a11*(2*a21-gamma*theta)+gamma*theta*(-a12+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(2*a21+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(a11*(2*a21+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "1" && r01 == "0" && r10 == "0" && r11 == "1"){
                XX = (1/(2*gamma*(a22+gamma*theta)))*(a11*(a22+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(a11*(a22+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+gamma*theta)))*(a11*(a22-gamma*theta)+gamma*theta*(-a12+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(a11*(a22+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "1" && r01 == "0" && r10 == "1" && r11 == "1"){
                XX = (a11*(a21+a22)+(a11+a12)*gamma*theta)/(gamma*(a21+a22+gamma*theta));
                YY = (a21+a22)/gamma;
            }
            
            
            if(r00 == "1" && r01 == "0" && r10 == "2" && r11 == "1"){
                XX = (1/(2*gamma*(2*a21+a22+gamma*theta)))*(-gamma*theta*(-a12+a21+a22+gamma*theta)+a11*(2*a21+a22+gamma*theta)+sqrt(4*gamma*theta*(2*a21+a22+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(2*a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+gamma*theta)))*(gamma*theta*(-a12+a21+a22-gamma*theta)+a11*(2*a21+a22-gamma*theta)+sqrt(4*gamma*theta*(2*a21+a22+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(2*a21+a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "1" && r01 == "0" && r10 == "0" && r11 == "2"){
                XX = (1/(2*gamma*(2*a22+gamma*theta)))*(-gamma*theta*(-a12+a21+a22+gamma*theta)+a11*(2*a22+gamma*theta)+sqrt(4*gamma*theta*(2*a22+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(2*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+gamma*theta)))*(gamma*theta*(-a12+a21+a22-gamma*theta)+a11*(2*a22-gamma*theta)+sqrt(4*gamma*theta*(2*a22+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(2*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "1" && r01 == "0" && r10 == "1" && r11 == "2"){
                XX = (1/(2*gamma*(a21+2*a22+gamma*theta)))*(-gamma*theta*(-a12+a21+a22+gamma*theta)+a11*(a21+2*a22+gamma*theta)+sqrt(4*gamma*theta*(a21+2*a22+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(a21+2*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+gamma*theta)))*(gamma*theta*(-a12+a21+a22-gamma*theta)+a11*(a21+2*a22-gamma*theta)+sqrt(4*gamma*theta*(a21+2*a22+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(a21+2*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "1" && r01 == "0" && r10 == "2" && r11 == "2"){
                XX = (1/(2*gamma*(2*(a21+a22)+gamma*theta)))*(2*a11*(a21+a22)+a11*gamma*theta-gamma*theta*(-a12+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(2*(a21+a22)+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(2*(a21+a22)+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+gamma*theta)))*(2*a11*(a21+a22)-a11*gamma*theta+gamma*theta*(-a12+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(2*(a21+a22)+gamma*theta)*(a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(2*(a21+a22)+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "0" && r10 == "0" && r11 == "0"){
                XX = ((a11+a12-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a21+a22)*(2*a11+gamma*theta)+(a11+a12-a21-a22+gamma*theta)*exp(2))))/(2*gamma*theta);
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a21+a22)*(2*a11+gamma*theta)+(a11+a12-a21-a22+gamma*theta)*exp(2))))/(4*a11+2*gamma*theta);
            }
            
            
            if(r00 == "2" && r01 == "0" && r10 == "1" && r11 == "0"){
                XX = (1/(2*gamma*(a21+gamma*theta)))*(a11*(2*a21+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(2*a11*(a21+a22)+(a11+a12)*gamma*theta)+(a11*(2*a21+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a11+gamma*theta)))*(a11*(2*a21-gamma*theta)+gamma*theta*(-a12+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(2*a11*(a21+a22)+(a11+a12)*gamma*theta)+(a11*(2*a21+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "0" && r10 == "2" && r11 == "0"){
                XX = (1/(2*gamma*(2*a21+gamma*theta)))*(a11*(4*a21+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(2*a21+gamma*theta)*(2*a11*(a21+a22)+(a11+a12)*gamma*theta)+(a11*(4*a21+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a11+gamma*theta)))*(a11*(4*a21-gamma*theta)+gamma*theta*(-a12+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(2*a21+gamma*theta)*(2*a11*(a21+a22)+(a11+a12)*gamma*theta)+(a11*(4*a21+gamma*theta)-gamma*theta*(-a12+a21+a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "0" && r10 == "0" && r11 == "1"){
                XX = (1/(2*gamma*(a22+gamma*theta)))*(-gamma*theta*(-a12+a21+a22+gamma*theta)+a11*(2*a22+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(2*a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(2*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a11+gamma*theta)))*(gamma*theta*(-a12+a21+a22-gamma*theta)+a11*(2*a22-gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(2*a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(2*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "0" && r10 == "1" && r11 == "1"){
                XX = (2*a11*(a21+a22)+(a11+a12)*gamma*theta)/(gamma*(a21+a22+gamma*theta));
                YY = (a21+a22)/gamma;
            }
            
            
            if(r00 == "2" && r01 == "0" && r10 == "2" && r11 == "1"){
                XX = (1/(2*gamma*(2*a21+a22+gamma*theta)))*(-gamma*theta*(-a12+a21+a22+gamma*theta)+a11*(4*a21+2*a22+gamma*theta)+sqrt(4*gamma*theta*(2*a21+a22+gamma*theta)*(2*a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(4*a21+2*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a11+gamma*theta)))*(gamma*theta*(-a12+a21+a22-gamma*theta)+a11*(4*a21+2*a22-gamma*theta)+sqrt(4*gamma*theta*(2*a21+a22+gamma*theta)*(2*a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(4*a21+2*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "0" && r10 == "0" && r11 == "2"){
                XX = (1/(2*gamma*(2*a22+gamma*theta)))*(-gamma*theta*(-a12+a21+a22+gamma*theta)+a11*(4*a22+gamma*theta)+sqrt(4*gamma*theta*(2*a22+gamma*theta)*(2*a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(4*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a11+gamma*theta)))*(gamma*theta*(-a12+a21+a22-gamma*theta)+a11*(4*a22-gamma*theta)+sqrt(4*gamma*theta*(2*a22+gamma*theta)*(2*a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(4*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "0" && r10 == "1" && r11 == "2"){
                XX = (1/(2*gamma*(a21+2*a22+gamma*theta)))*(-gamma*theta*(-a12+a21+a22+gamma*theta)+a11*(2*a21+4*a22+gamma*theta)+sqrt(4*gamma*theta*(a21+2*a22+gamma*theta)*(2*a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(2*a21+4*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a11+gamma*theta)))*(gamma*theta*(-a12+a21+a22-gamma*theta)+a11*(2*a21+4*a22-gamma*theta)+sqrt(4*gamma*theta*(a21+2*a22+gamma*theta)*(2*a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(2*a21+4*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "0" && r10 == "2" && r11 == "2"){
                XX = (1/(2*gamma*(2*(a21+a22)+gamma*theta)))*(4*a11*(a21+a22)+a11*gamma*theta-gamma*theta*(-a12+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(2*(a21+a22)+gamma*theta)*(2*a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(4*(a21+a22)+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a11+gamma*theta)))*(4*a11*(a21+a22)-a11*gamma*theta+gamma*theta*(-a12+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(2*(a21+a22)+gamma*theta)*(2*a11*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a12+a21+a22+gamma*theta)-a11*(4*(a21+a22)+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "0" && r01 == "1" && r10 == "0" && r11 == "0"){
                XX = (a11*theta-theta*(-a12+a21+a22+gamma*theta)+sqrt(theta*exp(2)*(a11*exp(2)+2*a11*(a12-a21-a22+gamma*theta)+(a12+a21+a22+gamma*theta)*exp(2))))/(2*gamma*theta);
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a21+a22)*(a12+gamma*theta)+(a11+a12-a21-a22+gamma*theta)*exp(2))))/(2*(a12+gamma*theta));
            }
            
            
            if(r00 == "0" && r01 == "1" && r10 == "1" && r11 == "0"){
                XX = (1/(2*gamma*(a21+gamma*theta)))*(a12*(a21+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(a12*(a21+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a12+gamma*theta)))*(a12*(a21-gamma*theta)+gamma*theta*(-a11+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(a12*(a21+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "0" && r01 == "1" && r10 == "2" && r11 == "0"){
                XX = (1/(2*gamma*(2*a21+gamma*theta)))*(a12*(2*a21+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(2*a21+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(a12*(2*a21+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a12+gamma*theta)))*(a12*(2*a21-gamma*theta)+gamma*theta*(-a11+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(2*a21+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(a12*(2*a21+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "0" && r01 == "1" && r10 == "0" && r11 == "1"){
                XX = (1/(2*gamma*(a22+gamma*theta)))*(a12*(a22+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(a12*(a22+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a12+gamma*theta)))*(a12*(a22-gamma*theta)+gamma*theta*(-a11+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(a12*(a22+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "0" && r01 == "1" && r10 == "1" && r11 == "1"){
                XX = (a12*(a21+a22)+(a11+a12)*gamma*theta)/(gamma*(a21+a22+gamma*theta));
                YY = (a21+a22)/gamma;
            }
            
            
            if(r00 == "0" && r01 == "1" && r10 == "2" && r11 == "1"){
                XX = (1/(2*gamma*(2*a21+a22+gamma*theta)))*(-gamma*theta*(-a11+a21+a22+gamma*theta)+a12*(2*a21+a22+gamma*theta)+sqrt(4*gamma*theta*(2*a21+a22+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(2*a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a12+gamma*theta)))*(gamma*theta*(-a11+a21+a22-gamma*theta)+a12*(2*a21+a22-gamma*theta)+sqrt(4*gamma*theta*(2*a21+a22+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(2*a21+a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "0" && r01 == "1" && r10 == "0" && r11 == "2"){
                XX = (1/(2*gamma*(2*a22+gamma*theta)))*(-gamma*theta*(-a11+a21+a22+gamma*theta)+a12*(2*a22+gamma*theta)+sqrt(4*gamma*theta*(2*a22+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(2*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a12+gamma*theta)))*(gamma*theta*(-a11+a21+a22-gamma*theta)+a12*(2*a22-gamma*theta)+sqrt(4*gamma*theta*(2*a22+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(2*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "0" && r01 == "1" && r10 == "1" && r11 == "2"){
                XX = (1/(2*gamma*(a21+2*a22+gamma*theta)))*(-gamma*theta*(-a11+a21+a22+gamma*theta)+a12*(a21+2*a22+gamma*theta)+sqrt(4*gamma*theta*(a21+2*a22+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(a21+2*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a12+gamma*theta)))*(gamma*theta*(-a11+a21+a22-gamma*theta)+a12*(a21+2*a22-gamma*theta)+sqrt(4*gamma*theta*(a21+2*a22+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(a21+2*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "0" && r01 == "1" && r10 == "2" && r11 == "2"){
                XX = (1/(2*gamma*(2*(a21+a22)+gamma*theta)))*(2*a12*(a21+a22)+a12*gamma*theta-gamma*theta*(-a11+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(2*(a21+a22)+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(2*(a21+a22)+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a12+gamma*theta)))*(2*a12*(a21+a22)-a12*gamma*theta+gamma*theta*(-a11+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(2*(a21+a22)+gamma*theta)*(a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(2*(a21+a22)+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "1" && r01 == "1" && r10 == "0" && r11 == "0"){
                XX = (a11+a12)/gamma;
                YY = ((a21+a22)*theta)/(a11+a12+gamma*theta);
            }
            
            
            if(r00 == "1" && r01 == "1" && r10 == "1" && r11 == "0"){
                XX = (a11+a12)/gamma;
                YY = ((a11+a12)*a21+(a21+a22)*gamma*theta)/(gamma*(a11+a12+gamma*theta));
            }
            
            
            if(r00 == "1" && r01 == "1" && r10 == "2" && r11 == "0"){
                XX = (a11+a12)/gamma;
                YY = (2*(a11+a12)*a21+(a21+a22)*gamma*theta)/(gamma*(a11+a12+gamma*theta));
            }
            
            
            if(r00 == "1" && r01 == "1" && r10 == "0" && r11 == "1"){
                XX = (a11+a12)/gamma;
                YY = ((a11+a12)*a22+(a21+a22)*gamma*theta)/(gamma*(a11+a12+gamma*theta));
            }
            
            
            if(r00 == "1" && r01 == "1" && r10 == "1" && r11 == "1"){
                XX = (a11+a12)/gamma;
                YY = (a21+a22)/gamma;
            }
            
            
            if(r00 == "1" && r01 == "1" && r10 == "2" && r11 == "1"){
                XX = (a11+a12)/gamma;
                YY = ((a11+a12)*(2*a21+a22)+(a21+a22)*gamma*theta)/(gamma*(a11+a12+gamma*theta));
            }
            
            
            if(r00 == "1" && r01 == "1" && r10 == "0" && r11 == "2"){
                XX = (a11+a12)/gamma;
                YY = (2*(a11+a12)*a22+(a21+a22)*gamma*theta)/(gamma*(a11+a12+gamma*theta));
            }
            
            
            if(r00 == "1" && r01 == "1" && r10 == "1" && r11 == "2"){
                XX = (a11+a12)/gamma;
                YY = ((a11+a12)*(a21+2*a22)+(a21+a22)*gamma*theta)/(gamma*(a11+a12+gamma*theta));
            }
            
            
            if(r00 == "1" && r01 == "1" && r10 == "2" && r11 == "2"){
                XX = (a11+a12)/gamma;
                YY = ((a21+a22)*(2*(a11+a12)+gamma*theta))/(gamma*(a11+a12+gamma*theta));
            }
            
            
            if(r00 == "2" && r01 == "1" && r10 == "0" && r11 == "0"){
                XX = (a11*theta-theta*(-a12+a21+a22+gamma*theta)+sqrt(theta*exp(2)*(a11*exp(2)+(a12+a21+a22+gamma*theta)*exp(2)+2*a11*(a12+3*(a21+a22)+gamma*theta))))/(2*gamma*theta);
                YY = (-theta*(a11+a12-a21-a22+gamma*theta)+sqrt(theta*exp(2)*(a11*exp(2)+(a12+a21+a22+gamma*theta)*exp(2)+2*a11*(a12+3*(a21+a22)+gamma*theta))))/(2*(2*a11+a12+gamma*theta));
            }
            
            
            if(r00 == "2" && r01 == "1" && r10 == "1" && r11 == "0"){
                XX = (1/(2*gamma*(a21+gamma*theta)))*(2*a11*a21+a12*a21+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*a11+a12+gamma*theta)*((a11+a12)*a21+(a21+a22)*gamma*theta)+((2*a11+a12)*a21+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
                YY = (1/(2*gamma*(2*a11+a12+gamma*theta)))*(2*a11*a21+a12*a21-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*a11+a12+gamma*theta)*((a11+a12)*a21+(a21+a22)*gamma*theta)+((2*a11+a12)*a21+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "1" && r10 == "2" && r11 == "0"){
                XX = (1/(2*gamma*(2*a21+gamma*theta)))*(4*a11*a21+2*a12*a21+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*a11+a12+gamma*theta)*(2*(a11+a12)*a21+(a21+a22)*gamma*theta)+(2*(2*a11+a12)*a21+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
                YY = (1/(2*gamma*(2*a11+a12+gamma*theta)))*(4*a11*a21+2*a12*a21-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*a11+a12+gamma*theta)*(2*(a11+a12)*a21+(a21+a22)*gamma*theta)+(2*(2*a11+a12)*a21+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "1" && r10 == "0" && r11 == "1"){
                XX = (1/(2*gamma*(a22+gamma*theta)))*(2*a11*a22+a12*a22+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*a11+a12+gamma*theta)*((a11+a12)*a22+(a21+a22)*gamma*theta)+((2*a11+a12)*a22+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
                YY = (1/(2*gamma*(2*a11+a12+gamma*theta)))*(2*a11*a22+a12*a22-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*a11+a12+gamma*theta)*((a11+a12)*a22+(a21+a22)*gamma*theta)+((2*a11+a12)*a22+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "1" && r10 == "1" && r11 == "1"){
                XX = ((2*a11+a12)*(a21+a22)+(a11+a12)*gamma*theta)/(gamma*(a21+a22+gamma*theta));
                YY = (a21+a22)/gamma;
            }
            
            
            if(r00 == "2" && r01 == "1" && r10 == "2" && r11 == "1"){
                XX = (1/(2*gamma*(2*a21+a22+gamma*theta)))*(4*a11*a21+2*a12*a21+2*a11*a22+a12*a22+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*a21+a22+gamma*theta)*((2*a11+a12)*(a21+a22)+(a11+a12)*gamma*theta)+(-gamma*theta*(a21+a22+gamma*theta)+a12*(2*a21+a22+gamma*theta)+a11*(4*a21+2*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a11+a12+gamma*theta)))*(4*a11*a21+2*a12*a21+2*a11*a22+a12*a22-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*a21+a22+gamma*theta)*((2*a11+a12)*(a21+a22)+(a11+a12)*gamma*theta)+(-gamma*theta*(a21+a22+gamma*theta)+a12*(2*a21+a22+gamma*theta)+a11*(4*a21+2*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "1" && r10 == "0" && r11 == "2"){
                XX = (1/(2*gamma*(2*a22+gamma*theta)))*(4*a11*a22+2*a12*a22+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*a11+a12+gamma*theta)*(2*(a11+a12)*a22+(a21+a22)*gamma*theta)+(2*(2*a11+a12)*a22+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
                YY = (1/(2*gamma*(2*a11+a12+gamma*theta)))*(4*a11*a22+2*a12*a22-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*a11+a12+gamma*theta)*(2*(a11+a12)*a22+(a21+a22)*gamma*theta)+(2*(2*a11+a12)*a22+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "1" && r10 == "1" && r11 == "2"){
                XX = (1/(2*gamma*(a21+2*a22+gamma*theta)))*(2*a11*a21+a12*a21+4*a11*a22+2*a12*a22+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a21+2*a22+gamma*theta)*((2*a11+a12)*(a21+a22)+(a11+a12)*gamma*theta)+(-gamma*theta*(a21+a22+gamma*theta)+a12*(a21+2*a22+gamma*theta)+a11*(2*a21+4*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a11+a12+gamma*theta)))*(2*a11*a21+a12*a21+4*a11*a22+2*a12*a22-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a21+2*a22+gamma*theta)*((2*a11+a12)*(a21+a22)+(a11+a12)*gamma*theta)+(-gamma*theta*(a21+a22+gamma*theta)+a12*(a21+2*a22+gamma*theta)+a11*(2*a21+4*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "1" && r10 == "2" && r11 == "2"){
                XX = (1/(2*gamma*(2*(a21+a22)+gamma*theta)))*(4*a11*a21+2*a12*a21+4*a11*a22+2*a12*a22+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*(a21+a22)+gamma*theta)*((2*a11+a12)*(a21+a22)+(a11+a12)*gamma*theta)+(4*a11*(a21+a22)+2*a12*(a21+a22)+a11*gamma*theta+a12*gamma*theta-gamma*theta*(a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a11+a12+gamma*theta)))*(4*a11*a21+2*a12*a21+4*a11*a22+2*a12*a22-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*(a21+a22)+gamma*theta)*((2*a11+a12)*(a21+a22)+(a11+a12)*gamma*theta)+(4*a11*(a21+a22)+2*a12*(a21+a22)+a11*gamma*theta+a12*gamma*theta-gamma*theta*(a21+a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "0" && r01 == "2" && r10 == "0" && r11 == "0"){
                XX = ((a11+a12-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a21+a22)*(2*a12+gamma*theta)+(a11+a12-a21-a22+gamma*theta)*exp(2))))/(2*gamma*theta);
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a21+a22)*(2*a12+gamma*theta)+(a11+a12-a21-a22+gamma*theta)*exp(2))))/(4*a12+2*gamma*theta);
            }
            
            
            if(r00 == "0" && r01 == "2" && r10 == "1" && r11 == "0"){
                XX = (1/(2*gamma*(a21+gamma*theta)))*(a12*(2*a21+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(2*a12*(a21+a22)+(a11+a12)*gamma*theta)+(a12*(2*a21+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a12+gamma*theta)))*(a12*(2*a21-gamma*theta)+gamma*theta*(-a11+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(a21+gamma*theta)*(2*a12*(a21+a22)+(a11+a12)*gamma*theta)+(a12*(2*a21+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "0" && r01 == "2" && r10 == "2" && r11 == "0"){
                XX = (1/(2*gamma*(2*a21+gamma*theta)))*(a12*(4*a21+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(2*a21+gamma*theta)*(2*a12*(a21+a22)+(a11+a12)*gamma*theta)+(a12*(4*a21+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a12+gamma*theta)))*(a12*(4*a21-gamma*theta)+gamma*theta*(-a11+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(2*a21+gamma*theta)*(2*a12*(a21+a22)+(a11+a12)*gamma*theta)+(a12*(4*a21+gamma*theta)-gamma*theta*(-a11+a21+a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "0" && r01 == "2" && r10 == "0" && r11 == "1"){
                XX = (1/(2*gamma*(a22+gamma*theta)))*(-gamma*theta*(-a11+a21+a22+gamma*theta)+a12*(2*a22+gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(2*a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(2*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a12+gamma*theta)))*(gamma*theta*(-a11+a21+a22-gamma*theta)+a12*(2*a22-gamma*theta)+sqrt(4*gamma*theta*(a22+gamma*theta)*(2*a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(2*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "0" && r01 == "2" && r10 == "1" && r11 == "1"){
                XX = (2*a12*(a21+a22)+(a11+a12)*gamma*theta)/(gamma*(a21+a22+gamma*theta));
                YY = (a21+a22)/gamma;
            }
            
            
            if(r00 == "0" && r01 == "2" && r10 == "2" && r11 == "1"){
                XX = (1/(2*gamma*(2*a21+a22+gamma*theta)))*(-gamma*theta*(-a11+a21+a22+gamma*theta)+a12*(4*a21+2*a22+gamma*theta)+sqrt(4*gamma*theta*(2*a21+a22+gamma*theta)*(2*a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(4*a21+2*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a12+gamma*theta)))*(gamma*theta*(-a11+a21+a22-gamma*theta)+a12*(4*a21+2*a22-gamma*theta)+sqrt(4*gamma*theta*(2*a21+a22+gamma*theta)*(2*a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(4*a21+2*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "0" && r01 == "2" && r10 == "0" && r11 == "2"){
                XX = (1/(2*gamma*(2*a22+gamma*theta)))*(-gamma*theta*(-a11+a21+a22+gamma*theta)+a12*(4*a22+gamma*theta)+sqrt(4*gamma*theta*(2*a22+gamma*theta)*(2*a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(4*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a12+gamma*theta)))*(gamma*theta*(-a11+a21+a22-gamma*theta)+a12*(4*a22-gamma*theta)+sqrt(4*gamma*theta*(2*a22+gamma*theta)*(2*a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(4*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "0" && r01 == "2" && r10 == "1" && r11 == "2"){
                XX = (1/(2*gamma*(a21+2*a22+gamma*theta)))*(-gamma*theta*(-a11+a21+a22+gamma*theta)+a12*(2*a21+4*a22+gamma*theta)+sqrt(4*gamma*theta*(a21+2*a22+gamma*theta)*(2*a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(2*a21+4*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a12+gamma*theta)))*(gamma*theta*(-a11+a21+a22-gamma*theta)+a12*(2*a21+4*a22-gamma*theta)+sqrt(4*gamma*theta*(a21+2*a22+gamma*theta)*(2*a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(2*a21+4*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "0" && r01 == "2" && r10 == "2" && r11 == "2"){
                XX = (1/(2*gamma*(2*(a21+a22)+gamma*theta)))*(4*a12*(a21+a22)+a12*gamma*theta-gamma*theta*(-a11+a21+a22+gamma*theta)+sqrt(4*gamma*theta*(2*(a21+a22)+gamma*theta)*(2*a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(4*(a21+a22)+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*a12+gamma*theta)))*(4*a12*(a21+a22)-a12*gamma*theta+gamma*theta*(-a11+a21+a22-gamma*theta)+sqrt(4*gamma*theta*(2*(a21+a22)+gamma*theta)*(2*a12*(a21+a22)+(a11+a12)*gamma*theta)+(gamma*theta*(-a11+a21+a22+gamma*theta)-a12*(4*(a21+a22)+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "1" && r01 == "2" && r10 == "0" && r11 == "0"){
                XX = ((a11+a12-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a21+a22)*(a11+2*a12+gamma*theta)+(a11+a12-a21-a22+gamma*theta)*exp(2))))/(2*gamma*theta);
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a21+a22)*(a11+2*a12+gamma*theta)+(a11+a12-a21-a22+gamma*theta)*exp(2))))/(2*(a11+2*a12+gamma*theta));
            }
            
            
            if(r00 == "1" && r01 == "2" && r10 == "1" && r11 == "0"){
                XX = (1/(2*gamma*(a21+gamma*theta)))*(a11*a21+2*a12*a21+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+2*a12+gamma*theta)*((a11+a12)*a21+(a21+a22)*gamma*theta)+((a11+2*a12)*a21+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
                YY = (1/(2*gamma*(a11+2*a12+gamma*theta)))*(a11*a21+2*a12*a21-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+2*a12+gamma*theta)*((a11+a12)*a21+(a21+a22)*gamma*theta)+((a11+2*a12)*a21+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
            }
            
            
            if(r00 == "1" && r01 == "2" && r10 == "2" && r11 == "0"){
                XX = (1/(2*gamma*(2*a21+gamma*theta)))*(2*a11*a21+4*a12*a21+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+2*a12+gamma*theta)*(2*(a11+a12)*a21+(a21+a22)*gamma*theta)+(2*(a11+2*a12)*a21+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
                YY = (1/(2*gamma*(a11+2*a12+gamma*theta)))*(2*a11*a21+4*a12*a21-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+2*a12+gamma*theta)*(2*(a11+a12)*a21+(a21+a22)*gamma*theta)+(2*(a11+2*a12)*a21+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
            }
            
            
            if(r00 == "1" && r01 == "2" && r10 == "0" && r11 == "1"){
                XX = (1/(2*gamma*(a22+gamma*theta)))*(a11*a22+2*a12*a22+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+2*a12+gamma*theta)*((a11+a12)*a22+(a21+a22)*gamma*theta)+((a11+2*a12)*a22+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
                YY = (1/(2*gamma*(a11+2*a12+gamma*theta)))*(a11*a22+2*a12*a22-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+2*a12+gamma*theta)*((a11+a12)*a22+(a21+a22)*gamma*theta)+((a11+2*a12)*a22+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
            }
            
            
            if(r00 == "1" && r01 == "2" && r10 == "1" && r11 == "1"){
                XX = ((a11+2*a12)*(a21+a22)+(a11+a12)*gamma*theta)/(gamma*(a21+a22+gamma*theta));
                YY = (a21+a22)/gamma;
            }
            
            
            if(r00 == "1" && r01 == "2" && r10 == "2" && r11 == "1"){
                XX = (1/(2*gamma*(2*a21+a22+gamma*theta)))*(4*a12*a21+2*a12*a22+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+a11*(2*a21+a22+gamma*theta)+sqrt(4*gamma*theta*(2*a21+a22+gamma*theta)*((a11+2*a12)*(a21+a22)+(a11+a12)*gamma*theta)+(-gamma*theta*(a21+a22+gamma*theta)+a11*(2*a21+a22+gamma*theta)+a12*(4*a21+2*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+2*a12+gamma*theta)))*(4*a12*a21+2*a12*a22-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+a11*(2*a21+a22-gamma*theta)+sqrt(4*gamma*theta*(2*a21+a22+gamma*theta)*((a11+2*a12)*(a21+a22)+(a11+a12)*gamma*theta)+(-gamma*theta*(a21+a22+gamma*theta)+a11*(2*a21+a22+gamma*theta)+a12*(4*a21+2*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "1" && r01 == "2" && r10 == "0" && r11 == "2"){
                XX = (1/(2*gamma*(2*a22+gamma*theta)))*(2*a11*a22+4*a12*a22+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+2*a12+gamma*theta)*(2*(a11+a12)*a22+(a21+a22)*gamma*theta)+(2*(a11+2*a12)*a22+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
                YY = (1/(2*gamma*(a11+2*a12+gamma*theta)))*(2*a11*a22+4*a12*a22-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(a11+2*a12+gamma*theta)*(2*(a11+a12)*a22+(a21+a22)*gamma*theta)+(2*(a11+2*a12)*a22+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
            }
            
            
            if(r00 == "1" && r01 == "2" && r10 == "1" && r11 == "2"){
                XX = (1/(2*gamma*(a21+2*a22+gamma*theta)))*(2*a12*a21+4*a12*a22+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+a11*(a21+2*a22+gamma*theta)+sqrt(4*gamma*theta*(a21+2*a22+gamma*theta)*((a11+2*a12)*(a21+a22)+(a11+a12)*gamma*theta)+(-gamma*theta*(a21+a22+gamma*theta)+a11*(a21+2*a22+gamma*theta)+a12*(2*a21+4*a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+2*a12+gamma*theta)))*(2*a12*a21+4*a12*a22-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+a11*(a21+2*a22-gamma*theta)+sqrt(4*gamma*theta*(a21+2*a22+gamma*theta)*((a11+2*a12)*(a21+a22)+(a11+a12)*gamma*theta)+(-gamma*theta*(a21+a22+gamma*theta)+a11*(a21+2*a22+gamma*theta)+a12*(2*a21+4*a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "1" && r01 == "2" && r10 == "2" && r11 == "2"){
                XX = (1/(2*gamma*(2*(a21+a22)+gamma*theta)))*(2*a11*a21+4*a12*a21+2*a11*a22+4*a12*a22+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*(a21+a22)+gamma*theta)*((a11+2*a12)*(a21+a22)+(a11+a12)*gamma*theta)+(2*a11*(a21+a22)+4*a12*(a21+a22)+a11*gamma*theta+a12*gamma*theta-gamma*theta*(a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(a11+2*a12+gamma*theta)))*(2*a11*a21+4*a12*a21+2*a11*a22+4*a12*a22-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*(a21+a22)+gamma*theta)*((a11+2*a12)*(a21+a22)+(a11+a12)*gamma*theta)+(2*a11*(a21+a22)+4*a12*(a21+a22)+a11*gamma*theta+a12*gamma*theta-gamma*theta*(a21+a22+gamma*theta))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "2" && r10 == "0" && r11 == "0"){
                XX = ((a11+a12-a21-a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a21+a22)*(2*(a11+a12)+gamma*theta)+(a11+a12-a21-a22+gamma*theta)*exp(2))))/(2*gamma*theta);
                YY = ((-a11-a12+a21+a22)*theta-gamma*theta*exp(2)+sqrt(theta*exp(2)*(4*(a21+a22)*(2*(a11+a12)+gamma*theta)+(a11+a12-a21-a22+gamma*theta)*exp(2))))/(4*(a11+a12)+2*gamma*theta);
            }
            
            
            if(r00 == "2" && r01 == "2" && r10 == "1" && r11 == "0"){
                XX = (1/(2*gamma*(a21+gamma*theta)))*(2*a11*a21+2*a12*a21+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*(a11+a12)+gamma*theta)*((a11+a12)*a21+(a21+a22)*gamma*theta)+(2*(a11+a12)*a21+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
                YY = (1/(2*gamma*(2*(a11+a12)+gamma*theta)))*(2*a11*a21+2*a12*a21-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*(a11+a12)+gamma*theta)*((a11+a12)*a21+(a21+a22)*gamma*theta)+(2*(a11+a12)*a21+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "2" && r10 == "2" && r11 == "0"){
                XX = (1/(2*gamma*(2*a21+gamma*theta)))*(4*a11*a21+4*a12*a21+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*(a11+a12)+gamma*theta)*(2*(a11+a12)*a21+(a21+a22)*gamma*theta)+(4*(a11+a12)*a21+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
                YY = (1/(2*gamma*(2*(a11+a12)+gamma*theta)))*(4*a11*a21+4*a12*a21-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*(a11+a12)+gamma*theta)*(2*(a11+a12)*a21+(a21+a22)*gamma*theta)+(4*(a11+a12)*a21+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "2" && r10 == "0" && r11 == "1"){
                XX = (1/(2*gamma*(a22+gamma*theta)))*(2*a11*a22+2*a12*a22+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*(a11+a12)+gamma*theta)*((a11+a12)*a22+(a21+a22)*gamma*theta)+(2*(a11+a12)*a22+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
                YY = (1/(2*gamma*(2*(a11+a12)+gamma*theta)))*(2*a11*a22+2*a12*a22-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*(a11+a12)+gamma*theta)*((a11+a12)*a22+(a21+a22)*gamma*theta)+(2*(a11+a12)*a22+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "2" && r10 == "1" && r11 == "1"){
                XX = ((a11+a12)*(2*(a21+a22)+gamma*theta))/(gamma*(a21+a22+gamma*theta));
                YY = (a21+a22)/gamma;
            }
            
            
            if(r00 == "2" && r01 == "2" && r10 == "2" && r11 == "1"){
                XX = (1/(2*gamma*(2*a21+a22+gamma*theta)))*(4*a11*a21+4*a12*a21+2*a11*a22+2*a12*a22+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*(a11+a12)*gamma*theta*(2*a21+a22+gamma*theta)*(2*(a21+a22)+gamma*theta)+(2*(a11+a12)*(2*a21+a22)+(a11+a12-a21-a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
                YY = (1/(2*gamma*(2*(a11+a12)+gamma*theta)))*(4*a11*a21+4*a12*a21+2*a11*a22+2*a12*a22-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*(a11+a12)*gamma*theta*(2*a21+a22+gamma*theta)*(2*(a21+a22)+gamma*theta)+(2*(a11+a12)*(2*a21+a22)+(a11+a12-a21-a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "2" && r10 == "0" && r11 == "2"){
                XX = (1/(2*gamma*(2*a22+gamma*theta)))*(4*a11*a22+4*a12*a22+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*(a11+a12)+gamma*theta)*(2*(a11+a12)*a22+(a21+a22)*gamma*theta)+(4*(a11+a12)*a22+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
                YY = (1/(2*gamma*(2*(a11+a12)+gamma*theta)))*(4*a11*a22+4*a12*a22-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*gamma*theta*(2*(a11+a12)+gamma*theta)*(2*(a11+a12)*a22+(a21+a22)*gamma*theta)+(4*(a11+a12)*a22+(-a11-a12+a21+a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "2" && r10 == "1" && r11 == "2"){
                XX = (1/(2*gamma*(a21+2*a22+gamma*theta)))*(2*a11*a21+2*a12*a21+4*a11*a22+4*a12*a22+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*(a11+a12)*gamma*theta*(a21+2*a22+gamma*theta)*(2*(a21+a22)+gamma*theta)+(2*(a11+a12)*(a21+2*a22)+(a11+a12-a21-a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
                YY = (1/(2*gamma*(2*(a11+a12)+gamma*theta)))*(2*a11*a21+2*a12*a21+4*a11*a22+4*a12*a22-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*(a11+a12)*gamma*theta*(a21+2*a22+gamma*theta)*(2*(a21+a22)+gamma*theta)+(2*(a11+a12)*(a21+2*a22)+(a11+a12-a21-a22)*gamma*theta-gamma*exp(2)*theta*exp(2))*exp(2)));
            }
            
            
            if(r00 == "2" && r01 == "2" && r10 == "2" && r11 == "2"){
                XX = (1/(2*gamma*(2*(a21+a22)+gamma*theta)))*(4*a11*a21+4*a12*a21+4*a11*a22+4*a12*a22+a11*gamma*theta+a12*gamma*theta-a21*gamma*theta-a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*(a11+a12)*gamma*theta*(2*(a21+a22)+gamma*theta)*exp(2)+(4*a11*(a21+a22)+4*a12*(a21+a22)+a11*gamma*theta+a12*gamma*theta-gamma*theta*(a21+a22+gamma*theta))*exp(2)));
                YY = (1/(2*gamma*(2*(a11+a12)+gamma*theta)))*(4*a11*a21+4*a12*a21+4*a11*a22+4*a12*a22-a11*gamma*theta-a12*gamma*theta+a21*gamma*theta+a22*gamma*theta-gamma*exp(2)*theta*exp(2)+sqrt(4*(a11+a12)*gamma*theta*(2*(a21+a22)+gamma*theta)*exp(2)+(4*a11*(a21+a22)+4*a12*(a21+a22)+a11*gamma*theta+a12*gamma*theta-gamma*theta*(a21+a22+gamma*theta))*exp(2)));
                
            }
            
            break;
    }
    
}

void Populations::getFitness()
{
    for(int p=0; p<numPops; p++){
        for(int i=0; i<numInd; i++){
            pheno_to_fitness(xys[p][i].xx,
                             xys[p][i].yy,
                             opts[p].x_opt,
                             opts[p].y_opt,
                             opts[p].om11,
                             opts[p].om12,
                             xys[p][i].ww);
            //            cout << opts[p].x_opt << endl;
        }
    }
}

void Populations::pheno_to_fitness(double xx, double yy, double xopt, double yopt, double om11, double om12, double &W)
{
    
    
    //    A <- matrix(c(100,50,50,100), 2,2)
    //    A
    //    solve(A)
    double a(om11);
    double b(om12);
    double c(om12);
    double d(om11);
    
    double tmp;
    tmp = 1/(a*d -b*c);
    double OmegaInv[2][2];
    OmegaInv[0][0] = tmp*d;
    OmegaInv[0][1] = tmp*-b;
    OmegaInv[1][0] = tmp*-c;
    OmegaInv[1][1] = tmp*a;
    
    double dx(xx - xopt);
    double dy(yy - yopt);
    
    double r(OmegaInv[0][0]*dx + OmegaInv[1][0]*dy);
    
    double s(OmegaInv[0][1]*dx + OmegaInv[1][1]*dy);
    
    W = exp(-0.5*(r*dx + s*dy));
    
}

void Populations::printPop()
{
    for(int p=0; p<numPops; p++){
        for(int i=0; i<numInd; i++){
            for(int c=0;c<numLoci;c++){
                for(int l=0; l<numAlleles;l++){
                    cout<< "Pop=" << p << "\t" << "Ind=" << i << "\t" << "Locus=" << c << "\t"<< "allele=" << l << "\t(" << pop[p][i][c][l].regulatory<<" , "<<pop[p][i][c][l].coding<<")\t x=" << xys[p][i].xx << "\ty=" << xys[p][i].yy << "\tw=" << xys[p][i].ww << endl;
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
    (*popPtr).numLoci = 2;
    cout<<"# of loci are "<<(*popPtr).numLoci<<endl;
    (*popPtr).numAlleles =2;
    cout<<"# of alleles are "<<(*popPtr).numAlleles<<endl;
    
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
    normal_distribution<double> dist(0, allelic_stdev);
    
    //    cout << "Here's the result: " <<  dist(e2) << endl;
    return dist(e2);
}





















