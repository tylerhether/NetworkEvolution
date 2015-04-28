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
    int numPops,numInd,numLoci,numAlleles; // AS = after selection
    int *numIndAS;//[10];

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
    
    // For selection:
    locus**** pop_after_selection;
//    int nIndAS[];
    void selection();
    
    
    // For Recobmining, mutating, and mating:
    void recombine_mutate_matePop(int recomb_array[], int rolls);
    
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
    
    // Seed random numbers
    srand((unsigned int)time(NULL));
    
    int numPops(0);
    int numInds(1);
    double x1(1);
    double x2(1);
    char* reg_pattern = nullptr;
    double theta(300);
    double gamma(1);
    char *mod = nullptr;
    int num_generations(1);
    double allelic_Stdev(1);
    double XOPT(300);                // Trait 1 optimum
    double YOPT(300);                // Trait 2 optimum
    double OM11(10000);                // Variance in stabilizing selection (for all pops)
    double OM12(0);
    double rec(0.5);                  // This is the recombination rate between coding loci
    
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
            OM12 = atoi(argv[i+1]);
        }
        if(string(argv[i]) == "-r"){
            cout << "Recombination rate = " << argv[i+1] << endl;
            rec = static_cast<double>(atoi(argv[i+1]))/100;
        }
        
    }
    
    // Generate a random sequence of 1s and 0s for recombination:
    
    cout << "Building recombination sequences";
    
    int nrolls=numPops*numInds*num_generations < 50000 ? numPops*numInds*num_generations : 50000; // Is this enough?
    cout << " with " << nrolls << " elements." << endl;
    int *mRecombinationArray = new int[nrolls];
    int *mRecombinationArray2 = new int[nrolls];
//    int mRecombinationArray2[nrolls] = { 0 };
    random_device generator;
    mt19937 e2(generator());

    bernoulli_distribution distribution(rec);
    
    for (int i=0; i<nrolls; i++){
        
        if (distribution(generator)){
            
            mRecombinationArray[i]=1;
//             cout << mRecombinationArray[i] << " " << endl;
        } else {
            mRecombinationArray[i]=0;
//             cout << mRecombinationArray[i] << endl;
        }
    }
    
    for(int i=0; i<nrolls; i++){
        int sum(0);
        for(int j=0; j<i; j++){
            sum+=mRecombinationArray[j];
        }
        mRecombinationArray2[i] = sum % 2;
    }
    
    for(int i=0; i<nrolls; i++){
        cout << mRecombinationArray[i] << "\t" << mRecombinationArray2[i] << endl;
    }
    

    

    
    
    
    
    
    
    
    
    
    
    
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
    
    Pop.getFitness();
    Pop.printPop();
    // Recursion:
    for(int g=0; g<num_generations; g++){
        
        if(g % 100 == 0){
            cout << "Generation " << g << endl;
        }
        // Get new phenotypes
        Pop.getPheno(theta, gamma, *mod);
        // Fitness:
        Pop.getFitness();
        //            num_generations = num_generations + 1 -1;
                             // Troubleshooting: print populations
        // Need selection
       
        // Need mutation
        Pop.recombine_mutate_matePop(mRecombinationArray2, nrolls);
//        Pop.getPheno(theta, gamma, *mod);
//        Pop.getFitness();
        // Need mating

        // Need migration
        
        // Need to estimate hybrid fitness
        
        // Need output summary stats to file
        Pop.selection();
        
    }
    

    
//     Pop.printPop();
    // Cleaning up:
//    Pop.deletePops_XYs();
   
    /* Troubleshooting: run this block of code to determine the version of C++ that is used:
     * if( __cplusplus == 201103L ) std::cout << "C++11\n" ;
     * else if( __cplusplus == 199711L ) std::cout << "C++98\n" ;
     * else std::cout << "pre-standard C++\n" ;*/
    delete[] mRecombinationArray;
    delete[] mRecombinationArray2;
    return 0;
}




//Class Functions
Populations::Populations()
{
    pop=NULL;
}

void Populations::initilizePop(string reg_pattern, double theta, double gamma, char mod, double a1, double a2, double allelic_Stdev)
{
    numIndAS=new int[numPops];
    
    pop=new locus***[numPops];
    pop_after_selection=new locus***[numPops];
    for(int p=0; p<numPops; p++){
        numIndAS[p]=0;
        pop[p]=new locus**[numInd];
        pop_after_selection[p]=new locus**[numInd];
        for(int i=0; i<numInd; i++){
            pop[p][i]=new locus*[numLoci];
            pop_after_selection[p][i]=new locus*[numLoci];
            int counter(0);
            for(int c=0;c<numLoci;c++){
                pop[p][i][c]=new locus[numAlleles];
                pop_after_selection[p][i][c]=new locus[numAlleles];
                
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
                    pop_after_selection[p][i][c][l].coding = 0;
                    pop_after_selection[p][i][c][l].regulatory = "1";
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
                delete[] pop_after_selection[p][i][c];
            }
            delete pop[p][i];
            delete pop_after_selection[p][i];
        }
        delete pop[p];
        delete xys[p];
        delete pop_after_selection[p];
//        delete numIndAS;
    }
    delete pop;
    delete pop_after_selection;
    delete xys;
    delete opts;
    delete[] numIndAS;
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
    
    double a = om11;
    double b = om12;
    double c = om12;
    double d = om11;
    double tmp;
    tmp = 1.0/(a*d -b*c);
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

void Populations::selection(){
    for(int p=0; p<numPops; p++){
        numIndAS[p]=0;
        for(int i=0; i<numInd; i++){
            
            if(rand()/static_cast<double>(RAND_MAX)<xys[p][i].ww){ // Individual lives
         
                for(int c=0; c<numLoci; c++){
                    for(int a=0; a<numAlleles; a++){
                        // Put the surviving individual in the after selection bin:
                        pop_after_selection[p][numIndAS[p]][c][a].coding = pop[p][i][c][a].coding;
                        pop_after_selection[p][numIndAS[p]][c][a].regulatory = pop[p][i][c][a].regulatory;
                    }
                }
                numIndAS[p]++; // Keeps track of the number of individuals that survived selection
            }
            if(numIndAS[p]==0){ // If one population went extinct then a hole is torn in the universe.
                cout << "Wait, selection killed off all of population #"<< p << ". Aborting simulation now.\nConsider weaker selection or larger number of individuals." << endl;
                exit(1);
            }
            // Enhancement: check to make sure that the average fitness of the non-surviving individuals < w of surving inds.
        }
        // Debugging:
//        cout << "The number of individuals retained in population " << p << " = " << numIndAS[p] << endl;
    }

    
    
    
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








void Populations::recombine_mutate_matePop(int recomb_array[], int rolls)
{

}






//
//
//void Populations::recombine_mutate_matePop(int recomb_array[], int rolls)
//{
//    
//
//    /* Here is the random number generator that samples the starting position
//     * of the recombination array (2) */
//    random_device rd;
//    mt19937 gen(rd());
//    uniform_int_distribution<> dis(0, static_cast<int>(rolls-(4*numInd+1)));
////    for (int n=0; n<1000; ++n)            // For troubleshooting:
////        std::cout << dis(gen) << ' ';     // This block prints out the random integers
////    std::cout << '\n';
////    
//    
////    uniform_int_distribution<> pick_pairs(0, 20);
////    for (int n=0; n<1000; ++n)            // For troubleshooting:
////        std::cout << pick_pairs(gen) << ' ';     // This block prints out the random integers
////    std::cout << '\n';
////    
//
//
//    /* So for each population we need to mate the surviving individuals at random.
//     * The number of surviving individuals for each p populations is in the integer
//     * array 'numIndAS' */
//    for(int p=0; p<numPops; p++)
//    {
//        // Tailor the unif dist based on # survivers in population p
//        uniform_int_distribution<> pick_pairs(0, (numIndAS[p]-1)); cout << "The total number of individuals in pop " << p << " are " << numIndAS[p] << endl;
//        
//        
//        // Pick the starting index of the recombination rates array (2)
//        int index(dis(gen)); // cout << "The starting index is " << index << endl;
//        
//        /* From the surviving population we need to 'refill' a new one with 'numInd'
//         * individuals. These new individuals have half of their genomes derived
//         * from each of their parents' gametes */
//     
//        for (int i=0; i<numInd; i++)
//        {
//            random_device rd2;
//            mt19937 gen(rd2());
//            uniform_int_distribution<> pick_pairs(0, static_cast<int>(numIndAS[p]-1));
//            
//            int parent1(pick_pairs(gen));
//            int parent2(pick_pairs(gen));
//            
////            parent1 = pick_pairs(e2);
////            const int parent1 = static_cast<int>(pick_pairs(e2));
////            const int parent2 = pick_pairs(e2);
////
//            
////            cout << "The picked allele is # " << recomb_array[index] << endl;
//            
//            cout << "The first parent is # " << parent1 << " and its picked allele at locus one is " << recomb_array[index] << endl;
//            cout << "The first parent's index plus 1 is: " << parent1 << " + 1 = " << parent1+1 << endl;
////            pop[p][i][0][0].regulatory =
//            pop_after_selection[p][parent1][0][0].coding =1.0;// << endl;
//            
//            cout << "The first parent is # " << parent1 << " and its picked allele at locus two is " << recomb_array[index+1] << endl;
//            
//            index++; index++; // This moves the index over two spaces
//            
//            cout << "The second parent is # " << parent2 << " and its picked allele at locus one is " << recomb_array[index] << endl;
//            cout << "The second parent is # " << parent2 << " and its picked allele at locus two is " << recomb_array[index+1] << endl;
//            
//            index++; index++; // This moves the index over two spaces
//            
//            
//            
//            
////            cout << "Picked parent one = " << parent1 << ", parent two " << parent2 << endl;
//            
//            cout << pop[p][i][0][0].regulatory << endl;
////            cout << pop_after_selection[p][1][0][0].regulatory << endl;
//
//            
////            cout << pop_after_selection[p][parent1][0][recomb_array[index]].regulatory << endl;
//            
////            recomb_array[index]
//            
//            
//        }
////        delete[] surving_parents;
//    }
//    
//    
//    
//    
    // This section below needs work.
    
//    for(int p=0; p<numPops; p++){
//        // Where along the recombination array (2) should we start the indexing?
//        int index(dis(gen));
////        cout << "this is the starting index " << index << endl;
//        cout << " Population " << p << " has " << numIndAS[p] << " inds." << endl;
//        for(int ii=index; ii<(index+numInd); ii++){
//            cout << recomb_array[ii] << endl;
//            ii++; // This is on purpose
//        }
//        
//        for(int i=0; i<numInd; i++){
//            // Replace the value of pop[p][i] with a (mutated) recombinant from the after selection survivers
//            
//            
//            
//            for(int c=0;c<numLoci;c++){
//                for(int l=0; l<numAlleles;l++){
//                    pop[p][i][c][l].coding=pop[p][i][c][l].coding+0.00001;//a1/numAlleles+make_genos(a1/numAlleles, allelic_Stdev)
//                }
//                //                cout<<endl;
//            }
////            cout<<endl;
//        }
//    }
//}

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



















