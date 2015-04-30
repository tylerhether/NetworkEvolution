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
    
    // For mutating the regulatory alleles
//    void reg_mu(int indicator_int, char &network_char);
    
    // For Recobmining, mutating, and mating:
    void recombine_mutate_matePop(int *recomb_array, double *mutate_code_array, int *mutate_reg_array, int rolls);
    
    // For migrating:
    void migratePop(double migration_rate);
    
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

void reg_mu(int indicator_int, string &network_char);



/* * * * * * * * * * * * * * * * * * * * * * * *
 *                                              *
 *                                              *
 *                                              *
 *                                              *
 *                 MAIN  FUNCTION               *
 *                                              *
 *                                              *
 *                                              *
 * * * * * * * * * * * * * * * * * * * * * * * */
int main(int argc, char *argv[])
{
    
    // Seed random numbers
    srand((unsigned int)time(NULL));
    
    int numPops(0);
    int numInds(1);
    double mu(0.01);
    double reg_mu(0.001);
    double mu_var(0.01);             // Allow user to define this
    double x1(100);
    double x2(100);
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
    double m_rate(0.0);
//
    cout << "Reading in arguments" << endl <<  "\tNumber of arguments provided = " << (argc-1)/2 << endl;
    for(int i=1; i<argc; i++){
        if(string(argv[i]) == "-npops"){
            numPops = atoi(argv[i+1]);
            cout << "\tNumber of populations = " << numPops << endl;
        }
        if(string(argv[i]) == "-i"){
            numInds = atoi(argv[i+1]);
            cout << "\tNumber of individuals = " << numInds << endl;
        }
        if(string(argv[i]) == "-gen"){
            num_generations = atoi(argv[i+1]);
            cout << "\tNumber of generations = " << num_generations << endl;
        }
        if(string(argv[i]) == "-mu"){
            mu = atof(argv[i+1]);
            cout << "\tMutation rate of coding loci = " << mu << endl;
        }
        if(string(argv[i]) == "-mu_var"){
            mu_var = atof(argv[i+1]);
            cout << "\tAllelic variation introduced by (coding) mutations = " << mu << endl;
        }
        if(string(argv[i]) == "-reg_mu"){
            reg_mu = atof(argv[i+1]);
            cout << "\tMutation rate of regulatory loci = " << reg_mu << endl;
        }
        if(string(argv[i]) == "-m_rate"){
            m_rate = atof(argv[i+1]);
            cout << "\tMigration rate = " << m_rate << endl;
        }
        if(string(argv[i]) == "-x_start"){
            x1 = atoi(argv[i+1]);
            cout << "\tInitial mean trait x value = " << x1 << endl;
        }
        if(string(argv[i]) == "-y_start"){
            x2 = atoi(argv[i+1]);
            cout << "\tInitial mean trait y value = " << x2 << endl;
        }
        if(string(argv[i]) == "-start_network"){
            reg_pattern = argv[i+1];
            cout << "\tInitial network architecture = " << reg_pattern << endl;
        }
        if(string(argv[i]) == "-theta"){
            theta = atoi(argv[i+1]);
            cout << "\ttheta = " << theta << endl;
        }
        if(string(argv[i]) == "-gamma"){
            cout << "\tgamma = " << argv[i+1] << endl;
            theta = atoi(argv[i+1]);
        }
        if(string(argv[i]) == "-model"){
            mod = argv[i+1];
            cout << "\tRegulatory model = " << mod << endl;
        }
        if(string(argv[i]) == "-start_deviation"){
            allelic_Stdev = atoi(argv[i+1]);
            cout << "\tInitial standard deviation of allelic values = " << allelic_Stdev << endl;
        }
        if(string(argv[i]) == "-x_opt"){
            XOPT = atoi(argv[i+1]);
            cout << "\tTriat x's optimum = " << XOPT << endl;
        }
        if(string(argv[i]) == "-y_opt"){
            YOPT = atoi(argv[i+1]);
            cout << "\tTrait y's optimum = " << YOPT << endl;
        }
        if(string(argv[i]) == "-sel_var"){
            OM11 = atoi(argv[i+1]);
            cout << "\tVariance in stabilizing selection = " << OM11 << endl;
        }
        if(string(argv[i]) == "-sel_covar"){
            OM12 = atoi(argv[i+1]);
            cout << "\tCovariance in stabilizing selection = " << OM12 << endl;
        }
        if(string(argv[i]) == "-rec"){
            rec = atof(argv[i+1]);
            cout << "\tRecombination rate = " << rec << endl;
        }

    }
    
    // Generate a random sequence of 1s and 0s for recombination:
    cout << "Finished reading in arguments." << endl << endl << "Building recombination and mutation arrays";
    
    int nrolls=50000+(4*numInds);                        // There are 50000 possible starting points
    cout << " with " << nrolls << " elements..." << endl;
    int *mRecombinationArray = new int[nrolls];                 // This will store binary array (1=recombination, 0=no recomb.)
    int *mRecombinationArray2 = new int[nrolls];                // This will store the modulus 2 output

    double *mMutateCodingArray = new double[nrolls];                  // This array will store double of mutational effects
    int *mMutateRegulatoryArray = new int[nrolls];                 // For all 1s above, this will store the change in allelic value
    
    random_device generator;
    mt19937 e2(generator());

    // This section fills in the recombination array with whether (1) or
    // not (0) a recombinatino event occurs given the recombination rate, rec.
    bernoulli_distribution distribution(rec);
    for (int i=0; i<nrolls; i++){
        
        if (distribution(generator)){
            
            mRecombinationArray[i]=1;
           //  cout << mRecombinationArray[i] << " " << endl;
        } else {
            mRecombinationArray[i]=0;
            // cout << mRecombinationArray[i] << endl;
        }
    }
    
    // This performs the modulus 2 so that successive integers
    // indicate which of the two allles to use. Thanks AM for the tip.
    for(int i=0; i<nrolls; i++){
        int sum(0);
        for(int j=0; j<i; j++){
            sum+=mRecombinationArray[j];
        }
        mRecombinationArray2[i] = sum % 2;
    }
    
    /* Debugging:
     * for(int i=0; i<nrolls; i++){
     *     cout << mRecombinationArray[i] << "\t" << mRecombinationArray2[i] << endl;
     * }
     */

    /* For coding and regulatory mutations we do a similar thing
     * as the recombination rates except: 
     * 1) failure is coded as 1.0 and
     * 2) success (mutation) is given my a normal distribution with mean
     * of one and variance of mu_var
     * From 1 and 2 this array becomes one of scaling factors where values!=1 reflect the change
     * of allelic values from pre-mutated coding allelic values */
    
    random_device generator2;                  // this defines the dist for whether mutation happens
    mt19937 e22(generator2());
    bernoulli_distribution distributionMU(mu);
    
    random_device rd;                          // this defines the dist for the mutational effects
    mt19937 e222(rd());                        // if a mutation happened
    normal_distribution<double> dist2(1, mu_var);
    
  
    for (int i=0; i<nrolls; i++){
       // cout << "This is the value: " << dist2(e22) << endl;
        if (distributionMU(generator2)){
            
            mMutateCodingArray[i]= dist2(e222);
             //  cout << mMutateCodingArray[i] << " " << endl;
        } else {
            mMutateCodingArray[i]= 1;
             //  cout << mMutateCodingArray[i] << endl;
        }
    }
    
    // For regulatory mutations we also do a similar thing.
     
    random_device generator3;                       // This block is for whether a mutation happens at all
    mt19937 e2222(generator3());
    bernoulli_distribution m_rateDist(reg_mu);
    
    random_device zero_one_two;                     // This block is for the actual new value
    mt19937 e3(zero_one_two());
    uniform_int_distribution<> Dist0_2(0,2);
    
    
    //int tmpcounter(0);
    for (int i=0; i<nrolls; i++){
        if(m_rateDist(generator3)){
            mMutateRegulatoryArray[i] = Dist0_2(zero_one_two);
            //cout << mMutateRegulatoryArray[i] << endl;
            //tmpcounter++;
        } else {
            mMutateRegulatoryArray[i] = 4;
            //cout << mMutateRegulatoryArray[i] << endl;
        }

    }
    //cout << "here's the freq estimate of reg_mu " << static_cast<double>(tmpcounter)/static_cast<double>(nrolls) << endl;
    
    
    

    /* Running the simulations: **********************************************************
     * This initializes the populations with their starting genotypes and phenotypes and *
     * calculates their initial fitness values. It then runs through each generation and *
     * caluclates phenotypes from genotypes; calculates fitness from phenotypes; selects *
     * based on fitness; mates, recombines, and mutates; migrates; and outputs summary   *
     * statistics. Between the calls to fitness and selection hybrid incompatibility is  *
     * estimated as well.*****************************************************************/
    
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
        
        // Print Progress
        if(g % 100 == 0) cout << "Generation " << g << endl;
        
        
        // Get New Phenotypes
        Pop.getPheno(theta, gamma, *mod);
        
        // Calculate Fitness:
        Pop.getFitness();

        // Viability Selection:
        Pop.selection();
        
        // Recombine, mutate, and mate the survivers
        Pop.recombine_mutate_matePop(mRecombinationArray2, mMutateCodingArray, mMutateRegulatoryArray, nrolls);

        // Need migration
        Pop.migratePop(m_rate);
        // Need to estimate hybrid fitness
        
        // Need output summary stats to file
        
    }
//    char test = '2';
//    cout << "test is  " << test << endl;
//    Pop.reg_mu(1, test);
//    cout << "test is now " << test << endl;
    
    
    Pop.getPheno(theta, gamma, *mod);
    Pop.getFitness();
    Pop.printPop();
    // Cleaning up:
//    Pop.deletePops_XYs();
   
    /* Troubleshooting: run this block of code to determine the version of C++ that is used:
     * if( __cplusplus == 201103L ) std::cout << "C++11\n" ;
     * else if( __cplusplus == 199711L ) std::cout << "C++98\n" ;
     * else std::cout << "pre-standard C++\n" ;*/
    delete[] mRecombinationArray;
    delete[] mRecombinationArray2;
    delete[] mMutateCodingArray;
    delete[] mMutateRegulatoryArray;
    return 0;
}



/* * * * * * * * * * * * * * * * * * * * * * * *
 *                                              *
 *                                              *
 *                                              *
 *                                              *
 *                 CLASS FUNCTIONS              *
 *                                              *
 *                                              *
 *                                              *
 * * * * * * * * * * * * * * * * * * * * * * * */
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
            double TEST = rand()/static_cast<double>(RAND_MAX);

            if(TEST<xys[p][i].ww){ // Individual lives
         
                for(int c=0; c<numLoci; c++){
                    for(int a=0; a<numAlleles; a++){
                        // Put the surviving individual in the after selection bin:
                        pop_after_selection[p][numIndAS[p]][c][a].coding = pop[p][i][c][a].coding;
                        pop_after_selection[p][numIndAS[p]][c][a].regulatory = pop[p][i][c][a].regulatory;
                    }
                }
                numIndAS[p]++; // Keeps track of the number of individuals that survived selection
            } /* else
               * {
               * cout << "This individual didn't make it: " << TEST << " is less than " << xys[p][i].ww << endl;
               * }
               */
            // Enhancement: check to make sure that the average fitness of the non-surviving individuals < w of surving inds.
        }
        if(numIndAS[p]==0){ // If one population went extinct then a hole is torn in the universe.
            cout << "Wait, selection killed off all of population #"<< p << ". Aborting simulation now.\nConsider weaker selection or larger number of individuals." << endl;
            exit(1);
        }
        // Debugging:
//        cout << "The number of individuals retained in population " << p << " = " << numIndAS[p] << endl;
    }
}

void Populations::migratePop(double migration_rate)
{
    
    //cout << " before swap one is " << pop[0][0][0][0].regulatory << " and two is " << pop[0][0][1][0].regulatory << endl;
    //swap(pop[0][0][0][0].regulatory,pop[0][0][1][0].regulatory);
    //cout << " after swap one is " << pop[0][0][0][0].regulatory << " and two is " << pop[0][0][1][0].regulatory << endl;

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


void Populations::recombine_mutate_matePop(int *recomb_array, double *mutate_code_array, int *mutate_reg_array, int rolls)
{
    /* Here is the random number generator that samples the starting position
     * of the recombination array (2) */
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, static_cast<int>(rolls-(4*numInd+1)));
    
    // This generator is for picking the parents
    random_device rd2;
    mt19937 gen2(rd2());

    /* So for each population we need to mate the surviving individuals at random.
     * The number of surviving individuals for each p populations is in the integer
     * array 'numIndAS' */
    for(int p=0; p<numPops; p++)
    {
        // Tailor the unif dist based on # survivers in population p
        uniform_int_distribution<> pick_pairs(0, (numIndAS[p]-1)); //cout << "The total number of individuals in pop " << p << " are " << numIndAS[p] << endl;
        
        // Pick the starting index of the recombination rates array (2)
        int index(dis(gen)); // cout << "The starting index is " << index << endl;
        
        /* From the surviving population we need to 'refill' a new one with 'numInd'
         * individuals. These new individuals have half of their genomes derived
         * from each of their parents' gametes */
        for (int i=0; i<numInd; i++)
        {

            uniform_int_distribution<> pick_pairs(0, static_cast<int>(numIndAS[p]-1));
            
            int parent1(pick_pairs(gen2));
            int parent2(pick_pairs(gen2));
            
            //cout << "The first parent is # " << parent1 << " and its picked allele at the first locus is " << recomb_array[index] << " (" <<pop_after_selection[p][parent1][0][recomb_array[index]].coding << ")" << endl;
            pop[p][i][0][0].coding =     (pop_after_selection[p][parent1][0][recomb_array[index]].coding)*(mutate_code_array[index]);
            pop[p][i][0][0].regulatory = pop_after_selection[p][parent1][0][recomb_array[index]].regulatory;
            reg_mu(mutate_reg_array[index], pop[p][i][0][0].regulatory); // This function mutates the regulatory region

            
            //cout << "The first parent is # " << parent1 << " and its picked allele at the second locus is " << recomb_array[index+1] << " (" << pop_after_selection[p][parent1][1][recomb_array[index+1]].coding << ")" << endl;
            pop[p][i][1][0].coding =     pop_after_selection[p][parent1][1][recomb_array[index+1]].coding*(mutate_code_array[index+1]);
            pop[p][i][1][0].regulatory = pop_after_selection[p][parent1][1][recomb_array[index+1]].regulatory;
            reg_mu(mutate_reg_array[index+1],pop[p][i][1][0].regulatory); // This function mutates the regulatory region
            
            //cout << "The second parent is # " << parent2 << " and its picked allele at the first locus is " << recomb_array[index+2] << " (" << pop_after_selection[p][parent2][0][recomb_array[index+2]].coding<< ")" << endl;
            pop[p][i][0][1].coding =     pop_after_selection[p][parent2][0][recomb_array[index+2]].coding*(mutate_code_array[index+2]);
            pop[p][i][0][1].regulatory = pop_after_selection[p][parent2][0][recomb_array[index+2]].regulatory;
            reg_mu(mutate_reg_array[index+2],pop[p][i][0][1].regulatory); // This function mutates the regulatory region
            
            //cout << "The second parent is # " << parent2 << " and its picked allele at the second locus is " << recomb_array[index+3] << " (" << pop_after_selection[p][parent2][1][recomb_array[index+2]].coding << ")" << endl;
            pop[p][i][1][1].coding =     pop_after_selection[p][parent2][1][recomb_array[index+2]].coding*(mutate_code_array[index+3]);
            pop[p][i][1][1].regulatory = pop_after_selection[p][parent2][1][recomb_array[index+2]].regulatory;
            reg_mu(mutate_reg_array[index+3],pop[p][i][1][1].regulatory); // This function mutates the regulatory region
            
            index+=4;
        }
    }
}

/* * * * * * * * * * * * * * * * * * * * * * * *
 *                                              *
 *                                              *
 *                                              *
 *                                              *
 *               GENERAL FUNCTIONS              *
 *                                              *
 *                                              *
 *                                              *
 * * * * * * * * * * * * * * * * * * * * * * * */
void input(Populations *popPtr, int POPS, int INDS)
{
    
    (*popPtr).numPops = POPS;
    // cout<<POPS<<endl;
    //cout<<"# of populations are "<<(*popPtr).numPops<<endl;
    (*popPtr).numInd = INDS;
    //cout<<"# of individuals are "<<(*popPtr).numInd<<endl;
    (*popPtr).numLoci = 2;
    //cout<<"# of loci are "<<(*popPtr).numLoci<<endl;
    (*popPtr).numAlleles =2;
    //cout<<"# of alleles are "<<(*popPtr).numAlleles<<endl;
    
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
    
    
//    cout << "reg = " << reg_pattern << " x = " << x << " y = " << y << " theta = " << theta << " gamma= " << gamma << " model= " << mod << " a1= " << a1 << " a2= " << a2 << endl;
//    x = 20.3;
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

// This function actually does the mutating of the regulatory networks
void reg_mu(int indicator_int, string &network_char)
{
    if(indicator_int==0)
    {
        network_char = "0";
    } else if(indicator_int==1)
    {
        network_char = "1";
    } else if(indicator_int==2)
    {
        network_char = "2";
    } else {
        network_char = network_char;
    }
}













