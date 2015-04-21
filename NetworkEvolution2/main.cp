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
#include <cstdlib> // allows us to halt with exit() function


using namespace std; // This let's you use the shorthand the std library's functions

//Clasess and Strcutures
struct locus                                // Initialize the locus structure
{
    int regulatory;                         // regulatory alleles can only be the integers 0,1,2
    double coding;                          // coding alleles are continuous
};

class Populations                            // Initialize the population class
{
//private:                                  // Shouldn't we have some private members?
public:

    locus**** pop;                           // Define the structure of the population. Why 3 ***s?
    int numPops,numInd,numChromo,numLoci;   // Initialize variables that will make numPops populations
    Populations();                           // Start with a NULL population
    void initilizePop();                    // Initialize population
    void deletePop();                       // Destroy population. Can't we just used a destructor here? 
    void printPop();                        // Print population to screen (for debugging)
    
};

// Forward Declarations of Functions:       // Can be moved to a header file if gets too long.
void input(Populations *popPtr);             // These are function prototypes
void printLocus(locus Locus);               // These are function prototypes

// Main Function to run:
int main(int argc, char *argv[])
{
    srand((unsigned int)time(NULL));        // Seeding Random
    Populations ourPop1;                     // Initialize Populations
    
    // Running:
    input(&ourPop1);
    ourPop1.initilizePop();
    ourPop1.printPop();
    ourPop1.deletePop();
    
    
    
    
    
    
    return 0;
}

//Class Functions
Populations::Populations()
{
    pop=NULL;
}
void Populations::initilizePop()
{
    pop=new locus***[numInd];
    for(int p=0; p<numPops; p++){
        pop[p]=new locus**[numInd];
        for(int i=0; i<numInd; i++){
            pop[p][i]=new locus*[numChromo];
            for(int c=0;c<numChromo;c++){
                pop[p][i][c]=new locus[numLoci];
                for(int l=0; l<numLoci;l++){
                    pop[p][i][c][l].coding=0.01;           // Seed coding allele expression
                    pop[p][i][c][l].regulatory=rand() % 3; // Seed regulatory allele
                    //                pop[i][c][l].sex=rand() % 2;        // Seed sex
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
void input(Populations* popPtr)
{
    fstream inFile;
    inFile.open("/Users/tylerhether/Projects/09-NetworkEvolutioncpp/NetworkEvolution/NetworkEvolution2/Input.txt");
//    inFile.open(argv[1]);
    string myString;
    myString="Nothing";
    
    while(myString!= "populations:" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*popPtr).numPops;
    cout<<"# of populations are "<<(*popPtr).numPops<<endl;
    
    
    while(myString!= "individuals:" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*popPtr).numInd;
    cout<<"# of individuals are "<<(*popPtr).numInd<<endl;
    
    
    while(myString!= "chromosomes:" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*popPtr).numChromo;
    cout<<"# of chromosomes are "<<(*popPtr).numChromo<<endl;
    
    
    while(myString!= "loci:" && inFile.good())
    {
        inFile>>myString;
    }
    inFile>>(*popPtr).numLoci;
    // (*popPtr).numLoci = 0;
    cout<<"# of loci are "<<(*popPtr).numLoci<<endl;
    
}
void printLocus(locus Locus)
{
    cout<<"Locus"<<endl;
    cout<<"regulatory =" <<Locus.regulatory<<endl;
    cout<<"coding = " <<Locus.coding<<endl;
    
}





























