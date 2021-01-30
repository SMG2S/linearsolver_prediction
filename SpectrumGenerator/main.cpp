#include<iostream>
#include"generator.hpp"

int main(int argc, char **argv){

    int n = 23;

    std::string buffer = "Clustered";

    Spectrum *Sigma;

    if(buffer == "Concentred"){
        Sigma = new ConcentredSpectrum(n);
    }
    else if(buffer == "Classic"){
        Sigma = new ClassicSpectrum(n);
    }
    else if(buffer == "Elliptic"){
        Sigma = new EllipticSpectrum(n);
    }
    else if(buffer == "Clustered"){
        Sigma = new ClusteredSpectrum(n,5);
    }

    Sigma->printParameters();
    
    buffer = "Spectrum/"+buffer+"/spectrum_"+to_string(n)+".txt";
    //Sigma->setDensity(0.001);

    Sigma->randomPattern();

    Sigma->printSpectrum();

    

    Sigma->printSpectrum(buffer);
    Sigma->freeSpectrum();

    delete(Sigma);
    
    return 0;

}
