#ifndef GENEREATOR_H
#define GENERATOR_H

#include<iostream>
#include<fstream>
#include<time.h>
#include<math.h>
#include<string>

using namespace std;

class Eigenvalue
{
    public:
        double real_part;
        double imaginary_part;
    
        //constructor    
        Eigenvalue();
        Eigenvalue(double real_part, double imaginary_part);
        
        //methods
        Eigenvalue* getConjugate();
        
        void setEigenvalue(double real_part, double imaginary_part);
        void setRealPart(double real_part);
        void setImaginaryPart(double imaginary_part);
        
        double getRealPart();
        double getImaginaryPart();
        
        void printEigenvalue();
        
        string toString();
};

class Spectrum
{
    public:
        int size;
        double density;
        
        Eigenvalue *spectrum;
        
        //constructors
        Spectrum(int size);
        
        virtual ~Spectrum(){};
        
        //methods
        void printSpectrum();
        void printSpectrum(string path_to_file);
        
        void insertEigenvalue(int index, Eigenvalue eigenvalue);
        
        void setDensity(double density);
        
        virtual void randomPattern() = 0;
        virtual void printParameters() = 0;
        virtual void printType() = 0;
        virtual void freeSpectrum() = 0;

};

class ClassicSpectrum : public Spectrum
{
    public :
        
        ClassicSpectrum(int size);
        void randomPattern();
        void printType();
        void printParameters();
        void freeSpectrum();
};

class EllipticSpectrum : public Spectrum
{
    public:
        double ab_ratio;
        EllipticSpectrum(int size);
        void setAbRatio(double ratio);
        void randomPattern();
        void printType();
        void printParameters();
        void freeSpectrum();
};

class ConcentredSpectrum : public Spectrum
{
    public:
        double ratio = 2;
        double subDensity;
        
        int nb_points_c1;
        int nb_points_c2;
        
        ConcentredSpectrum(int size);
        void setRatio(double ratio);
        void setDensity(double ratio);
        void randomPattern();
        void printType();
        void printParameters();
        void freeSpectrum();
};

class ClusteredSpectrum : public Spectrum
{
    public :
        int n;
        int *cluster_sizes;
        double *densities;
        
        double distance_max = 100.0;
        double distance_min = 60.0;

        ClusteredSpectrum(int size, int n);
        void setClusterSizes(int *cluster_sizes);
        void setDensity(double *densities);
        void randomPattern();

        void printType();
        void printParameters();
        void freeSpectrum();
};

#endif
