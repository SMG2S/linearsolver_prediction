#include"generator.hpp"

Eigenvalue::Eigenvalue() : real_part(0), imaginary_part(0)
{}

Eigenvalue::Eigenvalue(double real_part, double imaginary_part) : real_part(real_part), imaginary_part(imaginary_part)
{}

void Eigenvalue::setEigenvalue(double real_part, double imaginary_part){
    this->real_part = real_part;
    this->imaginary_part = imaginary_part;
}

void Eigenvalue::setRealPart(double real_part){
    this->real_part = real_part;
}

void Eigenvalue::setImaginaryPart(double imaginary_part){
    this->imaginary_part = imaginary_part;
}

Eigenvalue* Eigenvalue::getConjugate(){
    return new Eigenvalue(this->real_part,this->imaginary_part);
}

double Eigenvalue::getRealPart(){
    return this->real_part;
}

double Eigenvalue::getImaginaryPart(){
    return this->imaginary_part;
}

void Eigenvalue::printEigenvalue(){
    cout << this->toString() << endl;
}

string Eigenvalue::toString(){
    
    string buffer = to_string(this->real_part);
    if(this->imaginary_part > 0) buffer += '+'+to_string(this->imaginary_part)+'i';
    else if(this->imaginary_part < 0) buffer += to_string(this->imaginary_part)+'i';
    
    return buffer;
}

Spectrum::Spectrum(int size) : size(size)
{
    this->spectrum = (Eigenvalue*)malloc(size*sizeof(Eigenvalue));
}

ClassicSpectrum::ClassicSpectrum(int size) : Spectrum(size)
{
    srand(time(0));
    this->setDensity((rand()%1000)/1000.0);
}

EllipticSpectrum::EllipticSpectrum(int size) : Spectrum(size)
{
    srand(time(0));
    this->setAbRatio((rand()%10000)/1000.0);
    this->setDensity((rand()%1000)/1000.0);
}

ConcentredSpectrum::ConcentredSpectrum(int size) : Spectrum(size)
{
    srand(time(0));
    this->setDensity((rand()%1000)/2000.0);
    this->setRatio((rand()%1000)/200.0);
    this->subDensity = 1 - this->density;
    
    this->nb_points_c2 = int(this->size/(1+this->ratio));
    this->nb_points_c1 = int(this->nb_points_c2 * this->ratio);
    this->nb_points_c2 += this->size - (this->nb_points_c1 + this->nb_points_c2);

}

ClusteredSpectrum::ClusteredSpectrum(int size, int n) : Spectrum(size)
{
    srand(time(0));
    this->n = n;
    this->densities = (double *)malloc(sizeof(double)*n);
    this->cluster_sizes   = (int *)malloc(sizeof(int)*n);
    this->distance_max = (rand()%50000)/100.0;
    this->distance_min = (rand()%(int(distance_max*100.0)))/100;
    for(int i = 0; i<n ; i++){
        this->densities[i] = (rand()%1000)/1000.0;
        this->cluster_sizes[i] = size/n;
    }
    this->cluster_sizes[0] += size%n;
}


void Spectrum::insertEigenvalue(int index, Eigenvalue eigenvalue){
    this->spectrum[index] = eigenvalue;
}

void Spectrum::printSpectrum(){
    for(int i = 0 ; i < this->size ; i++) cout << i << " , " << (this->spectrum[i]).toString() << endl;
}

void ClassicSpectrum::printParameters(){
    cout << "Type : Classic" << endl;
    cout << "Density : " << this->density << endl;
}

void EllipticSpectrum::printParameters(){
    cout << "Type : Elliptic" << endl;
    cout << "Density : " << this->density << endl;
    cout << "A/B ratio : " << this->ab_ratio << endl;
}

void ConcentredSpectrum::printParameters(){
    cout << "Type : Concentred" << endl;
    cout << "Classic density : " << this->density << endl;
    cout << "Concentred density : " << this->subDensity << endl;
}

void ClusteredSpectrum::printParameters(){

    cout << "Type : Clustered" << endl;
    cout << "Number of cluster : " << this->n << endl;
    cout << "Densities :";

    for(int i = 0; i < this->n-1; i++) cout << " " << this->densities[i];
    cout << " " << this->densities[this->n-1] << endl;

    cout << "Size of each cluster :";
    for(int i = 0; i < this->n-1; i++) cout << " " << this->cluster_sizes[i];
    cout << " " << this->cluster_sizes[this->n-1] << endl;
}

void ClassicSpectrum::printType(){
    cout << "Type : Classic of size "<< this->size << endl;
    //this->printSpectrum();   
}

void EllipticSpectrum::printType(){
    cout << "Type : Elliptic of size "<< this->size << endl;
    //this->printSpectrum();
}

void ConcentredSpectrum::printType(){
    cout << "Type : Concentred of size "<< this->size << endl;
    //this->printSpectrum();
}

void ClusteredSpectrum::printType(){
    cout <<"Type : Clustered of size "<< this->size << endl;
    //this->printSpectrum();
}

void Spectrum::printSpectrum(string path_to_file){
    ofstream file;
    file.open(path_to_file);
    file << "%%index real_part   imaginary_part" << endl;
    file << this->size << " " << this->size << " " << this->size << endl;
    for(int i = 0 ; i < this->size ; i++){
        file << i + 1 << " " << this->spectrum[i].getRealPart() << " " << this->spectrum[i].getImaginaryPart() << endl;
    }
}

void Spectrum::setDensity(double density){
    this->density = density;
}

void ClassicSpectrum::freeSpectrum(){
    free(this->spectrum);
}
void EllipticSpectrum::freeSpectrum(){
    free(this->spectrum);
}
void ConcentredSpectrum::freeSpectrum(){
    free(this->spectrum);
}
void ClusteredSpectrum::freeSpectrum(){
    free(this->spectrum);
    free(this->densities);
    free(this->cluster_sizes);
}

void ClassicSpectrum::randomPattern(){
    
    srand(time(0));
    double xy_length = sqrt(this->size/this->density);
    double real,imag;

    int sign;
    
    for(int i = 0; i < this->size; i++){
        sign = rand()%2; if(!sign) sign = -1;
        real = sign*(rand()%((int(xy_length/2))*100)/100.0);
        sign = rand()%2; if(!sign) sign = -1;
        imag = sign*(rand()%((int(xy_length/2))*100)/100.0);
        (this->spectrum[i]).setEigenvalue(real,imag);
    }
}

void EllipticSpectrum::setAbRatio(double ratio){
    this->ab_ratio = ratio;
}

void EllipticSpectrum::randomPattern(){
    
    srand(time(0));
    double p = this->size/this->density;
    double b = sqrt((p*p)/(2*(1+this->ab_ratio*this->ab_ratio)*M_PI*M_PI));
    double a = this->ab_ratio * b;
    double x,y;

    int sign;

    for(int i = 0; i < this->size ; i++){
        
        sign = rand()%2; if(!sign) sign = -1;
        x = sign*(rand()%(int(a)*100)+1)/100.0;
        sign = rand()%2; if(!sign) sign = -1;
        y = sign*b*sqrt(1.0-x*x/(a*a));
        (this->spectrum[i]).setEigenvalue(x,y);
    }    
}

void ConcentredSpectrum::setRatio(double ratio){
    this->ratio = ratio;
}

void ConcentredSpectrum::setDensity(double density)
{
    this->density = density;
    this->subDensity = 1 - density;
}

void ConcentredSpectrum::randomPattern(){
    
    srand(time(0));
    
    double xy_length = sqrt(this->nb_points_c1/this->density);

    double real,imag;
    
    int size_c1 = this->nb_points_c1;
    int size_c2 = this->nb_points_c2;

    int sign;
    
    for(int i = 0; i < size_c1; i++){
        sign = rand()%2; if(!sign) sign = -1;
        real = sign*(rand()%((int(xy_length/2))*100)/100.0);
        sign = rand()%2; if(!sign) sign = -1;
        imag = sign*(rand()%((int(xy_length/2))*100)/100.0);
        (this->spectrum[i]).setEigenvalue(real,imag);
    }
    
    double radius = sqrt(this->nb_points_c2/(this->subDensity * M_PI));
    double x,y;
    double x0 = xy_length/2;
    double y0 = xy_length/2;
    
    if(rand()%20 > 10) x0 *= -1;
    if(rand()%20 > 10) y0 *= -1;
    
    for(int i = size_c1; i < size_c1+size_c2; i++){
        sign = rand()%2; if(!sign) sign = -1;
        x = x0 + sign*(rand()%(int(radius)+1)*100)/100.0;
        sign = rand()%2; if(!sign) sign = -1;
        y = y0 + sign*(rand()%(int(sqrt(radius*radius-(x-x0)*(x-x0))+1)*100))/100.0;
        
        (this->spectrum[i]).setEigenvalue(x,y);
    }
    
}

void ClusteredSpectrum::setDensity(double *densities){
    for(int i = 0; i < n; i++) this->densities[i] = densities[i];
}

void ClusteredSpectrum::setClusterSizes(int *sizes){
    for(int i = 0; i < n; i++) this->cluster_sizes[i] = sizes[i];
}

void ClusteredSpectrum::randomPattern(){
    
    srand(time(0));

    double *x_center = (double *)malloc(n*sizeof(double));
    double *y_center = (double *)malloc(n*sizeof(double));

    double x,y;

    int shift;

    int sign = rand()%2; if(!sign) sign = -1;
    x_center[0] = sign*(rand()%10000)/1000.0;
    sign = rand()%2; if(!sign) sign = -1;
    y_center[0] = sign*(rand()%10000)/1000.0;

    double r = sqrt(cluster_sizes[0]/(densities[0]*M_PI));
    double r2 = r*r;

    for(int j = 0; j < cluster_sizes[0]; j++){
        sign = rand()%2; if(!sign) sign = -1;
        x = sign*rand()%(int(r)*1000)/1000.0 + x_center[0];
        sign = rand()%2; if(!sign) sign = -1;
        y = sign*rand()%((int(sqrt(r2 - (x-x_center[0])*(x-x_center[0])))+1)*1000)/1000.0 + y_center[0];
        (this->spectrum[j]).setEigenvalue(x,y);
        

    }

    shift = cluster_sizes[0];

    for(int i = 1; i < n ; i++){
        sign = rand()%2; if(!sign) sign = -1;
        x_center[i] = x_center[i-1]+sign*(distance_min + rand()%int(distance_max*1000)/1000.0);

        sign = rand()%2; if(!sign) sign = -1;
        y_center[i] = y_center[i-1]+sign*(distance_min + rand()%int(distance_max*1000)/1000.0);

        r = sqrt(cluster_sizes[i-1]/(densities[i-1]*M_PI));
        r2 = r*r;
        for(int j = shift; j < cluster_sizes[i]+shift; j++){
            sign = rand()%2; if(!sign) sign = -1;
            x = sign*rand()%(int(r)*1000)/1000.0 + x_center[i];
            sign = rand()%2; if(!sign) sign = -1;
            y = sign*rand()%((int(sqrt(r2 - (x-x_center[i])*(x-x_center[i])))+1)*1000)/1000.0 + y_center[i];
            (this->spectrum[j]).setEigenvalue(x,y);
        }


        shift += cluster_sizes[i];

    }

    free(x_center);
    free(y_center);
}

