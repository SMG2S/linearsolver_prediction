#include<iostream>
#include<time.h>
#include"solver.hpp"

static char help[] = "Solve a tridiagonal linear system with KSK.\n\n";

void InsertNewLine(ofstream *file, MatrixFeatures *Matrix, VectorFeatures *Vector){
	
	std::string buffer = "";
	char sep = ',';
	Matrix->PrintInputFeatures(&buffer);
    AddSep(&buffer,sep);
    Vector->PrintFeatures(&buffer);
    AddSep(&buffer,sep);
    Matrix->PrintOutputFeatures(&buffer);
    *file << buffer << endl;
}

void InsertFirstLine(ofstream *file){
    *file << "n" << "," << "nnz" << "," << "sparsity" << ",";
    *file << "trace_real" << "," << "trace_imag" << "," << "amplitude_mean" << "," << "amplitude_variance" << ",";
    *file << "mat_min.mod" << "," << "mat_max.mod" << "," << "diag_min.mod" << "," << "diag_max.mod" << ",";
    *file << "diag_mean" << "," << "diag_variance" << "," << "diag_gravity" << ",";
    *file << "lbandwidth" << "," << "rbandwidth" << "," << "solver" << ",";
    *file << "vec_mean" << "," << "vec_variance" << "," << "gravity_center" << "," << "vec_min_mod" << "," << "vec_max_mod";
    *file << "converged,iterations,elapsed_time" << endl;

}

int main(int argc, char **argv){

	srand(time(NULL));

	PetscInitialize(&argc,&argv,(char*)0,help);

    PetscViewer viewer;
	PetscViewerCreate(PETSC_COMM_WORLD,&viewer);


	//Initializing Parameters
    
    Vec x,u;
    Mat A;
    KSP kspgmres,kspbicg;

    std::string data_file_path = "sample.txt";
    ofstream file(data_file_path, ios::out | ios::app);

    int res;

    PetscReal norm,tolerance = 1.0e-6,mean_norm2_tol = 1.0e-1,density = (rand()%1000)/1000.0;
    
    int size, rank; MPI_Comm_size(PETSC_COMM_WORLD,&size); MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	MatrixFeatures *Matrix = new MatrixFeatures; Matrix->n = rand()%5001 + 20;
	VectorFeatures *Vector = new VectorFeatures; Vector->n = Matrix->n;

	PetscInt l = rand()%int(Matrix->n/1000)+20,d = l/2,its;
	
	std::string spectrum_type = "Classic";
	
	std::string spectrum_file_name = "SpectrumGenerator/Spectrum/" + spectrum_type + "/spectrum_" + to_string(Matrix->n) + ".txt";
        
    if(!rank) std::cout << "=========" << spectrum_type << "=========" << endl;

    //if(!rank) InsertFirstLine(&file);

    //Initializing Matrix and RHS

    if(spectrum_type == "Classic") InitClassicSpectrum(Matrix->n,density,spectrum_file_name);

    else if(spectrum_type == "Elliptic"){
    	float ab_ratio = ((rand()%50)+1)/10.0;	
    	InitEllipticSpectrum(Matrix->n,density,ab_ratio,spectrum_file_name);
    }

    else if(spectrum_type == "Concentred") InitConcentredSpectrum(Matrix->n,density,spectrum_file_name);
    
    else if(spectrum_type == "Clustered") InitClusteredSpectrum(Matrix->n,density,spectrum_file_name);

    Matrix->A = GenerateMatrixFromSpectrum(Matrix->n,d,l,spectrum_file_name);

    InitVectorsRandomly(Matrix->A,&x,&u,&(Vector->RHS),Matrix->n);

    Matrix->SetNnz();
    Matrix->SetBandwidth();
    Matrix->SetMinMax();
    Matrix->SetSparsity();
    Matrix->SetDiagMinMax();
    Matrix->SetTrace();
    Matrix->SetDiagVariance();
    Vector->SetFeatures();

    if(!rank){
    	Matrix->PrintFeatures();
    	Vector->PrintFeatures();
   	}

    //GMRES method

   	VecSet(x,(PetscScalar)1.0/Matrix->n);

    Matrix->solver = 0;

    Solve(KSPGMRES,Matrix->A,Vector->RHS,Matrix->n,(PetscReal)tolerance,x,&(Matrix->its),&(Matrix->elapsed_time),&kspgmres);

    FinalError(x,u,NORM_2,&norm);

    if((double)norm/Matrix->n <= mean_norm2_tol) Matrix->converged = 1;
    else Matrix->converged = 0;
    
    PetscPrintf(PETSC_COMM_WORLD,"Norm of error : %g, Iterations : %D, Elapsed Time : %g\n",(double)norm/Matrix->n,(Matrix->its),(Matrix->elapsed_time));

    if(!rank) InsertNewLine(&file,Matrix,Vector);

    //Bi Conugate Gradient method

    VecSet(x,(PetscScalar)1.0/Matrix->n);

    Matrix->solver = 1;
    
    Solve(KSPBICG,Matrix->A,Vector->RHS,Matrix->n,tolerance,x,&(Matrix->its),&(Matrix->elapsed_time),&kspbicg);

    FinalError(x,u,NORM_2,&norm);

    if((double)norm/Matrix->n <= mean_norm2_tol) Matrix->converged = 1;
    else Matrix->converged = 0;
    
    PetscPrintf(PETSC_COMM_WORLD,"Norm of error : %g, Iterations : %D, Elapsed Time : %g\n",(double)norm/Matrix->n,(Matrix->its),(Matrix->elapsed_time));

    if(!rank) InsertNewLine(&file,Matrix,Vector);

    //Cleaning ... 

    file.close();

    VecDestroy(&x);VecDestroy(&Vector->RHS);VecDestroy(&u);

    MatDestroy(&Matrix->A);

    KSPDestroy(&kspgmres); KSPDestroy(&kspbicg);

   	PetscFinalize();
 
    return 0;
}
