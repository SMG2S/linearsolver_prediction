#ifndef SOLVER_H
#define SOLVER_H

#include"structure.hpp"

#include"utilities.hpp"

inline void InitClassicSpectrum(int n,float density,std::string spectrum_file_path){

	int rank;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	if(!rank){
		ClassicSpectrum *Sigma = new ClassicSpectrum(n);
    	Sigma->setDensity(density);
    	Sigma->randomPattern();
    	Sigma->printParameters();
    	Sigma->printSpectrum(spectrum_file_path);

    	Sigma->freeSpectrum();
    	delete(Sigma);
    }
}

inline void InitEllipticSpectrum(int n,float density,float ab_ratio,std::string spectrum_file_path){

	int rank;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	if(!rank){
		EllipticSpectrum *Sigma = new EllipticSpectrum(n);
    	Sigma->setDensity(density);
    	Sigma->setAbRatio(ab_ratio);
    	Sigma->randomPattern();
    	Sigma->printParameters();
    	Sigma->printSpectrum(spectrum_file_path);

    	Sigma->freeSpectrum();
    	delete(Sigma);
    }
}

inline void InitConcentredSpectrum(int n,float density,std::string spectrum_file_path){

	int rank;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	if(!rank){
		ConcentredSpectrum *Sigma = new ConcentredSpectrum(n);
    	Sigma->setDensity(density);
    	Sigma->randomPattern();
    	Sigma->printParameters();
    	Sigma->printSpectrum(spectrum_file_path);

    	Sigma->freeSpectrum();
    	delete(Sigma);
    }
}

inline void InitClusteredSpectrum(int n,float density,std::string spectrum_file_path){

	int rank;

	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	if(!rank){
		std::cout << "ok" << endl;
		ClusteredSpectrum *Sigma = new ClusteredSpectrum(n,rand()%5 + 2);
    	//Sigma->setDensity(density);
    	Sigma->randomPattern();
    	Sigma->printParameters();
    	Sigma->printSpectrum(spectrum_file_path);

    	Sigma->freeSpectrum();
    	delete(Sigma);
    }
}


inline Mat GenerateMatrixFromSpectrum(PetscInt n,PetscInt d,PetscInt lbandwidth,std::string spectrum_file_path){

	std::string spectrum;
	spectrum.assign(spectrum_file_path);

	int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	int size; MPI_Comm_size(PETSC_COMM_WORLD,&size);

	int first,last; GetBounds(rank,size,n,&first,&last);
    
    Nilpotency<int> nilp;
	nilp.NilpType1(d,n);
    
    parMatrixSparse<std::complex<double>,int> *Mt;

    Mt = smg2s<std::complex<double>,int>(n,nilp,lbandwidth,spectrum,PETSC_COMM_WORLD);

	Mt->Loc_ConvertToCSR();

	Mat A;

	MatCreate(PETSC_COMM_WORLD,&A);

	MatSetSizes(A,last-first+1,n,n,n);

	A = ConvertToPETSCMat(Mt);
    
    return A;
}

inline PetscErrorCode InitVectorsRandomly(Mat A, Vec *x, Vec *u, Vec *b, int n){
	
	PetscErrorCode ierr;

	int local_m,local_n;

	MatGetLocalSize(A,&local_m,&local_n);

	//std::cout << local_m << "," << local_n << endl; 

	PetscRandom rctx; 

    ierr = VecCreate(PETSC_COMM_WORLD,x);

    ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
    PetscRandomSetType(rctx,PETSCRAND48);

    ierr = PetscObjectSetName((PetscObject)(*x),"x");
    ierr = VecSetSizes(*x,local_n,n);
    ierr = VecSetFromOptions(*x);
    ierr = VecDuplicate(*x,b);
    ierr = PetscObjectSetName((PetscObject)(*b),"RHS");

    ierr = VecCreate(PETSC_COMM_WORLD,u);
    ierr = VecSetSizes(*u,local_n,n);
	ierr = VecSetFromOptions(*u);

    ierr = VecSetRandom(*u,PETSC_NULL);

    ierr = PetscObjectSetName((PetscObject)(*u),"Solution");

    //VecView(*b,PETSC_VIEWER_STDOUT_WORLD);

    ierr = MatMult(A,*u,*b);

    ierr = PetscRandomDestroy(&rctx);

    return ierr;
}

inline PetscErrorCode Solve(KSPType ksptype,Mat A,Vec b,PetscInt n,PetscReal tol, Vec x,PetscInt *its, PetscReal *time,KSP *ksp){
	
	PetscErrorCode ierr;

	PetscReal begin = MPI_Wtime();

	ierr = KSPCreate(PETSC_COMM_WORLD,ksp);
	ierr = KSPSetType(*ksp,ksptype);
    ierr = KSPSetOperators(*ksp,A,A);
    ierr = KSPSetInitialGuessNonzero(*ksp,PETSC_TRUE);

    ierr = KSPSetTolerances(*ksp,PETSC_DEFAULT,tol,PETSC_DEFAULT,1000);
    ierr = KSPSetNormType(*ksp,KSP_NORM_UNPRECONDITIONED);
    ierr = KSPSetFromOptions(*ksp);
    ierr = KSPSetUp(*ksp);
    ierr = KSPSolve(*ksp,b,x);

    PetscReal end = MPI_Wtime();

    *time = end-begin;

    ierr = KSPGetIterationNumber(*ksp,its);

    return ierr;
}

inline PetscErrorCode FinalError(Vec x, Vec u, NormType p, PetscReal *error){
	
	PetscErrorCode ierr;
	Vec residual_vec;

	ierr = VecDuplicate(x,&residual_vec);
	ierr = VecCopy(x,residual_vec);

	ierr = VecAXPY(residual_vec,-1.0,u);
	ierr = VecNorm(residual_vec,p,error);

	ierr = VecDestroy(&residual_vec);

	return ierr;
}

#endif
