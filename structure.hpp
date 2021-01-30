#ifndef STRUCTURE_H
#define STRUCTURE_H

#include"utilities.hpp"

class VectorFeatures{

	public:
		Vec RHS;
		int n;
		double mean;
		double variance;
		double gravity_center;

		PetscReal min_mod;
		PetscReal max_mod;

		VectorFeatures();

		void PrintFeatures();
		void PrintFeatures(std::string *buffer);

		void SetFeatures();
};

inline VectorFeatures::VectorFeatures(){};

inline void VectorFeatures::PrintFeatures(){
	std::cout << "\nVector : " << endl;
	std::cout << "	-> Mean : " << mean << endl;
	std::cout << "	-> Variance : " << variance << endl;
	std::cout << "	-> Centre of gravity : " << gravity_center << endl;
	std::cout << "	-> Minimum module : " << min_mod << endl;
	std::cout << "	-> Maximum module : " << max_mod << "\n" << endl;
}

inline void VectorFeatures::PrintFeatures(std::string *buffer){
	*buffer += to_string(n) + "," + to_string(mean) + "," + to_string(variance) + "," + to_string(gravity_center) + "," + to_string(min_mod) + "," + to_string(max_mod); 
}

inline void VectorFeatures::SetFeatures(){
	
	int first, last; VecGetOwnershipRange(RHS,&first,&last);

	PetscScalar *v;

	VecGetArray(RHS,&v);

	PetscScalar local_min,local_max; local_min = local_max = v[0];

	PetscReal local_min_mod, local_max_mod; local_min_mod = local_max_mod = (PetscReal)sqrt(PetscSqr(PetscRealPart(local_max)) + PetscSqr(PetscImaginaryPart(local_max)));

	double local_mean = 0.0,local_variance = 0.0,local_mean_index = 0.0;

	PetscReal module_temp;

	for(int j = 1; j < last-first; j++){
		
		module_temp = (PetscReal)sqrt(PetscSqr(PetscRealPart(v[j])) + PetscSqr(PetscImaginaryPart(v[j])));

		local_mean += (double)module_temp;
		local_variance += (double)(module_temp*module_temp);
		local_mean_index += (double)(j+first)*module_temp;

		if(module_temp < local_min_mod){

			local_min = v[j];
			local_min_mod = module_temp;

		}else if(module_temp > local_max_mod){

			local_max = v[j];
			local_max_mod = module_temp;

		}
	}

	MPI_Allreduce(&local_max_mod,&max_mod,1,MPIU_REAL,MPI_MAX,PETSC_COMM_WORLD);
	MPI_Allreduce(&local_min_mod,&min_mod,1,MPIU_REAL,MPI_MIN,PETSC_COMM_WORLD);

	MPI_Allreduce(&local_mean,&mean,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
	MPI_Allreduce(&local_variance,&variance,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

	MPI_Allreduce(&local_mean_index,&gravity_center,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
	
	gravity_center/=(mean*n);

	mean /= n;
	variance = variance/n - mean*mean;

	VecRestoreArray(RHS,&v);

}

class MatrixFeatures{

	public:

		Mat A;

		//Input Features
		int n;
		int nnz;
		int ndr;

		double sparsity;
		double trace_real;
		double trace_imag;
		double amplitude_mean;
		double amplitude_variance;

		MatrixElement min;
		MatrixElement max;

		MatrixElement diag_min;
		MatrixElement diag_max;

		double diag_mean;
		double diag_variance;
		double diag_gravity;

		int lbandwidth;
		int rbandwidth;

		int solver; //0 for GMRES, 1 for BICG

		//Output Features

		int converged;
		int its;
		double elapsed_time;

		MatrixFeatures();

		void PrintFeatures();
		void PrintInputFeatures(std::string *buffer);
		void PrintOutputFeatures(std::string *buffer);

		void SetNnz();
		void SetTrace();
		void SetSparsity();
		void SetBandwidth();
		void SetMinMax();
		void SetMin();
		void SetMax();
		void SetDiagMinMax();
		void SetDiagVariance();
};

inline MatrixFeatures::MatrixFeatures(){};

inline void MatrixFeatures::PrintInputFeatures(std::string *buffer){
    *buffer += to_string(n) + "," + to_string(nnz) + "," + to_string(sparsity) + ",";
    *buffer += to_string(trace_real) + "," + to_string(trace_imag) + "," + to_string(amplitude_mean) + "," + to_string(amplitude_variance) + ",";
    *buffer += to_string(min.mod) + "," + to_string(max.mod) + "," + to_string(diag_min.mod) + "," + to_string(diag_max.mod) + ",";
    *buffer += to_string(diag_mean) + "," + to_string(diag_variance) + "," + to_string(diag_gravity) + ",";
    *buffer += to_string(lbandwidth) + "," + to_string(rbandwidth) + "," + to_string(solver);
}

inline void MatrixFeatures::PrintOutputFeatures(std::string *buffer){
	*buffer += to_string(converged) + "," + to_string(its) + "," + to_string(elapsed_time);
}

inline void MatrixFeatures::PrintFeatures(){
		std::cout << "\nMatrix " << n << "x" << n << " : " << endl;
	    std::cout << "	->Bandwidth : " << lbandwidth << "," << rbandwidth << endl;
    	std::cout << "	->Number of non-zero : " << nnz << endl;
    	std::cout << "	->Sparsity : " << sparsity << "%" << endl;
    	std::cout << "	->Number of dummy rows : " << ndr << endl;
    	std::cout << "	->Matrix's Minimum : " << min.val << ", Matrix's Maximum : " << max.val << endl;
    	std::cout << "	->Diag's Minimum : " << diag_min.val << ", Diag's Maximum : " << diag_max.val << endl;
    	std::cout << "	->Diag's Module (Mean) : " << diag_mean << ", Diag's Module (variance) : " << diag_variance << endl;
    	std::cout << "	->Trace (Real) : " << trace_real << ", Trace(Imaginary) : " << trace_imag << endl;
    	std::cout << "	->Amplitude (Mean) : " << amplitude_mean << ", Entropy : " << amplitude_variance << endl;
    	std::cout << "	->Diag Centre of gravity : " << diag_gravity  << "\n" << endl;
}

inline void MatrixFeatures::SetNnz(){

	int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	int size; MPI_Comm_size(PETSC_COMM_WORLD,&size);

	int first,last; MatGetOwnershipRange(A,&first,&last);

	//std::cout << rank << " : " << first << "," << last << endl;

	nnz = 0; this->ndr = 0; PetscInt row_nnz; const PetscInt *row_cols; const PetscScalar *row_vals;

	for(int i = first; i < last; i++){

		MatGetRow(A,i,&row_nnz,&row_cols,&row_vals);

		//std::cout << rank << ": " << i << "-> " << row_nnz << endl;

		if(row_nnz < 2) this->ndr += 1;

		this->nnz += row_nnz;

		MatRestoreRow(A,i,&row_nnz,&row_cols,&row_vals);

	}

	MPI_Allreduce(&nnz,&row_nnz,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);

	this->nnz = row_nnz;
}

inline void MatrixFeatures::SetSparsity(){
	this->sparsity = 100-(nnz*100.0/(n*n));
}

inline void MatrixFeatures::SetBandwidth(){

	int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	int size; MPI_Comm_size(PETSC_COMM_WORLD,&size);

	int first,last; MatGetOwnershipRange(A,&first,&last);
	int ldepth_max = 0,rdepth_max = 0;

	PetscInt row_nnz; const PetscInt *row_cols; const PetscScalar *row_vals;

	for(int i = first; i < last; i++){

		MatGetRow(A,i,&row_nnz,&row_cols,&row_vals);

		if(row_cols[0] < i         && (i-row_cols[0]) > ldepth_max)              ldepth_max = i-row_cols[0];
		if(row_cols[row_nnz-1] > i && (row_cols[row_nnz-1]-i) > rdepth_max) rdepth_max = row_cols[row_nnz-1]-i;

		MatRestoreRow(A,i,&row_nnz,&row_cols,&row_vals);
	}

	//printf("(%d) -> (%d,%d) | Left : %d , Right : %d\n",rank,first,last,ldepth_max,rdepth_max);

	MPI_Allreduce(&ldepth_max,&lbandwidth,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
	MPI_Allreduce(&rdepth_max,&rbandwidth,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
}

inline void MatrixFeatures::SetMin(){

	int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	int size; MPI_Comm_size(PETSC_COMM_WORLD,&size);

	int first,last; MatGetOwnershipRange(A,&first,&last);

	MatrixElement local_min,min;

	PetscInt row_nnz; const PetscInt *row_cols; const PetscScalar *row_vals;

	MPI_Datatype MyMPI_MatrixElement = InitMatrixElementMPIDataType(); MPI_Op MyMPI_MinMatrixReduction = InitMinMatrixMPIReduction();

	MatGetRow(A,first,&row_nnz,&row_cols,&row_vals);

	MinRow(row_nnz,row_cols,row_vals,&local_min);

	local_min.row = first;

	MatRestoreRow(A,first,&row_nnz,&row_cols,&row_vals);

	for(int i = first + 1; i < last; i++){

		MatGetRow(A,i,&row_nnz,&row_cols,&row_vals);

		MinRow(row_nnz,row_cols,row_vals,&min);

		if(min.mod < local_min.mod){

			local_min.mod = min.mod;
			local_min.val = min.val;
			local_min.col = min.col;

			local_min.row = i;
		}

		MatRestoreRow(A,i,&row_nnz,&row_cols,&row_vals);
	}

	MPI_Allreduce(&local_min,&min,1,MyMPI_MatrixElement,MyMPI_MinMatrixReduction,PETSC_COMM_WORLD);

	CleanMPIDatatype(&MyMPI_MatrixElement); CleanMPIOperation(&MyMPI_MinMatrixReduction);
}



inline void MatrixFeatures::SetMax(){

	int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	int size; MPI_Comm_size(PETSC_COMM_WORLD,&size);

	int first,last; MatGetOwnershipRange(A,&first,&last);

	MatrixElement local_max,max;

	PetscInt row_nnz; const PetscInt *row_cols; const PetscScalar *row_vals;

	MPI_Datatype MyMPI_MatrixElement = InitMatrixElementMPIDataType(); MPI_Op MyMPI_MaxMatrixReduction = InitMaxMatrixMPIReduction();

	MatGetRow(A,first,&row_nnz,&row_cols,&row_vals);

	MaxRow(row_nnz,row_cols,row_vals,&local_max);

	local_max.row = first;

	MatRestoreRow(A,first,&row_nnz,&row_cols,&row_vals);

	for(int i = first + 1; i < last; i++){

		MatGetRow(A,i,&row_nnz,&row_cols,&row_vals);

		MaxRow(row_nnz,row_cols,row_vals,&max);

		if(max.mod > local_max.mod){

			local_max.mod = max.mod;
			local_max.val = max.val;
			local_max.col = max.col;

			local_max.row = i;
		}

		MatRestoreRow(A,i,&row_nnz,&row_cols,&row_vals);
	}

	MPI_Allreduce(&local_max,&max,1,MyMPI_MatrixElement,MyMPI_MaxMatrixReduction,PETSC_COMM_WORLD);

	CleanMPIDatatype(&MyMPI_MatrixElement); CleanMPIOperation(&MyMPI_MaxMatrixReduction);
}

inline void MatrixFeatures::SetMinMax(){

	int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	int size; MPI_Comm_size(PETSC_COMM_WORLD,&size);

	int first,last; MatGetOwnershipRange(A,&first,&last);

	double local_amplitude_sum = 0, local_amplitude_variance = 0; PetscScalar diff;

	MatrixElement local_max,local_min;

	PetscInt row_nnz; const PetscInt *row_cols; const PetscScalar *row_vals;

	MPI_Datatype MyMPI_MatrixElement = InitMatrixElementMPIDataType();

	MPI_Op MyMPI_MinMatrixReduction = InitMinMatrixMPIReduction();
	MPI_Op MyMPI_MaxMatrixReduction = InitMaxMatrixMPIReduction();

	MatGetRow(A,first,&row_nnz,&row_cols,&row_vals);

	MinMaxRow(row_nnz,row_cols,row_vals,&local_min,&local_max);

	local_min.row = first; local_max.row = first;

	MatRestoreRow(A,first,&row_nnz,&row_cols,&row_vals);

	for(int i = first + 1; i < last; i++){

		MatGetRow(A,i,&row_nnz,&row_cols,&row_vals);

		MinMaxRow(row_nnz,row_cols,row_vals,&min,&max);

		diff = min.val - max.val;

		local_amplitude_sum += sqrt(PetscSqr(PetscRealPart(diff)) + PetscSqr(PetscImaginaryPart(diff)));
		local_amplitude_variance += PetscSqr(PetscRealPart(diff)) + PetscSqr(PetscImaginaryPart(diff));

		if(min.mod < local_min.mod){

			local_min.mod = min.mod;
			local_min.val = min.val;
			local_min.col = min.col;

			local_max.row = i;

		}else if(max.mod > local_max.mod){

			local_max.mod = max.mod;
			local_max.val = max.val;
			local_max.col = max.col;

			local_max.row = i;
		}

		MatRestoreRow(A,i,&row_nnz,&row_cols,&row_vals);
	}

	MPI_Allreduce(&local_min,&min,1,MyMPI_MatrixElement,MyMPI_MinMatrixReduction,PETSC_COMM_WORLD);
	MPI_Allreduce(&local_max,&max,1,MyMPI_MatrixElement,MyMPI_MaxMatrixReduction,PETSC_COMM_WORLD);
	MPI_Allreduce(&local_amplitude_sum,&amplitude_mean,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
	MPI_Allreduce(&local_amplitude_variance,&amplitude_variance,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

	amplitude_mean /= n;
	amplitude_variance = amplitude_variance/n - amplitude_mean*amplitude_mean;

	CleanMPIDatatype(&MyMPI_MatrixElement); 

	CleanMPIOperation(&MyMPI_MinMatrixReduction); CleanMPIOperation(&MyMPI_MaxMatrixReduction);
}

inline void MatrixFeatures::SetDiagMinMax(){

	int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	int size; MPI_Comm_size(PETSC_COMM_WORLD,&size);

	int first,last,j; MatGetOwnershipRange(A,&first,&last);

	MatrixElement local_max,local_min,diag;

	PetscInt row_nnz; const PetscInt *row_cols; const PetscScalar *row_vals;

	MPI_Datatype MyMPI_MatrixElement = InitMatrixElementMPIDataType();

	MPI_Op MyMPI_MinMatrixReduction = InitMinMatrixMPIReduction();
	MPI_Op MyMPI_MaxMatrixReduction = InitMaxMatrixMPIReduction();

	local_min.mod = PETSC_INFINITY; local_max.mod = 0;

	for(int i = first; i < last; i++){

		MatGetRow(A,i,&row_nnz,&row_cols,&row_vals);

		j = 0;

		while(j < row_nnz && row_cols[j] < i){
			j++;
		}

		if(row_cols[j] == i){

			diag.val = row_vals[j];
			diag.mod = (PetscReal)sqrt(PetscSqr(PetscRealPart(diag.val)) + PetscSqr(PetscImaginaryPart(diag.val)));
			diag.row = i;
			diag.col = i;

			if(diag.mod < local_min.mod){

				local_min.mod = diag.mod;
				local_min.val = diag.val;
				local_min.col = diag.col;

				local_max.row = i;

			}if(diag.mod > local_max.mod){

				local_max.mod = diag.mod;
				local_max.val = diag.val;
				local_max.col = diag.col;

				local_max.row = i;
			}
		}
		MatRestoreRow(A,i,&row_nnz,&row_cols,&row_vals);
	}

	MPI_Allreduce(&local_min,&diag_min,1,MyMPI_MatrixElement,MyMPI_MinMatrixReduction,PETSC_COMM_WORLD);
	MPI_Allreduce(&local_max,&diag_max,1,MyMPI_MatrixElement,MyMPI_MaxMatrixReduction,PETSC_COMM_WORLD);

	CleanMPIDatatype(&MyMPI_MatrixElement); 

	CleanMPIOperation(&MyMPI_MinMatrixReduction); CleanMPIOperation(&MyMPI_MaxMatrixReduction);
}

inline void MatrixFeatures::SetTrace(){

	int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	int size; MPI_Comm_size(PETSC_COMM_WORLD,&size);

	int first,last,j; MatGetOwnershipRange(A,&first,&last);

	PetscInt row_nnz; const PetscInt *row_cols; const PetscScalar *row_vals;

	PetscScalar trace,local_trace = (PetscScalar)0.0;

	for(int i = first; i < last; i++){
		
		MatGetRow(A,i,&row_nnz,&row_cols,&row_vals);

		j = 0;

		while(j < row_nnz && row_cols[j] < i){
			j++;
		}

		if(row_cols[j] == i){
			local_trace += row_vals[j];
		}

		MatRestoreRow(A,i,&row_nnz,&row_cols,&row_vals);
	}


	MPI_Allreduce(&local_trace,&trace,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);
	trace_real = (double)PetscRealPart(trace);
	trace_imag = (double)PetscImaginaryPart(trace);
}

inline void MatrixFeatures::SetDiagVariance(){

	int rank; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	int size; MPI_Comm_size(PETSC_COMM_WORLD,&size);

	int first,last,j; MatGetOwnershipRange(A,&first,&last);

	PetscInt row_nnz; const PetscInt *row_cols; const PetscScalar *row_vals;

	PetscScalar trace,local_trace = (PetscScalar)0.0;

	double local_mean, local_variance,local_gravity;

	for(int i = first; i < last; i++){
		
		MatGetRow(A,i,&row_nnz,&row_cols,&row_vals);

		j = 0;

		while(j < row_nnz && row_cols[j] < i){
			j++;
		}

		if(row_cols[j] == i){
			local_mean +=  sqrt(PetscSqr(PetscRealPart(row_vals[j])) + PetscSqr(PetscImaginaryPart(row_vals[j])));
			local_variance += PetscSqr(PetscRealPart(row_vals[j])) + PetscSqr(PetscImaginaryPart(row_vals[j]));
			local_gravity += i*sqrt(PetscSqr(PetscRealPart(row_vals[j])) + PetscSqr(PetscImaginaryPart(row_vals[j])));
		}

		MatRestoreRow(A,i,&row_nnz,&row_cols,&row_vals);
	}


	MPI_Allreduce(&local_mean,&diag_mean,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
	MPI_Allreduce(&local_variance,&diag_variance,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
	MPI_Allreduce(&local_gravity,&diag_gravity,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);

	diag_gravity /= (n*diag_mean);
	diag_mean /= n;
	diag_variance = (diag_variance)/n - diag_mean*diag_mean;

}


#endif