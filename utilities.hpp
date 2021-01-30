#ifndef UTILITIES_H
#define UTILITIES_H

#include"string"
#include"petscksp.h"
#include"petsc_interface.h"

#include"../spectrum_generator/generator.hpp"


#include<mpi.h>

struct MatrixElement {
	PetscScalar val;
	PetscReal   mod;
	PetscInt    col;
	PetscInt    row;
};

inline void AddSep(std::string *buffer,char sep){
	*buffer += sep;
}

inline void ReduceMaxMatrix(void *in, void *out, int *len, int *dptr){

	MatrixElement *input  = (MatrixElement *)(in);
	MatrixElement *output = (MatrixElement *)(out);
	for(int i = 0; i < *len; i++){
		if(input[i].mod > output[i].mod){
			output[i].mod = input[i].mod;
			output[i].val = input[i].val;
			output[i].col = input[i].col;
			output[i].row = input[i].row;
		}
	}
}

inline void ReduceMinMatrix(void *in, void *out, int *len, int *dptr){

	MatrixElement *input  = (MatrixElement *)(in);
	MatrixElement *output = (MatrixElement *)(out);
	for(int i = 0; i < *len; i++){
		if(input[i].mod < output[i].mod){
			output[i].mod = input[i].mod;
			output[i].val = input[i].val;
			output[i].col = input[i].col;
			output[i].row = input[i].row;
		}
	}
}

inline MPI_Datatype InitMatrixElementMPIDataType(){

	MPI_Datatype MPI_MatrixElement;

	int block_lengths[4];
	MPI_Aint displacements[4];
	MPI_Datatype typelist[4] = {MPIU_SCALAR, MPIU_REAL, MPIU_INT, MPIU_INT};

	block_lengths[0] = 1; displacements[0] = (MPI_Aint)offsetof(struct MatrixElement, val);
	block_lengths[1] = 1; displacements[1] = (MPI_Aint)offsetof(struct MatrixElement, mod);
	block_lengths[2] = 1; displacements[2] = (MPI_Aint)offsetof(struct MatrixElement, col);
	block_lengths[3] = 1; displacements[3] = (MPI_Aint)offsetof(struct MatrixElement, row);

	MPI_Type_create_struct(4, block_lengths, displacements, typelist,&MPI_MatrixElement);
    MPI_Type_commit(&MPI_MatrixElement);

    return MPI_MatrixElement;
}

inline void CleanMPIDatatype(MPI_Datatype *data){
	MPI_Type_free(data);
}


inline MPI_Op InitMinMatrixMPIReduction(){

	MPI_Op MyMPI_Reduction;

	MPI_Op_create(&ReduceMinMatrix,1,&MyMPI_Reduction);

    return MyMPI_Reduction;
}

inline MPI_Op InitMaxMatrixMPIReduction(){

	MPI_Op MyMPI_Reduction;

	MPI_Op_create(&ReduceMaxMatrix,1,&MyMPI_Reduction);

    return MyMPI_Reduction;
}

inline void CleanMPIOperation(MPI_Op *operation){
	MPI_Op_free(operation);
}

inline void GetBounds(int rank, int size, int n, int *first, int *last){

	*last = n/size,first;

	*first = rank*(*last);

	if(rank < n%size){

		(*last)  ++;
		(*first) += rank;

	} else (*first) += n%size;

	(*last) += (*first)-1;
}

inline void MinRow(PetscInt row_nnz, const PetscInt *row_cols, const PetscScalar *row_vals, MatrixElement *m){
	
	m->col = row_cols[0]; m->val = row_vals[0]; m->mod = (PetscReal)sqrt(PetscSqr(PetscRealPart(m->val)) + PetscSqr(PetscImaginaryPart(m->val)));
	
	PetscReal module_temp;

	for(int j = 1; j < row_nnz; j++){
		module_temp = (PetscReal)sqrt(PetscSqr(PetscRealPart(row_vals[j])) + PetscSqr(PetscImaginaryPart(row_vals[j])));
		if(module_temp < m->mod){
			m->mod = module_temp;
			m->val = row_vals[j];
			m->col = row_cols[j];
		}
	}
}

inline void MaxRow(PetscInt row_nnz, const PetscInt *row_cols, const PetscScalar *row_vals, MatrixElement *m){
	
	m->col = row_cols[0]; m->val = row_vals[0]; m->mod = (PetscReal)sqrt(PetscSqr(PetscRealPart(m->val)) + PetscSqr(PetscImaginaryPart(m->val)));
	
	PetscReal module_temp;

	for(int j = 1; j < row_nnz; j++){
		module_temp = (PetscReal)sqrt(PetscSqr(PetscRealPart(row_vals[j])) + PetscSqr(PetscImaginaryPart(row_vals[j])));
		if(module_temp > m->mod){
			m->mod = module_temp;
			m->val = row_vals[j];
			m->col = row_cols[j];
		}
	}
}

inline void MinMaxRow(PetscInt row_nnz, const PetscInt *row_cols, const PetscScalar *row_vals, MatrixElement *min, MatrixElement *max){
	
	max->col = row_cols[0]; max->val = row_vals[0]; max->mod = (PetscReal)sqrt(PetscSqr(PetscRealPart(max->val)) + PetscSqr(PetscImaginaryPart(max->val)));

	min->col = max->col; min->val = max->val; min->mod = max->mod;
	
	PetscReal module_temp;

	for(int j = 1; j < row_nnz; j++){
		module_temp = (PetscReal)sqrt(PetscSqr(PetscRealPart(row_vals[j])) + PetscSqr(PetscImaginaryPart(row_vals[j])));
		if(module_temp < min->mod){
			min->mod = module_temp;
			min->val = row_vals[j];
			min->col = row_cols[j];
		}else if(module_temp > max->mod){
			max->mod = module_temp;
			max->val = row_vals[j];
			max->col = row_cols[j];
		}
	}
}


#endif