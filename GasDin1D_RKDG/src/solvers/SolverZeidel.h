#pragma once

#include "MatrixSolver.h"
#include <vector>
#include <list>


class SolverZeidel : public MatrixSolver
{
public:
	//virtual void init(int cellsCount, int blockDimension);

	virtual int solve(double eps, int& maxIter);
	virtual char* getName() { return "ZEIDEL"; }
	//virtual void zero();
	//virtual void setMatrElement(int i, int j, double** matrDim);
	//virtual void addMatrElement(int i, int j, double** matrDim);




private:
	//
	//typedef std::pair<int, double> Element;
	//typedef std::vector< Element > Row;

	//struct sort_row_class
	//{
	//	bool operator() (Element i, Element j)
	//	{
	//		return (i.first < j.first);
	//	}
	//} sort_row;	
	//
	//int N;
	//Row *rows;
	//void set(int i, int j, double a);
	//void add(int i, int j, double a);
	//void assemble();
};

