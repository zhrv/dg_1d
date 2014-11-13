#ifndef _CSR_
#define _CSR_

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <map>

typedef std::map<int, double> Row;

struct CSRMatrix
{
	CSRMatrix(int N);
	~CSRMatrix();

	void	zero();
	void	set(int i, int j, double aa);
	double	get(int i, int j);
	void	add(int i, int j, double aa);
	void	printToFile(const char *fileName);
	void	init(int i, int j);
	void	assemble();

	double *a;
	int	   *ia;
	int    *ja;
	int     na;
	int     _na;
	int     n;

	static int DELTA;

	Row		*rows;
};

#endif