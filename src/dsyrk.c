#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/libextern.h>

#include "dsyrk.h"

// computes a * transpose(a) for the symmetric, nXn matrix a
void aat(double * x, int * n, double * result) {
	double alpha = 1.0;
	double beta = 0.0;
	dsyrk_("L", "N", n, n, &alpha, x, n, &beta, result, n);
}

// int main(int argc, char *argv[]) { return 0;}

// just compile, no link:
// gcc -I"\program files\r\r-2.13.1\include" -c dsyrk.c
// compile AND link.  Doesn't work with MinGW's gcc if 'i386' is replaced by 'x64'
// gcc -I"\Program Files\R\R-2.13.1\include" -L"\Program Files\R\R-2.13.1\bin\i386" -lRblas -o dsyrk.exe dsyrk.c