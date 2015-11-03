#ifndef UTILITIES_H
#define UTILITIES_H

#include <admodel.h>
#include <fvar.hpp>

adstring stripExtension(adstring fileName);
ivector getIndex(const dvector& a, const dvector& b);
void write_proj_headers(ofstream &ofsP, int syr, int nyr);
void write_proj_output(ofstream &ofsP, int syr, int nyr, double tac, int pyr, dvector p_sbt, dmatrix p_ft, dvar_matrix ft, double bo, dmatrix fmsy, dvector bmsy);

#endif
