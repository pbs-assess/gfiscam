#ifndef __UTILITIES_H
#define __UTILITIES_H

#include <admodel.h>
adstring stripExtension(adstring fileName);
ivector getIndex(const dvector& a, const dvector& b);
void write_proj_headers(ofstream &ofsP, int syr, int nyr, bool include_msy);
void write_proj_output(ofstream &ofsP, int syr, int nyr, double tac, int pyr, dvector p_sbt, dmatrix p_ft, dvar_matrix ft, double bo, dmatrix fmsy, dvector bmsy, bool include_msy);

#endif
