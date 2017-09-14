#ifndef __UTILITIES_H
#define __UTILITIES_H

#include <admodel.h>

adstring stripExtension(adstring fileName);
ivector getIndex(const dvector& a, const dvector& b);
void write_proj_headers(ofstream &ofsP,
                        int syr,
                        int nyr,
                        bool include_msy,
                        bool include_bo);
void write_proj_output(ofstream &ofsP,
                       int syr,
                       int nyr,
                       int nage,
                       double tac,
                       int pyr,
                       dvector p_sbt,
                       dmatrix p_ft,
                       dmatrix p_N,
                       dvar3_array M,
                       dmatrix dWt_bar,
                       dvar_matrix ft,
                       double bo,
                       dmatrix fmsy,
                       dvector bmsy,
                       bool include_msy,
                       bool include_bo);

#endif
