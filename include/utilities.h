
#ifndef __UTILITIES_H
#define __UTILITIES_H

#include <admodel.h>

adstring stripExtension(adstring fileName);

ivector getIndex(const dvector& a, const dvector& b);

void write_proj_headers(ofstream &ofsP,
                        int syr,
                        int nyr,
                        int nfleet,
                        int n_ags,
                        int ngroup,
                        bool include_msy,
                        bool include_sbo);

void write_proj_output(ofstream &ofsP,
                       int syr,
                       int nyr,
                       int nage,
                       int nfleet,
                       int n_ags,
                       int ngroup,
                       double tac,
                       int pyr,
                       dvector p_sbt,
                       dvector p_rt,
                       d3_array p_ft,
                       d3_array p_N,
                       dvar3_array M,
                       dmatrix ma,
                       dmatrix dWt_bar,
                       dvar3_array ft,
                       double sbo,
                       dmatrix fmsy,
                       dvector bmsy,
                       bool include_msy,
                       bool include_sbo);

#endif
