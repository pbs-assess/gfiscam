#ifndef __UTILITIES_H
#define __UTILITIES_H

#include <admodel.h>

adstring stripExtension(adstring fileName);

ivector getIndex(const dvector& a, const dvector& b);

void write_proj_headers(ofstream &ofsP,
                        int syr,
                        int nyr,
                        bool include_msy,
                        bool include_sbo);

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
                       dmatrix ma,
                       dmatrix dWt_bar,
                       dvar_matrix ft,
                       double bo,
                       dmatrix fmsy,
                       dvector bmsy,
                       bool include_msy,
                       bool include_bo);

void write_proj_headers_dd(ofstream &ofsP,
                           int nyr,
                           int pyr);

void write_proj_output_dd(ofstream &ofsP,
                          double tac,
                          int nyr,
                          int pyr,
                          dvector p_bt,
                          dvector p_ft,
                          dmatrix fmsy,
                          dvector bmsy,
                          double bmin,
                          double meanbshort,
                          double meanblong,
                          double meanfshort,
                          double meanflong);

#endif
