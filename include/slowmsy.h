#ifndef __SLOW_MSY_H
#define __SLOW_MSY_H

#include <admodel.h>
void slow_msy(dvector& ftest,
              dvector& ye,
              dvector& be,
              double& msy,
              double& fmsy,
              double& bmsy,
              double& bo,
              int sage,
              int nage,
              int nyr,
              dvar3_array M,
              dmatrix dWt_bar,
              dmatrix ma,
              dvar_vector ro,
              dvar_vector kappa,
              dvar4_array log_sel,
              dvector d_iscamCntrl,
              dvector pf_cntrl);
#endif
