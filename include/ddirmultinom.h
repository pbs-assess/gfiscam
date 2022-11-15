/**
 * \brief Class for implementing the Dirichlet Multinomial negative log-likelihood.
 * \author Chris Grandin
 * \addtogroup Likelihoods [group-title]
 */

/**
Dirichlet Multinomial likelihood for size-age composition data.

Author: Chris Grandin
Date  : Sept 10, 2021

Descrition:
**/

#ifndef __DIRICHLET_MULTINOMIAL_H
#define __DIRICHLET_MULTINOMIAL_H

#include <admodel.h>

dvariable ddirmultinom(const dvar_vector& obs,
		       const dvar_vector& p,
		       const dvariable& log_phi,
		       const dvariable& samp_size,
		       bool linear);


#endif
