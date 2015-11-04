#ifndef __LOGISTIC_STUDENT_T_H
#define __LOGISTIC_STUDENT_T_H

#include <admodel.h>
#include "LogisticNormal.h"

class logistic_student_t: public logistic_normal
{
private:
	dvariable m_v;  // degrees of freedom.
	dvariable negative_log_likelihood();

public:
	logistic_student_t(const dmatrix& _O,const dvar_matrix& _E,
	                   const double _minp,const double _eps=0);

	dvariable operator() ();
	dvariable operator() (const dvariable& _df);

};

#endif
