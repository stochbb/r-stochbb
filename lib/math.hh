#ifndef __SBB_MATH_HH__
#define __SBB_MATH_HH__

#include <cmath>

namespace sbb {

/** Logarithm of the gamma function. */
double lgamma(double z);
/** The error function. */
double erfc(double x);
/** The lower-incomplete gamma function. */
double gamma_li(double s, double z);
/** The upper-incomplete gamma function. */
double gamma_ui(double s, double z);
/** The beta function. */
double betai(double a, double b, double x);

}

#endif // __SBB_MATH_HH__
