#pragma once 
#include <cmath>

namespace numcpp::erf {

    inline double dint(const double x) { return( (x>0) ? floor(x) : -floor(-x) ); }

    inline double calErf(double x, const int jint) {

        static const double a[5] = { 3.1611237438705656,113.864154151050156,377.485237685302021,3209.37758913846947,.185777706184603153 };
        static const double b[4] = { 23.6012909523441209,244.024637934444173,1282.61652607737228,2844.23683343917062 };
        static const double c__[9] = { .564188496988670089,8.88314979438837594,66.1191906371416295,298.635138197400131,881.95222124176909,1712.04761263407058,2051.07837782607147,1230.33935479799725,2.15311535474403846e-8 };
        static const double d__[8] = { 15.7449261107098347,117.693950891312499,537.181101862009858,1621.38957456669019,3290.79923573345963,4362.61909014324716,3439.36767414372164,1230.33935480374942 };
        static const double p[6] = { .305326634961232344,.360344899949804439,.125781726111229246,.0160837851487422766,6.58749161529837803e-4,.0163153871373020978 };
        static const double q[5] = { 2.56852019228982242,1.87295284992346047,.527905102951428412,.0605183413124413191,.00233520497626869185 };

        static const double zero = 0.;
        static const double half = .5;
        static const double one = 1.;
        static const double two = 2.;
        static const double four = 4.;
        static const double sqrpi = 0.56418958354775628695;
        static const double thresh = .46875;
        static const double sixten = 16.;

        double y, del, ysq, xden, xnum, result;

        static const double xinf = 1.79e308;
        static const double xneg = -26.628;
        static const double xsmall = 1.11e-16;
        static const double xbig = 26.543;
        static const double xhuge = 6.71e7;
        static const double xmax = 2.53e307;
        
        y = std::fabs(x);
        if (y <= thresh) {

            ysq = zero;

            if (y > xsmall) {
                ysq = y * y;
            }

            xnum = a[4] * ysq;

            xden = ysq;

            for (int i__ = 1; i__ <= 3; ++i__) {
                xnum = (xnum + a[i__ - 1]) * ysq;
                xden = (xden + b[i__ - 1]) * ysq;
            }
            result = x * (xnum + a[3]) / (xden + b[3]);
            if (jint != 0) {
                result = one - result;
            }
            if (jint == 2) {
                result = exp(ysq) * result;
            }
            goto L800;

        } else if (y <= four) {
            xnum = c__[8] * y;
            xden = y;
            for (int i__ = 1; i__ <= 7; ++i__) {
                xnum = (xnum + c__[i__ - 1]) * y;
                xden = (xden + d__[i__ - 1]) * y;
            }
            result = (xnum + c__[7]) / (xden + d__[7]);
            if (jint != 2) {
                double d__1 = y * sixten;
                ysq = dint(d__1) / sixten;
                del = (y - ysq) * (y + ysq);
                d__1 = exp(-ysq * ysq) * exp(-del);
                result = d__1 * result;
            }

        } else {
            result = zero;
            if (y >= xbig) {
                if (jint != 2 || y >= xmax) {
                    goto L300;
                }
                if (y >= xhuge) {
                    result = sqrpi / y;
                    goto L300;
                }
            }
            ysq = one / (y * y);
            xnum = p[5] * ysq;
            xden = ysq;
            for (int i__ = 1; i__ <= 4; ++i__) {
                xnum = (xnum + p[i__ - 1]) * ysq;
                xden = (xden + q[i__ - 1]) * ysq;
            }
            result = ysq * (xnum + p[4]) / (xden + q[4]);
            result = (sqrpi - result) / y;
            if (jint != 2) {
                double d__1 = y * sixten;
                ysq = dint(d__1) / sixten;
                del = (y - ysq) * (y + ysq);
                d__1 = exp(-ysq * ysq) * exp(-del);
                result = d__1 * result;
            }
        }
        L300:
        if (jint == 0) {
            result = (half - result) + half;
            if (x < zero) {
                result = -(result);
            }
        } else if (jint == 1) {
            if (x < zero) {
                result = two - result;
            }
        } else {
            if (x < zero) {
                if (x < xneg) {
                    result = xinf;
                } else {
                    double d__1 = x * sixten;
                    ysq = dint(d__1) / sixten;
                    del = (x - ysq) * (x + ysq);
                    y = exp(ysq * ysq) * exp(del);
                    result = y + y - result;
                }
            }
        }
        L800:
        return result;
    }

    inline double erf(double x) {return calErf(x, 0);}

    inline double erfc(double x) {return calErf(x, 1);}

    inline double erfcx(double x) {return calErf(x, 2);}

}