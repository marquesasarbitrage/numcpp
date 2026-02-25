#pragma once
#include <limits>

namespace numcpp {

    namespace constants {

        constexpr double ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399460599343818684758586311649;    
        
        constexpr double PI = 3.14159265358979323846;

        constexpr double SQRT_TWO_PI = 2.50662827463100050242;

        constexpr double ONE_OVER_SQRT_TWO = 0.7071067811865475244008443621048490392848359376887;

        constexpr double DOUBLE_NAN = std::numeric_limits<double>::quiet_NaN();

    }
}