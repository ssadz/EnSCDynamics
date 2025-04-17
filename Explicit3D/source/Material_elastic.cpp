#include "../include/Material_elastic.h"
#include "../include/types.h"
#include <cmath> // For std::sqrt

namespace EnSC {

using namespace Types;

void Material_elastic::update() {
    G = E / ((one + v) * two);
    K = E / (one - two * v) / three;
    lambda = two * G * v / (one - two * v);
    WOS = std::sqrt((lambda + two * G) / rho);
}

} // namespace EnSC 