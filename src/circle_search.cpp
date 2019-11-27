#include "circle_search.h"
#include "flow/gradient_constraint_types.h"

namespace LWS {

    CircleSearch::CircleSearch(PolyCurveNetwork* c, double a, double b, double eps) {
        curves = c;
        alpha = a;
        beta = b;
        epsilon = eps;
    }
}