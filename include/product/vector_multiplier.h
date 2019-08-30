#pragma once

#include <vector>

namespace LWS {
    template <typename T>
    class VectorMultiplier {
        public:
        void Multiply(std::vector<double> &v, std::vector<double> &b) const {
            static_cast<T const&>(*this).Multiply(v, b);
        }
    };
}
