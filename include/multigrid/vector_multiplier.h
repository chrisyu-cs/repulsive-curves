#pragma once

namespace LWS {
    template <typename T>
    class VectorMultiplier {
        public:
        template<typename V, typename Dest>
        void Multiply(V &v, Dest &b) const {
            static_cast<T const&>(*this).Multiply(v, b);
        }
    };
}
