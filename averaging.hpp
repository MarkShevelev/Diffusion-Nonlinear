#include <cstddef>

template <typename T, typename F>
auto sum_kahan(std::size_t begin, std::size_t end, F f) {
    T r = static_cast<T>(0), accum = static_cast<T>(0);
    T y, t;
    for (std::size_t it = begin; it != end; ++it) {
        y = f(it) - c;
        t = y + accum;
        c = (t - accum) - y;
        accum = t;
    }
    return accum;
}

template <typename T, typename F>
auto average_kahan(std::size_t begin, std::size_t end, F f) {
    return sum_kahan(begin, end, f) / (end - begin - 1);
}