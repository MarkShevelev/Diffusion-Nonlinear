#include <iostream>

#include <thomas_sweep.hpp>

int main() {
    std::size_t const size = 128;
    float a[size], b[size], c[size], d[size], x[size];
    x[0] = x[size - 1] = 1.f;

    for (std::size_t elm_idx = 1; elm_idx != size - 1; ++elm_idx) {
        a[elm_idx] = 1.0f;
        b[elm_idx] = 3.0f;
        c[elm_idx] = 1.0f;
        d[elm_idx] = 5.0f;
    }
    d[1] = d[size - 2] = 4.0f;
    a[1] = 0.f;
    c[size - 2] = 0.f;

    thomas_sweep(a, b, c, d, x, size);

    for (std::size_t elm_idx = 1; elm_idx != size - 1; ++elm_idx)
        std::cout << x[elm_idx] << '\n';
    std::cout << std::endl;

    return 0;
}