#include <cstddef>
#include <cmath>

template <typename T>
void thomas_sweep(T *a, T *b, T *c, T *d, T *x_curr, std::size_t size) {
    for (std::size_t elm_idx = 2; elm_idx != size - 1; ++elm_idx) {
			T w = a[elm_idx] / b[elm_idx - 1];
			b[elm_idx] = std::fma(-w, c[elm_idx - 1], b[elm_idx]);
			d[elm_idx] = std::fma(-w, d[elm_idx - 1], d[elm_idx]);
		}
		x_curr[size - 2] = d[size - 2] / b[size - 2];

		for (size_t elm_idx = size - 3; elm_idx != 0; --elm_idx) {
			x_curr[elm_idx] = fma(-c[elm_idx], x_curr[elm_idx + 1], d[elm_idx]) / b[elm_idx];
		}
}