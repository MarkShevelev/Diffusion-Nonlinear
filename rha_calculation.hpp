#include <cstddef>

// explicit part of the right hand side
// the size should be even to sustain the rule: dV/dI > 0 F0 = dV_dI * f(-I)
template<typename T>
void F0_calculation(T *F0, T const *f, T const *dV_dI, std::size_t size) {
    for (std::size_t elm_idx = 1; elm_idx != size / 2; ++elm_idx)
        F0[elm_idx] = dV_dI[elm_idx] * f[elm_idx];

    for (std::size_t elm_idx = size / 2; elm_idx != size - 1; ++elm_idx) 
        F0[elm_idx] = -F0[size - elm_idx - 1]; 
}

template<typename T>
void F1_calculation(T *F1, T const *f, T const *Dc, T const *V, std::size_t size, T rdI) {
    for (std::size_t elm_idx = 1; elm_idx != size - 1; ++elm_idx)
        F1[elm_idx] = 
            rdI * rdI * (
                Dc[elm_idx - 1] * f[elm_idx - 1] + 
                Dc[elm_idx] * f[elm_idx + 1] - 
                (Dc[elm_idx - 1] + Dc[elm_idx]) * f[elm_idx]
            )
            - static_cast<T>(0.5) * rdI * (
                V[elm_idx + 1] * f[elm_idx + 1] - V[elm_idx - 1] * f[elm_idx - 1]
            );
}