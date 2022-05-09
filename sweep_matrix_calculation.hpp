#include <cstddef>

template <typename T>
void lha_sweep_matrix(T *a, T *b, T *c, T const *Dc, T const *V, std::size_t size, T rdI, T dt) {
    a[1] = static_cast<T>(0);
    b[1] = static_cast<T>(1) + dt * rdI * rdI * (Dc[0] + Dc[1]);
    c[1] = - dt * rdI * (rdI * Dc[1] - T(0.5) * V[1]);

    for (std::size_t elm_idx = 2; elm_idx != size - 2; ++elm_idx) {
        a[elm_idx] = - dt * rdI * (rdI * Dc[elm_idx - 1] + T(0.5) * V[elm_idx]);
        b[elm_idx] = static_cast<T>(1) + dt * rdI * rdI * (Dc[elm_idx] + Dc[elm_idx - 1]);
        c[elm_idx] = - dt * rdI * (rdI * Dc[elm_idx] - T(0.5) * V[elm_idx]);
    } 

    a[size - 2] =  - dt * rdI * (rdI * Dc[size - 3] + T(0.5) * V[size - 2]);
    b[size - 2] = static_cast<T>(1) + dt * rdI * rdI * (Dc[size - 3] + Dc[size - 2]);
    c[size - 2] = static_cast<T>(0);
}

template <typename T>
void prediction_step_rha(
    T *d, T const *f_curr, T const *Dc, T const *V, T rdI,
    T const *F0_curr, T const *F1_curr, 
    std::size_t size, T dt) {
    
    for (std::size_t elm_idx = 1; elm_idx != size - 1; ++elm_idx)
        d[elm_idx] = f_curr[elm_idx] + dt * F0_curr[elm_idx];

    d[1] += f_curr[0] * dt * rdI * (rdI * Dc[0] + static_cast<T>(0.5) * V[1]);
    d[size - 2] += f_curr[size - 1] * dt * rdI * (rdI * Dc[size - 2] - static_cast<T>(0.5) * V[size - 2]);
}

template <typename T>
void correction_step_rha(
    T *d, T const *f_curr, T const *Dc, T const *V, T rdI,
    T const *F0_curr, T const *F1_curr, T const *F0_pred, T const *F1_pred,
    std::size_t size, T dt) {
    
    for (std::size_t elm_idx = 1; elm_idx != size - 1; ++elm_idx)
        d[elm_idx] = f_curr[elm_idx] + 
            static_cast<T>(0.5) * dt * (F0_curr[elm_idx] + F0_pred[elm_idx]) +
            static_cast<T>(0.5) * dt * (F1_curr[elm_idx] - F1_pred[elm_idx]);

    d[1] += f_curr[0] * dt * rdI * (rdI * Dc[0] + static_cast<T>(0.5) * V[1]);
    d[size - 2] += f_curr[size - 1] * dt * rdI * (rdI * Dc[size - 2] - static_cast<T>(0.5) * V[size - 2]);
}