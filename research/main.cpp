#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

#include <sweep_matrix_calculation.hpp>
#include <rha_calculation.hpp>
#include <thomas_sweep.hpp>

template <typename T>
void step(
    T *f_curr, T *f_pred,
    T *a, T *b, T *c, T *d, 
    T *F0_curr, T *F1_curr, T *F0_pred, T *F1_pred, 
    T const *Dc, T const *V, T const *dV_dI, std::size_t size, T rdI, T dt
) {
    lha_sweep_matrix(a, b, c, Dc, V, size, rdI, dt);
    F0_calculation(F0_curr, f_curr, dV_dI, size);
    F1_calculation(F1_curr, f_curr, Dc, V, size, rdI);
    prediction_step_rha(d, f_curr, Dc, V, rdI, F0_curr, F1_curr, size, dt);
    thomas_sweep(a, b, c, d, f_pred, size);

    f_pred[0] = f_pred[1];
    f_pred[size - 1] = f_pred[size - 2];

    lha_sweep_matrix(a, b, c, Dc, V, size, rdI, dt);
    F0_calculation(F0_pred, f_pred, dV_dI, size);
    F1_calculation(F1_pred, f_pred, Dc, V, size, rdI);
    correction_step_rha(d, f_curr, Dc, V, rdI, F0_curr, F1_curr, F0_pred, F1_pred, size, dt);
    thomas_sweep(a, b, c, d, f_curr, size);
    
    f_curr[0] = f_curr[1];
    f_curr[size - 1] = f_curr[size - 2];
}

template <typename T>
struct GetI {
    std::size_t const size;

    GetI(std::size_t size): size(size) { }

    T operator()(std::size_t elm_idx) const {
        auto dI = static_cast<T>(2) / static_cast<T>(size - 1);
        return elm_idx < size / 2 ?
            dI * static_cast<T>(elm_idx) - static_cast<T>(1) :
            - (dI * static_cast<T>(size - elm_idx - 1) - static_cast<T>(1));
    }
};

template <typename T>
void f_initial(T *f, std::size_t size, std::size_t l) {
    auto fPI2 = static_cast<T>(6.28318530718);
    for (std::size_t elm_idx = 1; elm_idx != size / 2; ++elm_idx)
        f[elm_idx] = f[size - elm_idx - 1] = std::sin(fPI2 / l * elm_idx);
    f[0] = f[size - 1] = static_cast<T>(0);
}

template <typename T>
void Dc_initial(T *Dc, std::size_t size, T eps, GetI<T> const &I) {
    for (std::size_t elm_idx = 1; elm_idx != size - 2; ++elm_idx) {
        auto hI = static_cast<T>(0.5) * (I(elm_idx) + I(elm_idx + 1));
        Dc[elm_idx] = eps * eps * std::pow(static_cast<T>(1) - hI * hI, static_cast<T>(2));
    }
    Dc[0] = static_cast<T>(0);
    Dc[size - 2] = static_cast<T>(0);
}

template <typename T>
void V_initial(T *V, std::size_t size, T eps, T kappa, GetI<T> const &I) {
    for (std::size_t elm_idx = 1; elm_idx != size / 2; ++elm_idx) {
        V[elm_idx] = - kappa * std::sqrt(eps) * std::pow(static_cast<T>(1) - I(elm_idx) * I(elm_idx), static_cast<T>(1.25));
        V[size - elm_idx - 1] = V[elm_idx];
    }
    V[0] = V[size - 1] = static_cast<T>(0);
}

template <typename T>
void dV_dI_initial(T *dV_dI, std::size_t size, T eps, T kappa, GetI<T> const &I) {
    for (std::size_t elm_idx = 1; elm_idx != size / 2; ++elm_idx) {
        dV_dI[elm_idx] = kappa * std::sqrt(eps) * static_cast<T>(2.5) * I(elm_idx) * std::pow(static_cast<T>(1) - I(elm_idx) * I(elm_idx), static_cast<T>(0.25));
        dV_dI[size - elm_idx - 1] = - dV_dI[elm_idx];
    }
    
    dV_dI[0] = dV_dI[size - 1] = static_cast<T>(0);
}

int main() {
    std::size_t const size = 32;
    std::vector<double> 
        f_curr(size), f_pred(size),
        Dc(size), V(size), dV_dI(size),
        F0_curr(size), F1_curr(size), F0_pred(size), F1_pred(size),
        a(size), b(size), c(size), d(size);
    double dI = 2.f / (size - 1), dt = 1.0e-2f, rdI = 1.0f / dI, eps = 1.f, kappa = 0.0;
    GetI<double> I(size);
    std::cout << dI << " " << rdI << " " << dt << " " << kappa << " " << rdI * dt << " " << rdI * rdI * dt << std::endl;

    f_initial(f_curr.data(), size, (size - 1) * 2);
    Dc_initial(Dc.data(), size, eps, I);
    V_initial(V.data(), size, eps, kappa, I);
    dV_dI_initial(dV_dI.data(), size, eps, kappa, I);

    for (unsigned iter_cnt = 0; iter_cnt != 300; ++iter_cnt) {
        step(
            f_curr.data(), f_pred.data(),
            a.data(), b.data(), c.data(), d.data(),
            F0_curr.data(), F1_curr.data(), F0_pred.data(), F1_pred.data(),
            Dc.data(), V.data(), dV_dI.data(), size, rdI, dt
        );
    }

    // export
    {
        std::ofstream fout("./data.txt");
        if (!fout)
            std::cerr << "Can't open file" << std::endl;
        else {
            fout << std::setprecision(8) << std::scientific;
            for (std::size_t elm_idx = 0; elm_idx != size; ++elm_idx) {
                fout << I(elm_idx) << " " << f_curr[elm_idx] << '\n';
            }
            fout << std::endl;
        }
    }

    return 0;
}