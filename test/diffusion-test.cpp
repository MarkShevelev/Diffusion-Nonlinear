#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <thomas_sweep.hpp>
#include <rha_calculation.hpp>
#include <sweep_matrix_calculation.hpp>
#include <averaging.hpp>

float fPI2 = 6.28318530718f;

template <typename T>
void f_initial(T *f, std::size_t size, std::size_t l) {
    for (std::size_t elm_idx = 1; elm_idx != size / 2; ++elm_idx)
        f[elm_idx] = f[size - elm_idx - 1] = std::sin(fPI2 / l * elm_idx);
    f[0] = f[size - 1] = static_cast<T>(0);
}

template <typename T>
void Dc_initial(T *Dc, std::size_t size) {
    for (std::size_t elm_idx = 0; elm_idx != size; ++elm_idx)
        Dc[elm_idx] = static_cast<T>(1);
    Dc[0] = static_cast<T>(0);
    Dc[size - 2] = static_cast<T>(0);
}

template <typename T>
void V_initial(T *V, std::size_t size, T dI, T kappa) {
    for (std::size_t elm_idx = 1; elm_idx != size / 2; ++elm_idx) {
        auto I = -((size / 2 - elm_idx) * dI - dI / 2);
        V[elm_idx] = - kappa * std::pow(static_cast<T>(1) - I * I, 5.f / 4.f);
        V[size - elm_idx - 1] = V[elm_idx];
    }
    V[0] = V[size - 1] = static_cast<T>(0);
}

template <typename T>
void dV_dI_initial(T *dV_dI, std::size_t size, T dI, T kappa) {
    for (std::size_t elm_idx = 1; elm_idx != size / 2; ++elm_idx) {
        auto I = -((size / 2 - elm_idx) * dI - dI / 2);
        dV_dI[elm_idx] = kappa * 5.f / 4.f * 2.f * I * std::pow(static_cast<T>(1) - I * I, 1.f / 4.f);
        dV_dI[size - elm_idx - 1] = - dV_dI[elm_idx];
    }
    
    dV_dI[0] = dV_dI[size - 1] = static_cast<T>(0);
}

int main() {
    std::size_t const size = 32;
    std::vector<double> f_curr(size), f_pred(size);
    std::vector<double> Dc(size), V(size), dV_dI(size);
    std::vector<double> F0_curr(size), F1_curr(size), F0_pred(size), F1_pred(size);
    std::vector<double> a(size), b(size), c(size), d(size);
    double dI = 2.f / (size - 1), dt = 1.0e-2f, rdI = 1.0f / dI, kappa = 1.0e-2f;
    std::cout << dI << " " << rdI << " " << dt << " " << kappa << " " << rdI * dt << " " << rdI * rdI * dt << std::endl;

    f_initial(f_curr.data(), size, (size - 1) * 2);
    Dc_initial(Dc.data(), size);
    V_initial(V.data(), size, dI, kappa);
    dV_dI_initial(dV_dI.data(), size, dI, kappa);

    std::size_t iter_size = 200;
    std::vector<std::pair<double,double>> avg_vs_time;

    auto f_avg = average_kahan<double>(1, size - 1, [&f_curr](std::size_t idx) {return f_curr[idx];});

    // step
    for (unsigned iter_cnt = 0; iter_cnt != iter_size; ++iter_cnt) {
        /*avg_vs_time.emplace_back(
            dt * iter_cnt,
            std::sqrt(
                average_kahan<double>(1, size - 1, [&f_curr, f_avg](std::size_t idx) {return (f_curr[idx] - f_avg) * (f_curr[idx] - f_avg);})
            )
        );*/

        avg_vs_time.emplace_back(
            dt * iter_cnt,
            average_kahan<double>(1, size - 1, [&f_curr, f_avg] (std::size_t idx) {return f_curr[idx];})
        );

        lha_sweep_matrix(a.data(), b.data(), c.data(), Dc.data(), V.data(), size, rdI, dt);
        F0_calculation(F0_curr.data(), f_curr.data(), dV_dI.data(), size);
        F1_calculation(F1_curr.data(), f_curr.data(), Dc.data(), V.data(), size, rdI);
        prediction_step_rha(d.data(), f_curr.data(), Dc.data(), V.data(), rdI, F0_curr.data(), F1_curr.data(), size, dt);
        thomas_sweep(a.data(), b.data(), c.data(), d.data(), f_pred.data(), size);

        f_pred[0] = f_pred[1];
        f_pred[size - 1] = f_pred[size - 2];

        lha_sweep_matrix(a.data(), b.data(), c.data(), Dc.data(), V.data(), size, rdI, dt);
        F0_calculation(F0_pred.data(), f_pred.data(), dV_dI.data(), size);
        F1_calculation(F1_pred.data(), f_pred.data(), Dc.data(), V.data(), size, rdI);
        correction_step_rha(d.data(), f_curr.data(), Dc.data(), V.data(), rdI, F0_curr.data(), F1_curr.data(), F0_pred.data(), F1_pred.data(), size, dt);
        thomas_sweep(a.data(), b.data(), c.data(), d.data(), f_curr.data(), size);
        
        f_curr[0] = f_curr[1];
        f_curr[size - 1] = f_curr[size - 2];
    }

    {
        std::ofstream fout("./data-two-step.dat");
        fout << std::setprecision(8) << std::scientific;
        for (std::size_t elm_idx = 0; elm_idx != size; ++elm_idx) {
            auto I = elm_idx >= size / 2 ? 
            (elm_idx - size / 2) * dI + dI / 2 : 
            -((size / 2 - elm_idx) * dI - dI / 2);
            fout << I << " " << f_curr[elm_idx] << '\n';
        }
        fout << std::endl;
    }

    {
        std::ofstream fout("./avg-vs-time.dat");
        fout << std::setprecision(8) << std::scientific;
        for (auto a : avg_vs_time) 
            fout << a.first << " " << a.second << '\n';
        fout << std::endl; 
    }

    return 0;
}