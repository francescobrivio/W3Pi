#ifndef ISOLATION_H
#define ISOLATION_H

#include "data.h"

inline dr2_t float_dr_to_ap_dr2(float dr) {
    return dr2_t(round(pow(dr/Puppi::ETAPHI_LSB,2)));
}
inline float ap_dr2_to_float_dr(dr2_t dr2) {
    return sqrt(dr2.to_int()*Puppi::ETAPHI_LSB*Puppi::ETAPHI_LSB);
}

void compute_isolated_l1t(const Puppi in[NPUPPI_MAX], Puppi out[NISO_MAX], Puppi::pt_t out_absiso[NISO_MAX]) ;
void compute_isolated_l1t_parametrized_cut(const Puppi in[NPUPPI_MAX], Puppi out[NISO_MAX], Puppi::pt_t out_absiso[NISO_MAX]) ;
void compute_isolated_l1t_no_cut(const Puppi in[NPUPPI_MAX], Puppi out[NISO_MAX], Puppi::pt_t out_absiso[NISO_MAX]) ;

void calculate_iso(const Puppi selected[NPUPPI_SEL], const Puppi particles[NPUPPI_MAX], Puppi::pt_t out_slected_iso[NPUPPI_SEL]);

#endif