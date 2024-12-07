#ifndef ISOLATION_H
#define ISOLATION_H

#include "data.h"

#define NPUPPI_MAX 216
#define NISO_MAX 12

inline dr2_t float_dr_to_ap_dr2(float dr) {
    return dr2_t(round(pow(dr/Puppi::ETAPHI_LSB,2)));
}
inline float ap_dr2_to_float_dr(dr2_t dr2) {
    return sqrt(dr2.to_int()*Puppi::ETAPHI_LSB*Puppi::ETAPHI_LSB);
}

void compute_isolated_l1t(const Puppi in[NPUPPI_MAX], Puppi out[NISO_MAX], Puppi::pt_t out_absiso[NISO_MAX]) ;

#endif