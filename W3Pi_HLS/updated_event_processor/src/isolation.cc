#include "isolation.h"
#include "sum_reduce.h"
#include "event_processor.h"
#include "data.h"

/// This function calculates the isolation of the given seed
void isolation_of(
    const Puppi &seed,
    const Puppi particles[NPUPPI_MAX],
    Puppi::pt_t &abs_iso
) {
    #pragma HLS ARRAY_PARTITION variable=particles complete
    #pragma HLS pipeline II=9

    const dr2_t dr2_max = float_dr_to_ap_dr2(0.4);
    const dr2_t dr2_veto = float_dr_to_ap_dr2(0.1);

    Puppi::pt_t tosum[NPUPPI_MAX];
    #pragma HLS ARRAY_PARTITION variable=tosum complete

    // 1 for all particles
    for (unsigned int i = 0; i < NPUPPI_MAX; ++i) {
        // 1.1 calculate ΔR
        dr2_t dr2 = deltaR2(seed, particles[i]);
        bool inside = (dr2 < dr2_max);
        // 1.2 add to sum particles that are `inside` and outside a "veto" radius
        tosum[i] =  inside && (dr2 > dr2_veto) ? particles[i].hwPt : Puppi::pt_t(0);
    }
    // 2 execute the sum
    abs_iso = SumReduceAll(tosum);
}

/// Calculates the isolation of each `selected` particle
void calculate_iso(
    const Puppi selected[NPUPPI_SEL],
    const Puppi particles[NPUPPI_MAX],
    Puppi::pt_t out_selected_iso[NPUPPI_SEL]
) {
    #pragma HLS ARRAY_PARTITION variable=selected complete
    #pragma HLS ARRAY_PARTITION variable=particles complete
    #pragma HLS ARRAY_PARTITION variable=out_selected_iso complete
    for (unsigned int i = 0; i < NPUPPI_SEL; ++i) {
        // the greater the factor, the better the latency, but it uses more resources
        // |factor|latency|lut%util|
        // |     1|     73|       2|
        // |     2|     55|       4|
        // |     3|     46|       7|
        // |     7|     16|      17|
        #pragma HLS unroll factor=3
        isolation_of(selected[i], particles, out_selected_iso[i]);
    }
}
