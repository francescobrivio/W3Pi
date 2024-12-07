#include "isolation.h"
#include "sum_reduce.h"
#include "best_seed.h"
#include "event_processor.h"
#include "data.h"

/// This function calculates the isolation of the highest non-masked pt particle (the "seed")
/// isolation is defined as the sum of the pts of all particles with a ΔR less then some threshold
void one_iteration(const Puppi in[NPUPPI_MAX], const bool masked[NPUPPI_MAX], bool masked_out[NPUPPI_MAX], Puppi & seed, Puppi::pt_t & absiso) {
    #pragma HLS ARRAY_PARTITION variable=in complete
    #pragma HLS ARRAY_PARTITION variable=masked complete
    #pragma HLS ARRAY_PARTITION variable=masked_out complete
    #pragma HLS pipeline II=9

    #pragma HLS inline off

    const dr2_t dr2_max = float_dr_to_ap_dr2(0.4);
    const dr2_t dr2_veto = float_dr_to_ap_dr2(0.1);

    Puppi::pt_t tosum[NPUPPI_MAX];
    #pragma HLS ARRAY_PARTITION variable=tosum complete

    // 1 find a seed
    seed = find_seed(in, masked);
    // 2 for all particles
    for (unsigned int i = 0; i < NPUPPI_MAX; ++i) {
        // 2.1 calculate ΔR
        dr2_t dr2 = deltaR2(seed, in[i]);
        bool inside = (dr2 < dr2_max);
        // 2.2 mask particles that are inside (including seed)
        masked_out[i] = masked[i] || inside;
        // 2.3 add to sum particles that are `inside` and outside a "veto" radius
        //     is this an important distiction? obviously we don't want to consider the seed
        //     in the isolation calculation, but do we want to exclude other particles that are too close?
        //     this could mean that we mask some particles that are not used in the actual calculation
        //     is this fine?
        tosum[i] =  inside && (dr2 > dr2_veto) ? in[i].hwPt : Puppi::pt_t(0);
    }
    // 3 execute the sum
    absiso = SumReduceAll(tosum);
}

/// This function calculates the isolation of at most `NISO_MAX` particles, chosen as the best seeds
void compute_isolated_l1t(const Puppi in[NPUPPI_MAX], Puppi out[NISO_MAX], Puppi::pt_t out_absiso[NISO_MAX])  {
    #pragma HLS ARRAY_PARTITION variable=in complete
    #pragma HLS ARRAY_PARTITION variable=out complete
    #pragma HLS ARRAY_PARTITION variable=out_absiso complete
    #pragma HLS pipeline II=9

    bool masked[NPUPPI_MAX], masked_out[NPUPPI_MAX];
    #pragma HLS ARRAY_PARTITION variable=masked complete
    #pragma HLS ARRAY_PARTITION variable=masked_out complete
    // 1 mask neutral particles
    for (unsigned int i = 0; i < NPUPPI_MAX; ++i) {
        masked[i] = in[i].hwID <= 1;
    }
    // 2 initialize outputs
    for (unsigned int j = 0; j < NISO_MAX; ++j) {
        out[j].clear();
        out_absiso[j] = 0;
    }
    ap_uint<4> niso = 0;
    // 3 for each output particle
    for (unsigned int j = 0; j < NISO_MAX; ++j) {
        Puppi seed;
        Puppi::pt_t iso_sum;
        // 3.1 calculate isolation of highest pt particle
        one_iteration(in, masked, masked_out, seed, iso_sum);
        // 3.2 copy over mask
        // idea: could modify mask directly in `one_iteration`
        for (unsigned int i = 0; i < NPUPPI_MAX; ++i) {
            masked[i] = masked_out[i];
        }
        // 3.3 save particle and isolation if isolation is less then pt
        //                                 ^ why?
        if (iso_sum < seed.hwPt) {
            out[niso] = seed;
            out_absiso[niso] = iso_sum;
            niso++;
        }
    }
}
