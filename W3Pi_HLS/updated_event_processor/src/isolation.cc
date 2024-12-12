#include "isolation.h"
#include "sum_reduce.h"
#include "best_seed.h"
#include "event_processor.h"
#include "data.h"

/// This function calculates the isolation of the highest non-masked pt particle (the "seed")
/// isolation is defined as the sum of the pts of all particles with a ΔR less then some threshold
void one_iteration(
    const Puppi in[NPUPPI_MAX],
    bool masked[NPUPPI_MAX],
    Puppi & seed,
    Puppi::pt_t & abs_iso
) {
    #pragma HLS ARRAY_PARTITION variable=in complete
    #pragma HLS ARRAY_PARTITION variable=masked complete
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
        masked[i] |= inside;
        // 2.3 add to sum particles that are `inside` and outside a "veto" radius
        tosum[i] =  inside && (dr2 > dr2_veto) ? in[i].hwPt : Puppi::pt_t(0);
    }
    // 3 execute the sum
    abs_iso = SumReduceAll(tosum);
}

/// like `one_iteration`, but we mask only the seed, and nothing else
void one_iteration_no_mask(
    const Puppi in[NPUPPI_MAX],
    bool masked[NPUPPI_MAX],
    Puppi & seed,
    Puppi::pt_t & abs_iso
) {
    #pragma HLS ARRAY_PARTITION variable=in complete
    #pragma HLS ARRAY_PARTITION variable=masked complete
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
        // 2.2 mask the seed
        masked[i] |= in[i].pack() == seed.pack();
        // 2.3 add to sum particles that are `inside` and outside a "veto" radius
        tosum[i] =  inside && (dr2 > dr2_veto) ? in[i].hwPt : Puppi::pt_t(0);
    }
    // 3 execute the sum
    abs_iso = SumReduceAll(tosum);
}

/// This function calculates the isolation of the given seed
void isolation_of(
    const Puppi & seed,
    const Puppi in[NPUPPI_MAX],
    Puppi::pt_t & abs_iso
) {
    #pragma HLS ARRAY_PARTITION variable=in complete
    #pragma HLS pipeline II=9

    const dr2_t dr2_max = float_dr_to_ap_dr2(0.4);
    const dr2_t dr2_veto = float_dr_to_ap_dr2(0.1);

    Puppi::pt_t tosum[NPUPPI_MAX];
    #pragma HLS ARRAY_PARTITION variable=tosum complete

    // 1 for all particles
    for (unsigned int i = 0; i < NPUPPI_MAX; ++i) {
        // 1.1 calculate ΔR
        dr2_t dr2 = deltaR2(seed, in[i]);
        bool inside = (dr2 < dr2_max);
        // 1.2 add to sum particles that are `inside` and outside a "veto" radius
        tosum[i] =  inside && (dr2 > dr2_veto) ? in[i].hwPt : Puppi::pt_t(0);
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
        #pragma HLS unroll factor=NPUPPI_SEL
        isolation_of(selected[i], particles, out_selected_iso[i]);
    }
}

/// This function calculates the isolation of at most `NISO_MAX` particles, chosen as the best seeds
void compute_isolated_l1t(
    const Puppi in[NPUPPI_MAX],
    Puppi out[NISO_MAX],
    Puppi::pt_t out_absiso[NISO_MAX]
) {
    #pragma HLS ARRAY_PARTITION variable=in complete
    #pragma HLS ARRAY_PARTITION variable=out complete
    #pragma HLS ARRAY_PARTITION variable=out_absiso complete
    #pragma HLS pipeline II=9

    bool masked[NPUPPI_MAX];
    #pragma HLS ARRAY_PARTITION variable=masked complete
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
        Puppi::pt_t abs_iso;
        // 3.1 calculate isolation of highest pt particle
        one_iteration(in, masked, seed, abs_iso);
        // 3.3 save particle and isolation if isolation is less then pt
        if (abs_iso < seed.hwPt) {
            out[niso] = seed;
            out_absiso[niso] = abs_iso;
            niso++;
        }
    }
}

/// like `compute_isolated_l1t`, but the cut is parametrized
void compute_isolated_l1t_parametrized_cut(
    const Puppi in[NPUPPI_MAX],
    Puppi out[NISO_MAX],
    Puppi::pt_t out_absiso[NISO_MAX]
) {
    #pragma HLS ARRAY_PARTITION variable=in complete
    #pragma HLS ARRAY_PARTITION variable=out complete
    #pragma HLS ARRAY_PARTITION variable=out_absiso complete
    #pragma HLS pipeline II=9

    bool masked[NPUPPI_MAX];
    #pragma HLS ARRAY_PARTITION variable=masked complete
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
        Puppi::pt_t abs_iso;
        // 3.1 calculate isolation of highest pt particle
        one_iteration_no_mask(in, masked, seed, abs_iso);
        // 3.3 save particle and isolation if isolation is less then pt * ISO_CUT_CONST
        if (abs_iso < ISO_CUT_CONST * seed.hwPt) {
            out[niso] = seed;
            out_absiso[niso] = abs_iso;
            niso++;
        }
    }
}

/// like `compute_isolated_l1t`, but there is no cut
void compute_isolated_l1t_no_cut(
    const Puppi in[NPUPPI_MAX],
    Puppi out[NISO_MAX],
    Puppi::pt_t out_absiso[NISO_MAX]
) {
    #pragma HLS ARRAY_PARTITION variable=in complete
    #pragma HLS ARRAY_PARTITION variable=out complete
    #pragma HLS ARRAY_PARTITION variable=out_absiso complete
    #pragma HLS pipeline II=9

    bool masked[NPUPPI_MAX];
    #pragma HLS ARRAY_PARTITION variable=masked complete
    // 1 mask neutral particles
    for (unsigned int i = 0; i < NPUPPI_MAX; ++i) {
        masked[i] = in[i].hwID <= 1;
    }
    // 2 for each output particle
    for (unsigned int j = 0; j < NISO_MAX; ++j) {
        Puppi seed;
        Puppi::pt_t abs_iso;
        // 2.1 calculate isolation of highest pt particle
        one_iteration_no_mask(in, masked, seed, abs_iso);
        // 2.2 save particle and isolation
        out[j] = seed;
        out_absiso[j] = abs_iso;
    }
}
