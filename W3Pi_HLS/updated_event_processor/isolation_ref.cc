#include "src/data.h"
#include "src/event_processor.h"
#include "src/isolation.h"

/// This function calculates the isolation of the given seed
void isolation_of_ref(
    const Puppi & seed,
    const Puppi in[NPUPPI_MAX],
    Puppi::pt_t & abs_iso
) {
    const dr2_t dr2_max = float_dr_to_ap_dr2(0.4);
    const dr2_t dr2_veto = float_dr_to_ap_dr2(0.1);

    abs_iso = Puppi::pt_t(0);
    for (unsigned int i = 0; i < NPUPPI_MAX; ++i) {
        dr2_t dr2 = deltaR2_ref(seed, in[i]);
        if (dr2 < dr2_max && dr2 > dr2_veto) {
            abs_iso += in[i].hwPt;
        }
    }
}

/// Calculates the isolation of each `selected` particle
void calculate_iso_ref(
    const Puppi selected[NPUPPI_SEL],
    const Puppi particles[NPUPPI_MAX],
    Puppi::pt_t out_selected_iso[NPUPPI_SEL]
) {
    for (unsigned int i = 0; i < NPUPPI_SEL; ++i) {
        isolation_of_ref(selected[i], particles, out_selected_iso[i]);
    }
}
