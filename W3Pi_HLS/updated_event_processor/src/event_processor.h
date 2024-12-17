#ifndef EVENT_PROCESSOR_H
#define EVENT_PROCESSOR_H

#include "data.h"
#include "../BDT/w3p_bdt.h"

// --------------------
// ----- FIRMWARE -----
// --------------------
void masker        (const Puppi input[NPUPPI_MAX], ap_uint<NPUPPI_MAX> & masked);
void slimmer       (const Puppi input[NPUPPI_MAX], const ap_uint<NPUPPI_MAX> masked, Puppi slimmed[NPUPPI_MAX]);
void orderer7      (const Puppi slimmed[NPUPPI_MAX],
                    Puppi ordered1 [NSPLITS], Puppi ordered2 [NSPLITS],
                    Puppi ordered3 [NSPLITS], Puppi ordered4 [NSPLITS],
                    Puppi ordered5 [NSPLITS], Puppi ordered6 [NSPLITS],
                    Puppi ordered7 [NSPLITS], Puppi ordered8 [NSPLITS],
                    Puppi ordered9 [NSPLITS], Puppi ordered10[NSPLITS],
                    Puppi ordered11[NSPLITS], Puppi ordered12[NSPLITS],
                    Puppi ordered13[NSPLITS], Puppi ordered14[NSPLITS],
                    Puppi ordered15[NSPLITS], Puppi ordered16[NSPLITS]);
void orderer7bis   (const Puppi slimmed[NPUPPI_MAX],
                    Puppi ordered1 [NSPLITS], Puppi ordered2 [NSPLITS],
                    Puppi ordered3 [NSPLITS], Puppi ordered4 [NSPLITS],
                    Puppi ordered5 [NSPLITS], Puppi ordered6 [NSPLITS],
                    Puppi ordered7 [NSPLITS], Puppi ordered8 [NSPLITS]);
void orderer7f     (const Puppi slimmed[NPUPPI_MAX], Puppi ordered[NSUBARR][NSPLITS]);
void merger        (Puppi ordered1 [NSPLITS], Puppi ordered2 [NSPLITS],
                    Puppi ordered3 [NSPLITS], Puppi ordered4 [NSPLITS],
                    Puppi ordered5 [NSPLITS], Puppi ordered6 [NSPLITS],
                    Puppi ordered7 [NSPLITS], Puppi ordered8 [NSPLITS],
                    Puppi ordered9 [NSPLITS], Puppi ordered10[NSPLITS],
                    Puppi ordered11[NSPLITS], Puppi ordered12[NSPLITS],
                    Puppi ordered13[NSPLITS], Puppi ordered14[NSPLITS],
                    Puppi ordered15[NSPLITS], Puppi ordered16[NSPLITS],
                    Puppi merged[NPUPPI_MAX]);
void merger7bis    (Puppi ordered1 [NSPLITS], Puppi ordered2 [NSPLITS],
                    Puppi ordered3 [NSPLITS], Puppi ordered4 [NSPLITS],
                    Puppi ordered5 [NSPLITS], Puppi ordered6 [NSPLITS],
                    Puppi ordered7 [NSPLITS], Puppi ordered8 [NSPLITS],
                    Puppi merged[NPUPPI_MAX]);
void merger7f      (Puppi ordered[NSUBARR][NSPLITS], Puppi merged[NPUPPI_MAX]);
void selector      (const Puppi merged[NPUPPI_MAX], Puppi selected[NPUPPI_SEL]);
void get_triplet_inputs(const Puppi selected[NPUPPI_SEL], idx_t idx0, idx_t idx1, idx_t idx2, w3p_bdt::input_t BDT_inputs[w3p_bdt::n_features]);
cos_t get_cos_phi      (Puppi::phi_t phi);
cosh_t get_cosh_eta    (Puppi::eta_t eta);
mass_t get_pair_mass   (const Puppi & p1, const Puppi & p2);
void get_event_inputs  (const Puppi selected[NPUPPI_SEL], w3p_bdt::input_t BDT_inputs[NTRIPLETS][w3p_bdt::n_features]);
void get_event_scores  (w3p_bdt::input_t BDT_inputs[NTRIPLETS][w3p_bdt::n_features], w3p_bdt::score_t BDT_scores[NTRIPLETS]);
void get_highest_score (w3p_bdt::score_t BDT_scores[NTRIPLETS], w3p_bdt::score_t & high_score);
void EventProcessor    (const Puppi input[NPUPPI_MAX], w3p_bdt::score_t & max_score);
void EventProcessor7bis(const Puppi input[NPUPPI_MAX], w3p_bdt::score_t & max_score);
void EventProcessor7f  (const Puppi input[NPUPPI_MAX], w3p_bdt::score_t & max_score);

// ---------------------
// ----- REFERENCE -----
// ---------------------
void masker_ref        (const Puppi input[NPUPPI_MAX], ap_uint<NPUPPI_MAX> & masked);
void slimmer_ref       (const Puppi input[NPUPPI_MAX], const ap_uint<NPUPPI_MAX> masked, Puppi slimmed[NPUPPI_MAX]);
void orderer7_ref      (const Puppi slimmed[NPUPPI_MAX],
                        Puppi ordered1 [NSPLITS], Puppi ordered2 [NSPLITS],
                        Puppi ordered3 [NSPLITS], Puppi ordered4 [NSPLITS],
                        Puppi ordered5 [NSPLITS], Puppi ordered6 [NSPLITS],
                        Puppi ordered7 [NSPLITS], Puppi ordered8 [NSPLITS],
                        Puppi ordered9 [NSPLITS], Puppi ordered10[NSPLITS],
                        Puppi ordered11[NSPLITS], Puppi ordered12[NSPLITS],
                        Puppi ordered13[NSPLITS], Puppi ordered14[NSPLITS],
                        Puppi ordered15[NSPLITS], Puppi ordered16[NSPLITS]);
void orderer7bis_ref   (const Puppi slimmed[NPUPPI_MAX],
                        Puppi ordered1 [NSPLITS], Puppi ordered2 [NSPLITS],
                        Puppi ordered3 [NSPLITS], Puppi ordered4 [NSPLITS],
                        Puppi ordered5 [NSPLITS], Puppi ordered6 [NSPLITS],
                        Puppi ordered7 [NSPLITS], Puppi ordered8 [NSPLITS]);
void merger_ref        (Puppi slimmed[NPUPPI_MAX], Puppi merged[NPUPPI_MAX]);
void selector_ref      (const Puppi merged[NPUPPI_MAX], Puppi selected[NPUPPI_SEL]);
mass_t get_pair_mass_ref   (const Puppi & p1, const Puppi & p2);
void get_triplet_inputs_ref(const Puppi selected[NPUPPI_SEL], idx_t idx0, idx_t idx1, idx_t idx2, w3p_bdt::input_t BDT_inputs[w3p_bdt::n_features]);
void get_event_inputs_ref  (const Puppi selected[NPUPPI_SEL], w3p_bdt::input_t BDT_inputs[NTRIPLETS][w3p_bdt::n_features]);
void get_event_scores_ref  (w3p_bdt::input_t BDT_inputs[NTRIPLETS][w3p_bdt::n_features], w3p_bdt::score_t BDT_scores[NTRIPLETS]);
void get_highest_score_ref (w3p_bdt::score_t BDT_scores[NTRIPLETS], w3p_bdt::score_t & high_score);
void EventProcessor_ref(const Puppi input[NPUPPI_MAX], w3p_bdt::score_t & max_score);

#endif
