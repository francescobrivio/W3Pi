#include "src/event_processor.h"
#include <cstdio>
#include <cstdlib>
#include <algorithm>

#include "BDT/conifer.h"

// ------------------------------------------------------------------
// Masker: mask L1Puppi objects that don't pass selections
void masker_ref (const Puppi input[NPUPPI_MAX], ap_uint<NPUPPI_MAX> & masked)
{
    masked = 0;

    for (unsigned int i = 0; i < NPUPPI_MAX; i++)
    {
        masked[i] = ( 
            //( input[i].hwPt < 3.25 )                      ||
            ( std::abs(input[i].hwEta) > Puppi::ETA_CUT ) ||
            ( input[i].hwID < 2 || input[i].hwID > 5 )
        );
    }
}

// ------------------------------------------------------------------
// Slimmer: replace masked candidates with empty/dummy puppi
void slimmer_ref (const Puppi input[NPUPPI_MAX], const ap_uint<NPUPPI_MAX> masked, Puppi slimmed[NPUPPI_MAX])
{
    for (int i = 0; i < NPUPPI_MAX; i++)
        slimmed[i].clear();

    for (int i = 0; i < NPUPPI_MAX; i++)
    {
        if (!masked[i])
        {
            slimmed[i] = input[i];
        }
    }
}

// ---------------------------------
// Sort slimmed candidates by pT using bitonicSort from bitonic_hybrid.h

bool puppiComparator (Puppi a, Puppi b)
{
    return (a.hwPt > b.hwPt);
}

// Split and order 16 arrays of 13 candidates (16 x 13 = 208)
void orderer7_ref (const Puppi slimmed[NPUPPI_MAX],
                   Puppi ordered1 [NSPLITS], Puppi ordered2 [NSPLITS],
                   Puppi ordered3 [NSPLITS], Puppi ordered4 [NSPLITS],
                   Puppi ordered5 [NSPLITS], Puppi ordered6 [NSPLITS],
                   Puppi ordered7 [NSPLITS], Puppi ordered8 [NSPLITS],
                   Puppi ordered9 [NSPLITS], Puppi ordered10[NSPLITS],
                   Puppi ordered11[NSPLITS], Puppi ordered12[NSPLITS],
                   Puppi ordered13[NSPLITS], Puppi ordered14[NSPLITS],
                   Puppi ordered15[NSPLITS], Puppi ordered16[NSPLITS])
{
    // Split in 16 arrays
    for (unsigned int i=0; i < NSPLITS; i++)
    {
        #pragma HLS UNROLL
        ordered1 [i] = slimmed[i];
        ordered2 [i] = slimmed[i+NSPLITS];
        ordered3 [i] = slimmed[i+2*NSPLITS];
        ordered4 [i] = slimmed[i+3*NSPLITS];
        ordered5 [i] = slimmed[i+4*NSPLITS];
        ordered6 [i] = slimmed[i+5*NSPLITS];
        ordered7 [i] = slimmed[i+6*NSPLITS];
        ordered8 [i] = slimmed[i+7*NSPLITS];
        ordered9 [i] = slimmed[i+8*NSPLITS];
        ordered10[i] = slimmed[i+9*NSPLITS];
        ordered11[i] = slimmed[i+10*NSPLITS];
        ordered12[i] = slimmed[i+11*NSPLITS];
        ordered13[i] = slimmed[i+12*NSPLITS];
        ordered14[i] = slimmed[i+13*NSPLITS];
        ordered15[i] = slimmed[i+14*NSPLITS];
        ordered16[i] = slimmed[i+15*NSPLITS];
    }

    // Sort arrays
    std::stable_sort(ordered1 , ordered1  + NSPLITS, puppiComparator);
    std::stable_sort(ordered2 , ordered2  + NSPLITS, puppiComparator);
    std::stable_sort(ordered3 , ordered3  + NSPLITS, puppiComparator);
    std::stable_sort(ordered4 , ordered4  + NSPLITS, puppiComparator);
    std::stable_sort(ordered5 , ordered5  + NSPLITS, puppiComparator);
    std::stable_sort(ordered6 , ordered6  + NSPLITS, puppiComparator);
    std::stable_sort(ordered7 , ordered7  + NSPLITS, puppiComparator);
    std::stable_sort(ordered8 , ordered8  + NSPLITS, puppiComparator);
    std::stable_sort(ordered9 , ordered9  + NSPLITS, puppiComparator);
    std::stable_sort(ordered10, ordered10 + NSPLITS, puppiComparator);
    std::stable_sort(ordered11, ordered11 + NSPLITS, puppiComparator);
    std::stable_sort(ordered12, ordered12 + NSPLITS, puppiComparator);
    std::stable_sort(ordered13, ordered13 + NSPLITS, puppiComparator);
    std::stable_sort(ordered14, ordered14 + NSPLITS, puppiComparator);
    std::stable_sort(ordered15, ordered15 + NSPLITS, puppiComparator);
    std::stable_sort(ordered16, ordered16 + NSPLITS, puppiComparator);
}

// Split and order 8 arrays of 26 candidates (8 x 26 = 208)
void orderer7bis_ref (const Puppi slimmed[NPUPPI_MAX],
                   Puppi ordered1 [NSPLITS], Puppi ordered2 [NSPLITS],
                   Puppi ordered3 [NSPLITS], Puppi ordered4 [NSPLITS],
                   Puppi ordered5 [NSPLITS], Puppi ordered6 [NSPLITS],
                   Puppi ordered7 [NSPLITS], Puppi ordered8 [NSPLITS])
{
    // Split in 8 arrays
    for (unsigned int i=0; i < NSPLITS; i++)
    {
        #pragma HLS UNROLL
        ordered1[i] = slimmed[i];
        ordered2[i] = slimmed[i+NSPLITS];
        ordered3[i] = slimmed[i+2*NSPLITS];
        ordered4[i] = slimmed[i+3*NSPLITS];
        ordered5[i] = slimmed[i+4*NSPLITS];
        ordered6[i] = slimmed[i+5*NSPLITS];
        ordered7[i] = slimmed[i+6*NSPLITS];
        ordered8[i] = slimmed[i+7*NSPLITS];
    }

    // Sort arrays
    std::stable_sort(ordered1, ordered1 + NSPLITS, puppiComparator);
    std::stable_sort(ordered2, ordered2 + NSPLITS, puppiComparator);
    std::stable_sort(ordered3, ordered3 + NSPLITS, puppiComparator);
    std::stable_sort(ordered4, ordered4 + NSPLITS, puppiComparator);
    std::stable_sort(ordered5, ordered5 + NSPLITS, puppiComparator);
    std::stable_sort(ordered6, ordered6 + NSPLITS, puppiComparator);
    std::stable_sort(ordered7, ordered7 + NSPLITS, puppiComparator);
    std::stable_sort(ordered8, ordered8 + NSPLITS, puppiComparator);
}

// ------------------------------------------------------------------
// Merger
void merger_ref (Puppi slimmed[NPUPPI_MAX], Puppi merged[NPUPPI_MAX])
{
    // Clean initial values
    for (int i = 0; i < NPUPPI_MAX; i++)
        merged[i].clear();

    // Sort input candidates
    std::stable_sort(slimmed, slimmed + NPUPPI_MAX, puppiComparator);

    // Select highest pT candidates
    for (int i = 0; i < NPUPPI_MAX; i++)
        merged[i] = slimmed[i];
}

// ------------------------------------------------------------------
// Select first NPUPPI_SEL from the ordered list
void selector_ref(const Puppi merged[NPUPPI_MAX], Puppi selected[NPUPPI_SEL])
{
    for (unsigned int i = 0; i < NPUPPI_SEL; i++)
        selected[i] = merged[i];
}

// ------------------------------------------------------------------
// Get maximum deltaVz
Puppi::z0_t get_max_dVz_ref (Puppi::z0_t z0, Puppi::z0_t z1, Puppi::z0_t z2)
{
    Puppi::z0_t dVz_01 = z0 - z1;
    Puppi::z0_t dVz_02 = z0 - z2;
    Puppi::z0_t dVz_12 = z1 - z2;

    Puppi::z0_t max = std::max(dVz_01, std::max(dVz_02, dVz_12));

    return max;
}

// ------------------------------------------------------------------
// DeltaR methods
inline dr2_t deltaR2_ref(const Puppi & p1, const Puppi & p2)
{
    auto dphi = p1.hwPhi - p2.hwPhi;
    if (dphi > Puppi::INT_PI) dphi -= Puppi::INT_2PI;
    else if (dphi < -Puppi::INT_PI) dphi += Puppi::INT_2PI;
    auto deta = p1.hwEta - p2.hwEta;
    return dphi*dphi + deta*deta;
}

dr2_t get_min_deltaR2_ref (const Puppi p0, const Puppi p1, const Puppi p2)
{
    // Get the 3 dR values
    dr2_t dr2_01 = deltaR2_ref(p0, p1);
    dr2_t dr2_02 = deltaR2_ref(p0, p2);
    dr2_t dr2_12 = deltaR2_ref(p1, p2);

    // std library min
    dr2_t min = std::min(dr2_01, std::min(dr2_02, dr2_12));

    return min;
}

// ------------------------------------------------------------------
// Invariant mass quared of pair of particles
// m^2 = 2 * pT1 * pT2 * ( cosh(eta1 - eta2) - cos(phi1 - phi2) )
mass_t get_pair_mass_ref (const Puppi & p1, const Puppi & p2)
{
    // Get dPhi and dEta
    auto dphi = p1.hwPhi - p2.hwPhi;
    if (dphi > Puppi::INT_PI) dphi -= Puppi::INT_2PI;
    else if (dphi < -Puppi::INT_PI) dphi += Puppi::INT_2PI;
    auto deta = p1.hwEta - p2.hwEta;

    // Get cos and cosh
    cos_t cosdPhi   = std::cos(Puppi::floatPhi(dphi));
    cosh_t coshdEta = std::cosh(Puppi::floatEta(deta));

    // Compute and return final mass
    mass_t pair_m = 2 * p1.hwPt * p2.hwPt * ( coshdEta - cosdPhi );

    return pair_m;
}

// ------------------------------------------------------------------
// Prepare inputs features:
//  0. pi2_pt
//  1. pi1_pt
//  2. m_01
//  3. triplet_charge
//  4. dVz_02
//  5. pi0_pt
//  6. triplet_pt
//  7. m_02
//  8. triplet_maxdVz
//  9. triplet_mindR
// 10. pi2_eta
void get_triplet_inputs_ref (const Puppi selected[NPUPPI_SEL], idx_t idx0, idx_t idx1, idx_t idx2, w3p_bdt::input_t BDT_inputs[w3p_bdt::n_features])
{
    // Fill BDT input for triplet
    BDT_inputs[0]  = selected[idx2].hwPt;
    BDT_inputs[1]  = selected[idx1].hwPt;
    BDT_inputs[2]  = get_pair_mass_ref(selected[idx0], selected[idx1]);
    BDT_inputs[4]  = selected[idx0].hwZ0 - selected[idx2].hwZ0;
    BDT_inputs[5]  = selected[idx0].hwPt;
    BDT_inputs[6]  = selected[idx0].hwPt + selected[idx1].hwPt + selected[idx2].hwPt;
    BDT_inputs[7]  = get_pair_mass_ref(selected[idx0], selected[idx2]);
    BDT_inputs[8]  = get_max_dVz_ref(selected[idx0].hwZ0, selected[idx1].hwZ0, selected[idx2].hwZ0);
    BDT_inputs[9]  = get_min_deltaR2_ref(selected[idx0], selected[idx1], selected[idx2]);
    BDT_inputs[10] = selected[idx2].hwEta;
}

// ------------------------------------------------------------------
// Get all event inputs
// - 8 triplets from: 5 from 1st+2nd, 2 from 1st+3rd and 1 from 2nd+3rd:
// (0,1,2)-(0,1,3)-(0,1,4)-(0,1,5)-(0,1,6)-(0,2,3)-(0,2,3)-(1,2,3)
void get_event_inputs_ref (const Puppi selected[NPUPPI_SEL], w3p_bdt::input_t BDT_inputs[NTRIPLETS][w3p_bdt::n_features])
{
    get_triplet_inputs_ref(selected, 0, 1, 2, BDT_inputs[0]);
    get_triplet_inputs_ref(selected, 0, 1, 3, BDT_inputs[1]);
    get_triplet_inputs_ref(selected, 0, 1, 4, BDT_inputs[2]);
    get_triplet_inputs_ref(selected, 0, 1, 5, BDT_inputs[3]);
    get_triplet_inputs_ref(selected, 0, 1, 6, BDT_inputs[4]);
    get_triplet_inputs_ref(selected, 0, 2, 3, BDT_inputs[5]);
    get_triplet_inputs_ref(selected, 0, 2, 4, BDT_inputs[6]);
    get_triplet_inputs_ref(selected, 1, 2, 3, BDT_inputs[7]);
}

// ------------------------------------------------------------------
// Transform array into vector
template<typename T>
std::vector<T> inputs_to_vec(T A[], int size)
{
    std::vector<T> vec_result;
    for (unsigned int i = 0; i < size; i++)
    {
        vec_result.push_back(A[i]);
    }
    return vec_result;
}

// ------------------------------------------------------------------
// Get BDT scores of the selected triplets
// from: https://github.com/cms-sw/cmssw/blob/master/L1Trigger/Phase2L1ParticleFlow/src/egamma/pftkegalgo_ref.cpp#L323-L325
void get_event_scores_ref (w3p_bdt::input_t BDT_inputs[NTRIPLETS][w3p_bdt::n_features], w3p_bdt::score_t BDT_scores[NTRIPLETS])
{
    // Declare BDT
    conifer::BDT<w3p_bdt::input_t, w3p_bdt::score_t, false> *composite_bdt_;

    // Load model
    std::string resolvedFileName = "conifer_binary_featV4_finalFit_v5.json";
    composite_bdt_ = new conifer::BDT<w3p_bdt::input_t, w3p_bdt::score_t, false>(resolvedFileName);

    // Get score for each triplet
    for (unsigned int i = 0; i < NTRIPLETS; i++)
    {
        // Make vec of input from array of inputs
        std::vector<w3p_bdt::input_t> triplet_vec_inputs = inputs_to_vec<w3p_bdt::input_t>(BDT_inputs[i], w3p_bdt::n_features);

        // Compure score
        std::vector<w3p_bdt::score_t> bdt_score = composite_bdt_->decision_function(triplet_vec_inputs);

        // Store score in output variable
        BDT_scores[i] = bdt_score.at(0);
    }
}

// ------------------------------------------------------------------
// Get highest BDT score
void get_highest_score_ref (w3p_bdt::score_t BDT_scores[NTRIPLETS], w3p_bdt::score_t & high_score)
{
    // Sort the BDT scores in descending order
    std::stable_sort(BDT_scores , BDT_scores  + NTRIPLETS, std::greater<w3p_bdt::score_t>());

    // Return highest score
    high_score = BDT_scores[0];
}

// ------------------------------------------------------------------
// Full EventProcessor
void EventProcessor_ref (const Puppi input[NPUPPI_MAX], w3p_bdt::score_t & max_score)
{
    // Mask candidates
    ap_uint<NPUPPI_MAX> masked;
    masker_ref(input, masked);

    // Replace masked candidates with dummy
    Puppi slimmed[NPUPPI_MAX];
    slimmer_ref(input, masked, slimmed);

    // No need for splitting and ordering in reference C++ code: it can be done in one go

    // "Merge" and sort the slimmed candidates
    Puppi merged[NPUPPI_MAX];
    merger_ref(slimmed, merged);

    // Select only highest pT ordered-candidates
    Puppi selected[NPUPPI_SEL];
    selector_ref(merged, selected);

    // Get inputs for each triplet
    w3p_bdt::input_t BDT_inputs[NTRIPLETS][w3p_bdt::n_features];
    get_event_inputs_ref(selected, BDT_inputs);

    // Get BDT score for each triplet
    w3p_bdt::score_t BDT_scores[NTRIPLETS];
    get_event_scores_ref(BDT_inputs, BDT_scores);

    // Get highest BDT score among triplets
    get_highest_score_ref(BDT_scores, max_score);
}
