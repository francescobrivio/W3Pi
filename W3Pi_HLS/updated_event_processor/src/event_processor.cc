#include "event_processor.h"
#include "bitonic_hybrid.h"
#include "isolation.h"
#include "data.h"

// ------------------------------------------------------------------
// Masker: mask L1Puppi objects that don't pass selections
void masker (const Puppi input[NPUPPI_MAX], ap_uint<NPUPPI_MAX> & masked)
{
    #pragma HLS ARRAY_PARTITION variable=input complete

    // Clear masked array
    masked = 0;

    // Apply selections
    LOOP_MASKER_FILL: for (unsigned int i = 0; i < NPUPPI_MAX; i++)
    {
        #pragma HLS UNROLL
        //bool badPt  = (input[i].hwPt < 3.25 );
        bool badEta = (input[i].hwEta < -Puppi::ETA_CUT || input[i].hwEta > Puppi::ETA_CUT);
        bool badID  = (input[i].hwID < 2 || input[i].hwID > 5);
        //masked[i] = (badPt || badEta || badID);
        masked[i] = (badEta || badID);
    }
}

// ------------------------------------------------------------------
// Slimmer: replace masked candidates with empty/dummy puppi
void slimmer (const Puppi input[NPUPPI_MAX], const ap_uint<NPUPPI_MAX> masked, Puppi slimmed[NPUPPI_MAX])
{
    #pragma HLS ARRAY_PARTITION variable=input complete
    #pragma HLS ARRAY_PARTITION variable=slimmed complete

    Puppi dummy;
    dummy.clear();

    LOOP_SLIMMER_FILL: for (unsigned int i = 0; i < NPUPPI_MAX; i++)
    {
        #pragma HLS UNROLL
        slimmed[i] = masked[i] ? dummy : input[i];
    }
}

// ------------------------------------------------------------------
// Sort slimmed candidates by pT using bitonicSort from bitonic_hybrid.h
// Split and order 16 arrays of 13 candidates (16 x 13 = 208)
void orderer7 (const Puppi slimmed[NPUPPI_MAX],
               Puppi ordered1 [NSPLITS], Puppi ordered2 [NSPLITS],
               Puppi ordered3 [NSPLITS], Puppi ordered4 [NSPLITS],
               Puppi ordered5 [NSPLITS], Puppi ordered6 [NSPLITS],
               Puppi ordered7 [NSPLITS], Puppi ordered8 [NSPLITS],
               Puppi ordered9 [NSPLITS], Puppi ordered10[NSPLITS],
               Puppi ordered11[NSPLITS], Puppi ordered12[NSPLITS],
               Puppi ordered13[NSPLITS], Puppi ordered14[NSPLITS],
               Puppi ordered15[NSPLITS], Puppi ordered16[NSPLITS])
{
    #pragma HLS ARRAY_PARTITION variable=slimmed   complete
    #pragma HLS ARRAY_PARTITION variable=ordered1  complete
    #pragma HLS ARRAY_PARTITION variable=ordered2  complete
    #pragma HLS ARRAY_PARTITION variable=ordered3  complete
    #pragma HLS ARRAY_PARTITION variable=ordered4  complete
    #pragma HLS ARRAY_PARTITION variable=ordered5  complete
    #pragma HLS ARRAY_PARTITION variable=ordered6  complete
    #pragma HLS ARRAY_PARTITION variable=ordered7  complete
    #pragma HLS ARRAY_PARTITION variable=ordered8  complete
    #pragma HLS ARRAY_PARTITION variable=ordered9  complete
    #pragma HLS ARRAY_PARTITION variable=ordered10 complete
    #pragma HLS ARRAY_PARTITION variable=ordered11 complete
    #pragma HLS ARRAY_PARTITION variable=ordered12 complete
    #pragma HLS ARRAY_PARTITION variable=ordered13 complete
    #pragma HLS ARRAY_PARTITION variable=ordered14 complete
    #pragma HLS ARRAY_PARTITION variable=ordered15 complete
    #pragma HLS ARRAY_PARTITION variable=ordered16 complete

    // Split in 8 arrays
    LOOP_ORDERER7_SPLIT: for (unsigned int i=0; i < NSPLITS; i++)
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
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered1  , 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered2  , 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered3  , 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered4  , 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered5  , 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered6  , 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered7  , 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered8  , 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, false>::run(ordered9 , 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, false>::run(ordered10, 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, false>::run(ordered11, 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, false>::run(ordered12, 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, false>::run(ordered13, 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, false>::run(ordered14, 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, false>::run(ordered15, 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, false>::run(ordered16, 0);
}

// Split and order 8 arrays of 26 candidates (8 x 26 = 208)
void orderer7bis (const Puppi slimmed[NPUPPI_MAX],
                  Puppi ordered1 [NSPLITS], Puppi ordered2 [NSPLITS],
                  Puppi ordered3 [NSPLITS], Puppi ordered4 [NSPLITS],
                  Puppi ordered5 [NSPLITS], Puppi ordered6 [NSPLITS],
                  Puppi ordered7 [NSPLITS], Puppi ordered8 [NSPLITS])
{
    #pragma HLS ARRAY_PARTITION variable=slimmed   complete
    #pragma HLS ARRAY_PARTITION variable=ordered1  complete
    #pragma HLS ARRAY_PARTITION variable=ordered2  complete
    #pragma HLS ARRAY_PARTITION variable=ordered3  complete
    #pragma HLS ARRAY_PARTITION variable=ordered4  complete
    #pragma HLS ARRAY_PARTITION variable=ordered5  complete
    #pragma HLS ARRAY_PARTITION variable=ordered6  complete
    #pragma HLS ARRAY_PARTITION variable=ordered7  complete
    #pragma HLS ARRAY_PARTITION variable=ordered8  complete

    // Split in 8 arrays
    LOOP_ORDERER7_SPLIT: for (unsigned int i=0; i < NSPLITS; i++)
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
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered1, 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered2, 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered3, 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered4, 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered5, 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered6, 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered7, 0);
    hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered8, 0);
}

// Sort in pairs the split arrays
void orderer7f (const Puppi slimmed[NPUPPI_MAX], Puppi ordered[NSUBARR][NSPLITS])
{
    #pragma HLS ARRAY_PARTITION variable=slimmed complete
    #pragma HLS ARRAY_PARTITION variable=ordered complete dim=2

    // Split in 8 arrays
    LOOP_ORDERER7E_FILL: for (unsigned int i=0; i < NSPLITS; i++)
    {
        #pragma HLS UNROLL
        ordered[0][i] = slimmed[i];
        ordered[1][i] = slimmed[i+NSPLITS];
        ordered[2][i] = slimmed[i+2*NSPLITS];
        ordered[3][i] = slimmed[i+3*NSPLITS];
        ordered[4][i] = slimmed[i+4*NSPLITS];
        ordered[5][i] = slimmed[i+5*NSPLITS];
        ordered[6][i] = slimmed[i+6*NSPLITS];
        ordered[7][i] = slimmed[i+7*NSPLITS];
    }

    // Sort the 8 arrays
    LOOP_ORDERER7E_SORT: for (unsigned int i=0; i < NSUBARR / 2; i++)
    {
        hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered[i]  , 0);
        hybridBitonicSort::bitonicSorter<Puppi, NSPLITS, 0, true>::run(ordered[i+1], 0);
    }
}

// ------------------------------------------------------------------
// Utils for merge-sorting

// get_bitonic_sequence
void get_bitonic_sequenceA(const Puppi in1[NSPLITS], const Puppi in2[NSPLITS], Puppi bitonic[2*NSPLITS]) {
    //#pragma HLS inline
    #pragma HLS array_partition variable=in1 complete
    #pragma HLS array_partition variable=in2 complete
    #pragma HLS array_partition variable=bitonic complete
    make_bitonic_loop: for(int id = 0, ia=NSPLITS-1; id<NSPLITS; id++, ia--)
    {
        #pragma HLS UNROLL
        bitonic[id] = in1[ia];
        bitonic[id+(NSPLITS)] = in2[id];
    }
}

void get_bitonic_sequenceB(const Puppi in1[2*NSPLITS], const Puppi in2[2*NSPLITS], Puppi bitonic[4*NSPLITS]) {
    //#pragma HLS inline
    #pragma HLS array_partition variable=in1 complete
    #pragma HLS array_partition variable=in2 complete
    #pragma HLS array_partition variable=bitonic complete
    make_bitonic_loop: for(int id = 0, ia=2*NSPLITS-1; id<2*NSPLITS; id++, ia--)
    {
        #pragma HLS UNROLL
        bitonic[id] = in1[ia];
        bitonic[id+(2*NSPLITS)] = in2[id];
    }
}

void get_bitonic_sequenceC(const Puppi in1[4*NSPLITS], const Puppi in2[4*NSPLITS], Puppi bitonic[8*NSPLITS]) {
    //#pragma HLS inline
    #pragma HLS array_partition variable=in1 complete
    #pragma HLS array_partition variable=in2 complete
    #pragma HLS array_partition variable=bitonic complete
    make_bitonic_loop: for(int id = 0, ia=4*NSPLITS-1; id<4*NSPLITS; id++, ia--)
    {
        #pragma HLS UNROLL
        bitonic[id] = in1[ia];
        bitonic[id+(4*NSPLITS)] = in2[id];
    }
}

void get_bitonic_sequenceD(const Puppi in1[8*NSPLITS], const Puppi in2[8*NSPLITS], Puppi bitonic[16*NSPLITS]) {
    #pragma HLS inline
    #pragma HLS array_partition variable=in1 complete
    #pragma HLS array_partition variable=in2 complete
    #pragma HLS array_partition variable=bitonic complete
    make_bitonic_loop: for(int id = 0, ia=8*NSPLITS-1; id<8*NSPLITS; id++, ia--)
    {
        #pragma HLS UNROLL
        bitonic[id] = in1[ia];
        bitonic[id+(8*NSPLITS)] = in2[id];
    }
}

// merge_sort
void merge_sortA(const Puppi in1[NSPLITS], const Puppi in2[NSPLITS], Puppi sorted_out[2*NSPLITS]) {
    //#pragma HLS inline
    #pragma HLS array_partition variable=in1 complete
    #pragma HLS array_partition variable=in2 complete
    #pragma HLS array_partition variable=sorted_out complete
    get_bitonic_sequenceA(in1, in2, sorted_out);
    hybridBitonicSort::bitonicMerger<Puppi, 2*NSPLITS, 0>::run(sorted_out, 0);
}

void merge_sortB(const Puppi in1[2*NSPLITS], const Puppi in2[2*NSPLITS], Puppi sorted_out[4*NSPLITS]) {
    //#pragma HLS inline
    #pragma HLS array_partition variable=in1 complete
    #pragma HLS array_partition variable=in2 complete
    #pragma HLS array_partition variable=sorted_out complete
    get_bitonic_sequenceB(in1, in2, sorted_out);
    hybridBitonicSort::bitonicMerger<Puppi, 4*NSPLITS, 0>::run(sorted_out, 0);
}

void merge_sortC(const Puppi in1[4*NSPLITS], const Puppi in2[4*NSPLITS], Puppi sorted_out[8*NSPLITS]) {
    //#pragma HLS inline
    #pragma HLS array_partition variable=in1 complete
    #pragma HLS array_partition variable=in2 complete
    #pragma HLS array_partition variable=sorted_out complete
    get_bitonic_sequenceC(in1, in2, sorted_out);
    hybridBitonicSort::bitonicMerger<Puppi, 8*NSPLITS, 0>::run(sorted_out, 0);
}

void merge_sortD(const Puppi in1[8*NSPLITS], const Puppi in2[8*NSPLITS], Puppi sorted_out[16*NSPLITS]) {
    #pragma HLS inline
    #pragma HLS array_partition variable=in1 complete
    #pragma HLS array_partition variable=in2 complete
    #pragma HLS array_partition variable=sorted_out complete
    get_bitonic_sequenceD(in1, in2, sorted_out);
    hybridBitonicSort::bitonicMerger<Puppi, 16*NSPLITS, 0>::run(sorted_out, 0);
}

// ------------------------------------------------------------------
// Merger: merge ordered arrays with bitonicMerger
void merger (Puppi ordered1 [NSPLITS], Puppi ordered2 [NSPLITS],
             Puppi ordered3 [NSPLITS], Puppi ordered4 [NSPLITS],
             Puppi ordered5 [NSPLITS], Puppi ordered6 [NSPLITS],
             Puppi ordered7 [NSPLITS], Puppi ordered8 [NSPLITS],
             Puppi ordered9 [NSPLITS], Puppi ordered10[NSPLITS],
             Puppi ordered11[NSPLITS], Puppi ordered12[NSPLITS],
             Puppi ordered13[NSPLITS], Puppi ordered14[NSPLITS],
             Puppi ordered15[NSPLITS], Puppi ordered16[NSPLITS],
             Puppi merged[NPUPPI_MAX])
{
    #pragma HLS ARRAY_PARTITION variable=ordered1  complete
    #pragma HLS ARRAY_PARTITION variable=ordered2  complete
    #pragma HLS ARRAY_PARTITION variable=ordered3  complete
    #pragma HLS ARRAY_PARTITION variable=ordered4  complete
    #pragma HLS ARRAY_PARTITION variable=ordered5  complete
    #pragma HLS ARRAY_PARTITION variable=ordered6  complete
    #pragma HLS ARRAY_PARTITION variable=ordered7  complete
    #pragma HLS ARRAY_PARTITION variable=ordered8  complete
    #pragma HLS ARRAY_PARTITION variable=ordered9  complete
    #pragma HLS ARRAY_PARTITION variable=ordered10 complete
    #pragma HLS ARRAY_PARTITION variable=ordered11 complete
    #pragma HLS ARRAY_PARTITION variable=ordered12 complete
    #pragma HLS ARRAY_PARTITION variable=ordered13 complete
    #pragma HLS ARRAY_PARTITION variable=ordered14 complete
    #pragma HLS ARRAY_PARTITION variable=ordered15 complete
    #pragma HLS ARRAY_PARTITION variable=ordered16 complete
    #pragma HLS ARRAY_PARTITION variable=merged    complete

    // Recursively merge and sort in pairs
    // Step 1
    Puppi merge1[2*NSPLITS];
    merge_sortA(ordered1, ordered2, merge1);

    Puppi merge2[2*NSPLITS];
    merge_sortA(ordered3, ordered4, merge2);

    Puppi merge3[2*NSPLITS];
    merge_sortA(ordered5, ordered6, merge3);

    Puppi merge4[2*NSPLITS];
    merge_sortA(ordered7, ordered8, merge4);

    Puppi merge5[2*NSPLITS];
    merge_sortA(ordered9, ordered10, merge5);

    Puppi merge6[2*NSPLITS];
    merge_sortA(ordered11, ordered12, merge6);

    Puppi merge7[2*NSPLITS];
    merge_sortA(ordered13, ordered14, merge7);

    Puppi merge8[2*NSPLITS];
    merge_sortA(ordered15, ordered16, merge8);

    // Step 2
    Puppi merge9[4*NSPLITS];
    merge_sortB(merge1, merge2, merge9);

    Puppi merge10[4*NSPLITS];
    merge_sortB(merge3, merge4, merge10);

    Puppi merge11[4*NSPLITS];
    merge_sortB(merge5, merge6, merge11);

    Puppi merge12[4*NSPLITS];
    merge_sortB(merge7, merge8, merge12);

    // Step 3
    Puppi merge13[8*NSPLITS];
    merge_sortC(merge9, merge10, merge13);

    Puppi merge14[8*NSPLITS];
    merge_sortC(merge11, merge12, merge14);

    // Final
    merge_sortD(merge13, merge14, merged);
}

// Same as merger but using 8 arrays of 26 elements
void merger7bis (Puppi ordered1 [NSPLITS], Puppi ordered2 [NSPLITS],
                 Puppi ordered3 [NSPLITS], Puppi ordered4 [NSPLITS],
                 Puppi ordered5 [NSPLITS], Puppi ordered6 [NSPLITS],
                 Puppi ordered7 [NSPLITS], Puppi ordered8 [NSPLITS],
                 Puppi merged[NPUPPI_MAX])
{
    #pragma HLS ARRAY_PARTITION variable=ordered1  complete
    #pragma HLS ARRAY_PARTITION variable=ordered2  complete
    #pragma HLS ARRAY_PARTITION variable=ordered3  complete
    #pragma HLS ARRAY_PARTITION variable=ordered4  complete
    #pragma HLS ARRAY_PARTITION variable=ordered5  complete
    #pragma HLS ARRAY_PARTITION variable=ordered6  complete
    #pragma HLS ARRAY_PARTITION variable=ordered7  complete
    #pragma HLS ARRAY_PARTITION variable=ordered8  complete
    #pragma HLS ARRAY_PARTITION variable=merged    complete

    // Recursively merge and sort in pairs
    // Step 1
    Puppi merge1[2*NSPLITS];
    merge_sortA(ordered1, ordered2, merge1);

    Puppi merge2[2*NSPLITS];
    merge_sortA(ordered3, ordered4, merge2);

    Puppi merge3[2*NSPLITS];
    merge_sortA(ordered5, ordered6, merge3);

    Puppi merge4[2*NSPLITS];
    merge_sortA(ordered7, ordered8, merge4);

    // Step 2
    Puppi merge5[4*NSPLITS];
    merge_sortB(merge1, merge2, merge5);

    Puppi merge6[4*NSPLITS];
    merge_sortB(merge3, merge4, merge6);

    // Final
    merge_sortC(merge5, merge6, merged);
}

// Same as merger7bis but starting from orderer7f
void merger7f (Puppi ordered[NSUBARR][NSPLITS], Puppi merged[NPUPPI_MAX])
{
    #pragma HLS ARRAY_PARTITION variable=ordered complete dim=2
    #pragma HLS ARRAY_PARTITION variable=merged  complete

    // original                                                       // --> v0
    //#pragma HLS ALLOCATION instances=merge_sortA limit=1 function   // --> v1 |
    //#pragma HLS ALLOCATION instances=merge_sortA limit=2 function   // --> v2 |
    //#pragma HLS ALLOCATION instances=merge_sortB limit=1 function   // --> v3 |-> all the same v1-4
    //#pragma HLS ALLOCATION instances=merge_sortA limit=1 function   // --> v4 |
    //#pragma HLS ALLOCATION instances=merge_sortB limit=1 function   //     v4 |
    // remove inlinings from merge_sortA                              // --> v5
    //#pragma HLS ALLOCATION instances=merge_sortA limit=1 function   // + remove inlinings from merge_sortA --> v6 // same as v5
    // remove inlinings from merge_sortB                              // --> v7
    // remove inlinings from merge_sortA and merge_sortB              // --> v8
    // remove inlinings from merge_sortC                              // --> v9  // almost identical to v1-4
    // remove inlinings from merge_sortA, merge_sortB and merge_sortC // --> v10 // almost identical to v8
    // remove inlinings from all merge_sort and all bitonic_sequence  // --> v11 // best so far
    // next:
    //   - test not inlining + ALLOCATION combination

    // Recursively merge and sort in pairs
    // Step 1
    Puppi merge1[2*NSPLITS];
    merge_sortA(ordered[0], ordered[1], merge1);

    Puppi merge2[2*NSPLITS];
    merge_sortA(ordered[2], ordered[3], merge2);

    Puppi merge3[2*NSPLITS];
    merge_sortA(ordered[4], ordered[5], merge3);

    Puppi merge4[2*NSPLITS];
    merge_sortA(ordered[6], ordered[7], merge4);

    // Step 2
    Puppi merge5[4*NSPLITS];
    merge_sortB(merge1, merge2, merge5);

    Puppi merge6[4*NSPLITS];
    merge_sortB(merge3, merge4, merge6);

    // Final
    merge_sortC(merge5, merge6, merged);
}

// ------------------------------------------------------------------
// Select first NPUPPI_SEL from the ordered list
void selector(const Puppi merged[NPUPPI_MAX], Puppi selected[NPUPPI_SEL])
{
    #pragma HLS ARRAY_PARTITION variable=merged complete
    #pragma HLS ARRAY_PARTITION variable=selected complete

    LOOP_SELECTOR: for (unsigned int i = 0; i < NPUPPI_SEL; i++)
    {
        #pragma HLS UNROLL
        selected[i] = merged[i];
    }
}


// ------------------------------------------------------------------
// Get maximum deltaVz
Puppi::z0_t get_max_dVz (Puppi::z0_t z0, Puppi::z0_t z1, Puppi::z0_t z2)
{
    #pragma HLS inline

    // Get all dVz values
    Puppi::z0_t dVz_01 = z0 - z1;
    Puppi::z0_t dVz_02 = z0 - z2;
    Puppi::z0_t dVz_12 = z1 - z2;

    // std library max
    Puppi::z0_t max = std::max(dVz_01, std::max(dVz_02, dVz_12));

    return max;
}

// ------------------------------------------------------------------
// DeltaR methods
inline dr2_t drToHwDr2 (float dr)
{
    return dr2_t(round(std::pow(dr/Puppi::ETAPHI_LSB,2)));
}

inline float hwDr2ToDr (dr2_t dr2)
{
    return std::sqrt(dr2.to_int()*Puppi::ETAPHI_LSB*Puppi::ETAPHI_LSB);
}

ap_int<Puppi::eta_t::width+1> deltaEta (Puppi::eta_t eta1, Puppi::eta_t eta2)
{
    #pragma HLS latency min=1
    #pragma HLS inline off
    return eta1 - eta2;
}

inline dr2_t deltaR2 (const Puppi & p1, const Puppi & p2)
{
    auto dphi = p1.hwPhi - p2.hwPhi;
    if (dphi > Puppi::INT_PI) dphi -= Puppi::INT_2PI;
    else if (dphi < -Puppi::INT_PI) dphi += Puppi::INT_2PI;
    auto deta = p1.hwEta - p2.hwEta;
    return dphi*dphi + deta*deta;
}

inline dr2_t deltaR2_slow (const Puppi & p1, const Puppi & p2)
{
    auto dphi = p1.hwPhi - p2.hwPhi;
    if (dphi > Puppi::INT_PI) dphi -= Puppi::INT_2PI;
    else if (dphi < -Puppi::INT_PI) dphi += Puppi::INT_2PI;
    auto deta = deltaEta(p1.hwEta, p2.hwEta);
    return dphi*dphi + deta*deta;
}

// Get min deltaR slow
dr2_t get_min_deltaR2_slow (const Puppi p0, const Puppi p1, const Puppi p2)
{
    // Get the 3 dR values
    dr2_t dr2_01 = deltaR2_slow(p0, p1);
    dr2_t dr2_02 = deltaR2_slow(p0, p2);
    dr2_t dr2_12 = deltaR2_slow(p1, p2);

    // std library min
    dr2_t min = std::min(dr2_01, std::min(dr2_02, dr2_12));

    return min;
}

// Get min deltaR
dr2_t get_min_deltaR2 (const Puppi p0, const Puppi p1, const Puppi p2)
{
    // Get the 3 dR values
    dr2_t dr2_01 = deltaR2(p0, p1);
    dr2_t dr2_02 = deltaR2(p0, p2);
    dr2_t dr2_12 = deltaR2(p1, p2);

    // std library min
    dr2_t min = std::min(dr2_01, std::min(dr2_02, dr2_12));

    return min;
}

// ------------------------------------------------------------------
// Trigonometric methods
// Declare and fill cosine LUT
void _lut_cos_init(cos_t table_cos[COSCOSH_LUT_SIZE])
{
    for (int i = 0; i < COSCOSH_LUT_SIZE; ++i)
    {
        // Get float value of phi
        float alpha = Puppi::ETAPHI_LSB * i;

        // Get cos values (times LSB)
        int ic = cos(alpha) * COSCOSH_LSB;

        // Fix boundaries depending on (half) table size
        if (ic >  (COSCOSH_LUT_SIZE/2-1)) ic =  (COSCOSH_LUT_SIZE/2-1);
        if (ic < -(COSCOSH_LUT_SIZE/2-1)) ic = -(COSCOSH_LUT_SIZE/2-1);

        // Fill table
        table_cos[i] = ic;
    }
}

// Get cos value from LUT
cos_t get_cos_phi (Puppi::phi_t phi)
{
    // Declare and fill cos LUT
    cos_t _table_cos[COSCOSH_LUT_SIZE];
    _lut_cos_init(_table_cos);

    // Change phi to only positive value
    int iphi = phi;
    if (phi < 0) iphi = -iphi;

    // Return cos phi value
    return _table_cos[iphi];
}

// Declare and fill hyperbolic cosine LUT
void _lut_cosh_init(cosh_t table_cosh[COSCOSH_LUT_SIZE])
{
    for (int i = 0; i < COSCOSH_LUT_SIZE; ++i)
    {
        // Get float value of eta
        float alpha = Puppi::ETAPHI_LSB * i;

        // Get cosh values (times LSB)
        int ich = cosh(alpha) * COSCOSH_LSB;

        // Fix upper boundary
        if (ich > (COSCOSH_LUT_SIZE-1)) ich = (COSCOSH_LUT_SIZE-1);

        // Fill table
        table_cosh[i] = ich;
    }
}

// Get cosh value from LUT
cosh_t get_cosh_eta (Puppi::eta_t eta)
{
    // Declare and fill cos LUT
    cosh_t _table_cosh[COSCOSH_LUT_SIZE];
    _lut_cosh_init(_table_cosh);

    // Change eta to only positive value
    int ieta = eta;
    if (eta < 0) ieta = -ieta;

    // Return cosh eta value
    return _table_cosh[ieta];
}

// ------------------------------------------------------------------
// Invariant mass quared of pair of particles
// m^2 = 2 * pT1 * pT2 * ( cosh(eta1 - eta2) - cos(phi1 - phi2) )
// FIXME: this can possibly return a negative value??
mass_t get_pair_mass (const Puppi & p1, const Puppi & p2)
{
    // Get dPhi and dEta
    auto dphi = p1.hwPhi - p2.hwPhi;
    if (dphi > Puppi::INT_PI) dphi -= Puppi::INT_2PI;
    else if (dphi < -Puppi::INT_PI) dphi += Puppi::INT_2PI;
    auto deta = p1.hwEta - p2.hwEta;

    // Get cos and cosh
    cos_t cosdPhi   = get_cos_phi(dphi);
    cosh_t coshdEta = get_cosh_eta(deta);

    // Compute and return final mass
    mass_t pair_m = 2 * p1.hwPt * p2.hwPt * (coshdEta - cosdPhi);

    return pair_m;
}

// ------------------------------------------------------------------
// Prepare triplet inputs features:
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
void get_triplet_inputs (const Puppi selected[NPUPPI_SEL], idx_t idx0, idx_t idx1, idx_t idx2, w3p_bdt::input_t BDT_inputs[w3p_bdt::n_features])
{
    // Fill BDT input for triplet
    BDT_inputs[0]  = selected[idx2].hwPt;
    BDT_inputs[1]  = selected[idx1].hwPt;
    BDT_inputs[2]  = get_pair_mass(selected[idx0], selected[idx1]);
    BDT_inputs[3]  = selected[idx0].charge() + selected[idx1].charge() + selected[idx2].charge();
    BDT_inputs[4]  = selected[idx0].hwZ0 - selected[idx2].hwZ0;
    BDT_inputs[5]  = selected[idx0].hwPt;
    BDT_inputs[6]  = selected[idx0].hwPt + selected[idx1].hwPt + selected[idx2].hwPt;
    BDT_inputs[7]  = get_pair_mass(selected[idx0], selected[idx2]);
    BDT_inputs[8]  = get_max_dVz(selected[idx0].hwZ0, selected[idx1].hwZ0, selected[idx2].hwZ0);
    BDT_inputs[9]  = get_min_deltaR2_slow(selected[idx0], selected[idx1], selected[idx2]);
    BDT_inputs[10] = selected[idx2].hwEta;
}

// ------------------------------------------------------------------
// Get all event inputs
// - 8 triplets from: 5 from 1st+2nd, 2 from 1st+3rd and 1 from 2nd+3rd:
// (0,1,2)-(0,1,3)-(0,1,4)-(0,1,5)-(0,1,6)-(0,2,3)-(0,2,3)-(1,2,3)
void get_event_inputs (const Puppi selected[NPUPPI_SEL], w3p_bdt::input_t BDT_inputs[NTRIPLETS][w3p_bdt::n_features])
{
    #pragma HLS ARRAY_PARTITION variable=selected complete
    #pragma HLS ARRAY_PARTITION variable=BDT_inputs complete dim=0

    // 5 triplets from 1st+2nd pivots
    LOOP_EVENT_INPUTS: for (unsigned int i = 0; i < NTRIPLETS-3; i++)
    {
        #pragma HLS unroll
        idx_t third_idx = i + 2; // 2 pivots excluded
        get_triplet_inputs(selected, 0, 1, i, BDT_inputs[i]);
    }

    // 2 triplets from 1st+3rd pivots
    get_triplet_inputs(selected, 0, 2, 3, BDT_inputs[5]);
    get_triplet_inputs(selected, 0, 2, 4, BDT_inputs[6]);

    // 1 triplet from 2nd+3rd pivots
    get_triplet_inputs(selected, 1, 2, 3, BDT_inputs[7]);
}

// ------------------------------------------------------------------
// Get BDT scores of the selected triplets
void get_event_scores (w3p_bdt::input_t BDT_inputs[NTRIPLETS][w3p_bdt::n_features], w3p_bdt::score_t BDT_scores[NTRIPLETS])
{
    #pragma HLS ARRAY_PARTITION variable=BDT_inputs complete dim=0
    #pragma HLS ARRAY_PARTITION variable=BDT_scores complete

    // Get the BDT scores
    LOOP_EVENT_SCORES: for (unsigned int i = 0; i < NTRIPLETS; i++)
    {
        #pragma HLS unroll
        w3p_bdt::bdt.decision_function(BDT_inputs[i], &BDT_scores[i]);
    }
}

// ------------------------------------------------------------------
// Get highest BDT score
void get_highest_score (w3p_bdt::score_t BDT_scores[NTRIPLETS], w3p_bdt::score_t & high_score)
{
    #pragma HLS ARRAY_PARTITION variable=BDT_scores complete

    // Sort BDT scores in descending order
    hybridBitonicSort::bitonicSorter<w3p_bdt::score_t, NTRIPLETS, 0, true>::run(BDT_scores, 0);

    // Get the highest score
    high_score = BDT_scores[0];
}

// ------------------------------------------------------------------
// Full EventProcessor
void EventProcessor (const Puppi input[NPUPPI_MAX], w3p_bdt::score_t & max_score)
{
    #pragma HLS ARRAY_PARTITION variable=input complete

    // Mask candidates
    ap_uint<NPUPPI_MAX> masked;
    masker(input, masked);

    // Replace masked candidates with dummy
    Puppi slimmed[NPUPPI_MAX];
    slimmer(input, masked, slimmed);

    // Split in arrays and sort them according to pT
    Puppi ordered1[NSPLITS] , ordered2[NSPLITS] , ordered3[NSPLITS] , ordered4[NSPLITS] ,
          ordered5[NSPLITS] , ordered6[NSPLITS] , ordered7[NSPLITS] , ordered8[NSPLITS] ,
          ordered9[NSPLITS] , ordered10[NSPLITS], ordered11[NSPLITS], ordered12[NSPLITS],
          ordered13[NSPLITS], ordered14[NSPLITS], ordered15[NSPLITS], ordered16[NSPLITS];
    orderer7(slimmed, 
             ordered1 , ordered2 , ordered3 , ordered4 ,
             ordered5 , ordered6 , ordered7 , ordered8 ,
             ordered9 , ordered10, ordered11, ordered12,
             ordered13, ordered14, ordered15, ordered16);
    
    // Merge the sorted split-arrays
    Puppi merged[NPUPPI_MAX];
    merger(ordered1, ordered2, ordered3, ordered4,
           ordered5, ordered6, ordered7, ordered8,
           ordered9, ordered10, ordered11, ordered12,
           ordered13, ordered14, ordered15, ordered16,
           merged);

    // Select only highest pT ordered-candidates
    Puppi selected[NPUPPI_SEL];
    #pragma HLS ARRAY_PARTITION variable=selected complete
    selector(merged, selected);

    // Get inputs for each triplet
    w3p_bdt::input_t BDT_inputs[NTRIPLETS][w3p_bdt::n_features];
    get_event_inputs(selected, BDT_inputs);

    // Get BDT score for each triplet
    w3p_bdt::score_t BDT_scores[NTRIPLETS];
    get_event_scores(BDT_inputs, BDT_scores);

    // Get highest BDT score among triplets
    get_highest_score(BDT_scores, max_score);
}

// EventProcessor7bis - same as event processor, but with 8 arrays
void EventProcessor7bis (const Puppi input[NPUPPI_MAX], w3p_bdt::score_t & max_score)
{
    #pragma HLS ARRAY_PARTITION variable=input complete

    // Mask candidates
    ap_uint<NPUPPI_MAX> masked;
    masker(input, masked);

    // Replace masked candidates with dummy
    Puppi slimmed[NPUPPI_MAX];
    slimmer(input, masked, slimmed);

    // Split in arrays and sort them according to pT
    Puppi ordered1[NSPLITS], ordered2[NSPLITS], ordered3[NSPLITS], ordered4[NSPLITS],
          ordered5[NSPLITS], ordered6[NSPLITS], ordered7[NSPLITS], ordered8[NSPLITS];
    orderer7bis(slimmed,
                ordered1, ordered2, ordered3, ordered4,
                ordered5, ordered6, ordered7, ordered8);
    
    // Merge the sorted split-arrays
    Puppi merged[NPUPPI_MAX];
    merger7bis(ordered1, ordered2, ordered3, ordered4,
               ordered5, ordered6, ordered7, ordered8,
               merged);

    // Select only highest pT ordered-candidates
    Puppi selected[NPUPPI_SEL];
    #pragma HLS ARRAY_PARTITION variable=selected complete
    selector(merged, selected);

    // Get inputs for each triplet
    w3p_bdt::input_t BDT_inputs[NTRIPLETS][w3p_bdt::n_features];
    get_event_inputs(selected, BDT_inputs);

    // Get BDT score for each triplet
    w3p_bdt::score_t BDT_scores[NTRIPLETS];
    get_event_scores(BDT_inputs, BDT_scores);

    // Get highest BDT score among triplets
    get_highest_score(BDT_scores, max_score);
}

// EventProcessor7f - same as EventProcessor7bis, but using orderer7f and merger7f
void EventProcessor7f (const Puppi input[NPUPPI_MAX], w3p_bdt::score_t & max_score)
{
    #pragma HLS ARRAY_PARTITION variable=input complete

    // Mask candidates
    ap_uint<NPUPPI_MAX> masked;
    masker(input, masked);

    // Replace masked candidates with dummy
    Puppi slimmed[NPUPPI_MAX];
    slimmer(input, masked, slimmed);

    // Split in arrays and sort them according to pT
    // Make 2D array of [NSUBARR]x[NSPLITS] = [8]x[26] = 208 for ordering and merging
    Puppi ordered[NSUBARR][NSPLITS];
    orderer7f(slimmed, ordered);
    
    // Merge the sorted split-arrays
    Puppi merged[NPUPPI_MAX];
    merger7f(ordered, merged);

    // Select only highest pT ordered-candidates
    Puppi selected[NPUPPI_SEL];
    #pragma HLS ARRAY_PARTITION variable=selected complete
    selector(merged, selected);

    // Calculate isolation
    Puppi::pt_t selected_iso[NPUPPI_SEL];
    #pragma HLS ARRAY_PARTITION variable=selected_iso complete
    calculate_iso(selected, input, selected_iso);

    // Get inputs for each triplet
    w3p_bdt::input_t BDT_inputs[NTRIPLETS][w3p_bdt::n_features];
    get_event_inputs(selected, BDT_inputs);

    // Get BDT score for each triplet
    w3p_bdt::score_t BDT_scores[NTRIPLETS];
    get_event_scores(BDT_inputs, BDT_scores);

    // Get highest BDT score among triplets
    get_highest_score(BDT_scores, max_score);
}
