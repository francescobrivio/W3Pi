#include "best_seed.h"

Puppi BestSeed(Puppi a, bool a_masked, Puppi b, bool b_masked) {
    bool bestByPt = (a.hwPt >= b.hwPt);
    return (!a_masked && (bestByPt || b_masked)) ? a : b;
}

template<int width>
Puppi BestSeedReduce(const Puppi x[width], const bool masked[width]){
    // Tree reduce from https://github.com/definelicht/hlslib/blob/master/include/hlslib/xilinx/TreeReduce.h
    #pragma HLS inline
    static constexpr int halfWidth = width / 2;
    static constexpr int reducedSize = halfWidth + width % 2;
    Puppi reduced[reducedSize];
    bool maskreduced[reducedSize]; 
    #pragma HLS array_partition variable=reduced complete
    for(int i = 0; i < halfWidth; ++i) {
        #pragma HLS unroll
        reduced[i] = BestSeed(x[i*2], masked[i*2], x[i*2+1], masked[i*2+1]);
        maskreduced[i] = (masked[i*2] && masked[i*2+1]);
    }
    if(halfWidth != reducedSize){
        reduced[reducedSize - 1] = x[width - 1];
        maskreduced[reducedSize - 1] = masked[width - 1];
    }
    return BestSeedReduce<reducedSize>(reduced, maskreduced);
}

template<>
Puppi BestSeedReduce<2>(const Puppi x[2], const bool masked[2]){
    #pragma HLS inline
    return BestSeed(x[0], masked[0], x[1], masked[1]);
}

Puppi find_seed(const Puppi in[NPUPPI_MAX], const bool masked[NPUPPI_MAX]) {
    #pragma HLS inline off
	#pragma HLS pipeline ii=1
    return BestSeedReduce<NPUPPI_MAX>(in, masked);
}
