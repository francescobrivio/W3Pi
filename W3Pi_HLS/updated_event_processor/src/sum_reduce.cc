#include "sum_reduce.h"

template<int N>
Puppi::pt_t SumReduce(const Puppi::pt_t in[N]) {
       return SumReduce<N/2>(in) + SumReduce<N-N/2>(&in[N/2]);
}
template<>
Puppi::pt_t SumReduce<1>(const Puppi::pt_t in[1]) {
    return in[0];
}

Puppi::pt_t SumReduceAll(const Puppi::pt_t in[NPUPPI_MAX]) {
    #pragma HLS inline off
    #pragma HLS pipeline ii=1
    return SumReduce<NPUPPI_MAX>(in);
}