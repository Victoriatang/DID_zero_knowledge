#ifndef RANGEPROOF_H
#define RANGEPROOF_H

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>
#include <mcl/bls12_381.hpp>

using namespace mcl;



struct InnerProductArg {
    std::vector<G1> L;
    std::vector<G1> R;
    Fr a_tag;
    Fr b_tag;
};

struct RangeProof {
    G1 A, S, T1, T2;
    Fr tau_x, miu, tx;
    InnerProductArg inner_product_proof;
};

//compute the inner product
Fr inner_product(const std::vector<Fr>& a, const std::vector<Fr>& b);

// compute the challenge value
Fr get_x_challenge(const G1& L, const G1& R, const G1& ux);

InnerProductArg prove(
    const std::vector<G1>& G,
    const std::vector<G1>& H,
    const G1& ux,
    const G1& P,
    const std::vector<Fr>& a,
    const std::vector<Fr>& b,
    std::vector<G1>& L_vec,
    std::vector<G1>& R_vec
);
bool fast_verify(
    const InnerProductArg& proof,
    const std::vector<G1>& g_vec,
    const std::vector<G1>& hi_tag,
    const G1& ux,
    const G1& P
);

std::vector<Fr> power_vector(Fr y, size_t n);
RangeProof prove_range(
    const std::vector<G1>& g_vec_in,
    const std::vector<G1>& h_vec_in,
    const G1& G, 
    const G1& H, 
    std::vector<Fr> secret, 
    const std::vector<Fr>& blinding, 
    size_t bit_length
);
bool fast_verify_range(
    const RangeProof& proof,
    const std::vector<G1>& g_vec,
    const std::vector<G1>& h_vec,
    const G1& G,
    const G1& H,
    const std::vector<G1>& ped_com, 
    size_t bit_length
);

#endif