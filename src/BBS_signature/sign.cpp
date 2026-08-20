#include "sign.h"
#include <stdexcept>
#include <string>
#include <utility>
#include <cstdio>


//BBS implementation

BBS::BBS(BBSParams& params) : params_(params) {}

BBS::~BBS() = default;



void BBS::generate_keys(PublicKey& pk, SecretKey& sk) {
    const size_t message_count = params_.message_count;
    const mcl::G1& g1 = params_.g1;
    const mcl::G2& g2 = params_.g2;

    // 1. Pick x (Scalar field element)
    sk.x.setByCSPRNG();

    // 2. Generate pk = g2^x
    mcl::G2 X;
    mcl::G2::mul(X, g2, sk.x);
    pk.pk_g2 = X;
    std::vector<mcl::G1> bbs_bases(message_count + 1);
    for (mcl::G1& point : bbs_bases) {
        mcl::Fr scalar;
        scalar.setByCSPRNG();
        mcl::G1::mul(point, g1, scalar);
    }
    pk.pks_g1 = move(bbs_bases);
    //impl_->signature_public_key_g1 = move(bbs_bases);
}

mcl::G1 BBS::compute_B(const std::vector<mcl::Fr>& messages, const PublicKey& pk) {
    
    mcl::G1 B = params_.g1;

    for (size_t i = 0; i < messages.size(); ++i) {
        mcl::G1 term;
        // 执行点乘：term = h_i * m_i
        mcl::G1::mul(term, pk.pks_g1[i], messages[i]);
        // 执行点加：B = B + term
        B += term;
    }
    return B;
}

BBS::Signature BBS::sign(const std::vector<mcl::Fr>& messages, const PublicKey& pk, const SecretKey& sk) {
    mcl::Fr e;
    e.setByCSPRNG();


    mcl::G1 B = compute_B(messages, pk);

    // Compute A = B^{1 / (x + e)}
    mcl::Fr exp = sk.x + e;
    mcl::Fr::inv(exp,exp);

    mcl::G1 A;
    mcl::G1::mul(A, B, exp);

    return Signature(std::move(A), std::move(e));
}

bool BBS::verify(const std::vector<mcl::Fr>& messages, const Signature& sig, const PublicKey& pk) {
    mcl::G1 B = compute_B(messages, pk);
    mcl::G2 g2_e;
    mcl::G2::mul(g2_e, params_.g2, sig.e);
    mcl::G2 right_side = pk.pk_g2 + g2_e;
    mcl::GT pair1, pair2;
    mcl::pairing(pair1, sig.A, right_side);
    mcl::pairing(pair2, B, params_.g2);

    return pair1 == pair2;
}