#ifndef SIGN_H
#define SIGN_H

#include <cstddef>
#include <memory>
#include <utility>

#include <mcl/bls12_381.hpp>
struct BBSParams{
    mcl::CurveParam curve_param;
    size_t message_count;
    mcl::G1 g1;
    mcl::G2 g2;
    BBSParams(mcl::CurveParam curve_param,size_t message_count)
    {
        if (message_count == 0) {
            throw std::invalid_argument("message count must be > 0");
        }

        initPairing(curve_param);
        mcl::hashAndMapToG1(g1, "1");
        mcl::hashAndMapToG2(g2, "1");
    }
};



class BBS {
public:
    explicit BBS(BBSParams& params);
    ~BBS();

    
    struct PublicKey {
        mcl::G2 pk_g2;
        std::vector<mcl::G1> pks_g1;
    };

    struct SecretKey {
        mcl::Fr x;               
    };

    struct Signature {
       mcl::G1 A;
       mcl::Fr e;
    
       Signature() = default;
       Signature(mcl::G1 A, mcl::Fr e): A(std::move(A)), e(std::move(e)) {}
    };

    void generate_keys(PublicKey& pk, SecretKey& sk);

    
    Signature sign(const std::vector<mcl::Fr>& messages, const PublicKey& pk, const SecretKey& sk);

    bool verify(const std::vector<mcl::Fr>& messages, const Signature& sig, const PublicKey& pk);

private:   
    BBSParams &params_;
    // Compute B = g1 + \sum(h_i^{m_i})
    mcl::G1 compute_B(const std::vector<mcl::Fr>& messages, const PublicKey& pk);
};
#endif