#ifndef PROVE_H
#define PROVE_H

#include <cstddef>
#include <memory>
#include <utility>
#include "sign.h"
#include "range_proof.h"
#include <mcl/bls12_381.hpp>


struct publicParams
{
    BBSParams& params;
    BBS::PublicKey &pk;
    size_t N;
    mcl::G1 gepoch;
    mcl::G2 l1,l2,l3;
    publicParams(BBSParams& params,BBS::PublicKey &pk,
                    size_t N): params(params), pk(pk), N(N)
                    {
                        mcl::hashAndMapToG1(gepoch, "0");
                        mcl::hashAndMapToG2(l1, "1");
                        mcl::hashAndMapToG2(l2, "2");
                        mcl::hashAndMapToG2(l3, "3");
                            
                    }

};

class ACProof
{
public:
    explicit ACProof(publicParams& pp);
    ~ACProof();
    struct statement{
        mcl::G1 tag;
        mcl::G2 vk;
        mcl::G2 T;
        mcl::G1 T_prime;
    };
    struct witness{
        BBS::Signature sig;
        mcl::Fr sid;
        mcl::Fr m_r;
        mcl::Fr ctr;
        mcl::Fr r_T;
        mcl::Fr r_vk;
    };
    struct com{
        mcl::G1 A_bar,B_bar,U_B,G1_ctr;
        mcl::G2 vk_bar,U_1,U_2,G2_ctr;
        mcl::GT U_Q;
    };
    struct resp{
        mcl::Fr beta_r,beta_e,beta_v,beta_sid,beta_Q,beta_bf,beta_mr,beta_ctr,beta_rT;
    };
    struct proof{
        com cm;
        resp rp;
        RangeProof proof;
    };
    void GenMatrials(std::vector<mcl::Fr>& messages,BBS::Signature &sig,statement &st, witness &wt);
    proof Prove(const statement& st, const witness& wt);
    bool Verify(const statement& st, const proof& pf);
    
private:
    publicParams& pp_;
    mcl::Fr k_from_hash(const statement&st, const com &cm);
};




#endif
