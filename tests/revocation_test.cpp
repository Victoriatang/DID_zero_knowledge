#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <chrono>
#include <mcl/bls12_381.hpp>
#include <mcl/gmp_util.hpp>
#include <random>
#include "BBS_signature/threshold_CL.h"
#include "BBS_signature/prove_CL.h"

struct Timer {
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::string name;

    Timer(const std::string& n) : name(n) {
        start = std::chrono::high_resolution_clock::now();
    }
    ~Timer() {
        auto end = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(end - start).count();
        printf("[%s] %.3f ms\n", name.c_str(), ms);
    }
};
int main() {
    initPairing(mcl::BLS12_381);
    mcl::G1 G,H,Hs;
    mcl::hashAndMapToG1(G, "generator_G");
    mcl::hashAndMapToG1(H,"generator_H");
    mcl::hashAndMapToG1(Hs,"generator_Hs");
    mcl::G1 result,commit,temp;
    mcl::Fr x,r;
    x.setByCSPRNG();
    r.setByCSPRNG();

    {
    Timer t_total("Computing null");
    mcl::G1::mul(result,Hs,x);
    mcl::G1::mul(commit,G,x);
    mcl::G1::mul(temp,H,r);
    mcl::G1::add(commit,commit,temp);

    //proving
    mcl::Fr r0,x0;
    x0.setByCSPRNG();
    r0.setByCSPRNG();
    mcl::G1 result0,commit0,temp0;
    mcl::G1::mul(result0,Hs,x0);
    mcl::G1::mul(commit0,G,x0);
    mcl::G1::mul(temp0,H,r0);
    mcl::G1::add(commit0,commit0,temp0);
    mcl::Fr challenge;
    std::ostringstream output;
    output<<result0<<commit0<<result<<commit;
    challenge.setHashOf(output.str());

    mcl::Fr zr,zx;
    mcl::Fr::mul(zr,challenge,r);
    mcl::Fr::add(zr,zr,r0);
    mcl::Fr::mul(zx,challenge,x);
    mcl::Fr::add(zx,zx,x0);


    //verify
    for(int i = 0; i<20; i++)
    {
        std::ostringstream output0;
        output0<<result0<<commit0<<result<<commit;
        mcl::Fr challenge0;
        challenge0.setHashOf(output0.str());
        mcl::G1 lft,rft,temp;
        mcl::G1::mul(temp,G,zx);
        mcl::G1::mul(lft,H,zr);
        mcl::G1::add(lft,lft,temp);

        mcl::G1::mul(rft,commit,challenge0);
        mcl::G1::add(rft,rft,commit0);

        if(rft != lft)
            std::cout<<"checking the first equation failed"<<std::endl;
        
        mcl::G1::mul(lft,Hs,zx);
        mcl::G1::mul(rft,result,challenge0);
        mcl::G1::add(rft,rft,result0);

        if(rft != lft)
            std::cout<<"checking the second equation failed"<<std::endl;
    }
    mcl::G1 final = result;
    for(int i = 1; i<20; i++)
    {
        mcl::G1::add(final,final,result);
    }
}
  
        BICYCL::RandGen randgen;
        BICYCL::SecLevel seclevel = BICYCL::SecLevel::_128;
        // Group order for BLS12-381 scalar field (Fr)
        BICYCL::Mpz q("0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001");
        BICYCL::CL_HSMqk C(q, 1, seclevel, randgen);
        size_t party_count = 20;
        size_t threshold = 20;
        BICYCL::Mpz delta = factorial(party_count);
        BICYCL::Mpz coefficient_bound;
    const size_t bit_bound = C.secretkey_bound().nbits() - 124;
    BICYCL::Mpz::mulby2k(coefficient_bound, BICYCL::Mpz("1"), bit_bound);

    std::vector<BICYCL::CL_HSMqk::SecretKey> class_group_secret_shares;
    std::vector<BICYCL::Mpz> class_group_secret_share_values;
    std::vector<BICYCL::CL_HSMqk::PublicKey> class_group_public_shares;
    class_group_secret_shares.reserve(party_count);
    class_group_secret_share_values.reserve(party_count);
    class_group_public_shares.reserve(party_count);

    BICYCL::Mpz alpha(randgen.random_mpz(coefficient_bound));
    BICYCL::Mpz class_group_u;
    BICYCL::Mpz class_group_secret_key;
    BICYCL::Mpz::mul(class_group_u, alpha, delta);
    BICYCL::Mpz::mul(class_group_secret_key, class_group_u, delta);

    BICYCL::CL_HSMqk::SecretKey class_group_secret_key_delta(C, class_group_secret_key);
    BICYCL::CL_HSMqk::PublicKey class_group_public_key = C.keygen(class_group_secret_key_delta);

    std::vector<BICYCL::Mpz> class_group_coefficients;
    class_group_coefficients.reserve(threshold);
    for (size_t coefficient_index = 0; coefficient_index < threshold; ++coefficient_index) {
        class_group_coefficients.emplace_back(randgen.random_mpz(coefficient_bound));
    }

    for (size_t party_index = 0; party_index < party_count; ++party_index) {
        BICYCL::Mpz share_value = class_group_coefficients[threshold - 1];
        for (size_t coefficient_index = threshold - 1; coefficient_index > 0; --coefficient_index) {
            BICYCL::Mpz::mul(share_value, share_value, BICYCL::Mpz(party_index + 1));
            BICYCL::Mpz::add(share_value, share_value, class_group_coefficients[coefficient_index - 1]);
        }
        BICYCL::Mpz::mul(share_value, share_value, BICYCL::Mpz(party_index + 1));
        BICYCL::Mpz::add(share_value, share_value, class_group_u);
        class_group_secret_share_values.push_back(share_value);
        class_group_secret_shares.emplace_back(C, share_value);
        class_group_public_shares.emplace_back(C, class_group_secret_shares.back());
    }


    
    
        // Pick in the interval of [0, bound)
        BICYCL::Mpz bound(C.encrypt_randomness_bound());
        BICYCL::Mpz alpha_rho(randgen.random_mpz(bound));
    
        BICYCL::Mpz pt_bound = C.cleartext_bound();
        BICYCL::Mpz random_val = randgen.random_mpz(pt_bound);
        BICYCL::CL_HSMqk::CipherText ciphertext(C,class_group_public_key,BICYCL::CL_HSMqk::ClearText(C,random_val),  alpha_rho);
        std::vector<BICYCL::QFI> partial_decryption(party_count);
        
       BICYCL::HashAlgo hash_algorithm(seclevel);
       bool b;
    
        
         for(size_t i = 0; i<party_count - 1; i++)
         {
            partial_decrypt(C,class_group_secret_shares[i], ciphertext,partial_decryption[i]);
         }
    {
          Timer t_total("Threshold decryption");
         partial_decrypt(C,class_group_secret_shares[party_count-1], ciphertext,partial_decryption[party_count-1]);
         BICYCL::CL_HSMqk_Part_Dec_ZKProof proof(C,
                             hash_algorithm,
                             class_group_public_shares[party_count-1],
                             ciphertext,
                             partial_decryption[party_count-1],
                             class_group_secret_shares[party_count-1],
                             randgen);
       for(size_t i = 0; i< party_count; i++)
       {
            
            b = proof.verify(C,
                hash_algorithm,
                class_group_public_shares[party_count-1],
                ciphertext,
                partial_decryption[party_count-1]);
            if(!b)
                std::cout<<"Verify failed"<<std::endl;
       }
       auto result = aggregate_partial_ciphertext(C,threshold,partial_decryption,ciphertext,delta);
       
        

}

    return 0;
}