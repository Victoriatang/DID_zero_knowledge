#include <set>
#include <iostream>
#include <vector>
#include "BBS_signature/prove_CL.h"
#include <chrono> // Added for timing


int main() {
    try {
        std::cout << "=== Starting BBS+ with CL-Proof Test ===" << std::endl;
        
        // 1. Setup BBS Parameters and Keys
        size_t max_msg_count = 2;
        BBSParams params(mcl::BLS12_381, max_msg_count);
        BBS engine(params);
        
        BBS::PublicKey pk;
        BBS::SecretKey sk;
        size_t current_msg_count = 2; 
        engine.generate_keys(pk, sk);
        std::cout << "[+] Key pair generated successfully." << std::endl;

        // 2. Prepare Messages
        std::vector<mcl::bn::Fr> messages(current_msg_count);
        messages[0].setStr("123456789");
        messages[1].setStr("987654321");

        // 3. Basic BBS Signing and Verification
        auto sig = engine.sign(messages, pk, sk);
        bool is_valid = engine.verify(messages, sig, pk);

        if (is_valid) {
            std::cout << "[+] Basic BBS signature is VALID." << std::endl;
        } else {
            std::cerr << "[-] Basic BBS signature is INVALID." << std::endl;
            return 1;
        }

        // 4. Setup BICYCL CL-HSM Parameters
        BICYCL::RandGen randgen;
        BICYCL::SecLevel seclevel = BICYCL::SecLevel::_128;
        // Group order for BLS12-381 scalar field (Fr)
        BICYCL::Mpz q("0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001");
        BICYCL::CL_HSMqk C(q, 1, seclevel, randgen);

        BICYCL::CL_HSMqk::SecretKey cl_sk = C.keygen(randgen);
        BICYCL::CL_HSMqk::PublicKey cl_pk = C.keygen(cl_sk);
        
        IssueParams pp(params, pk, C, randgen, cl_pk);
        mcl_WT wt;
        mcl_ST st;
        BICYCL::Mpz rho;

        // --- Start Timing: Proof Generation ---
        std::cout << "[*] Generating CL-BBS Proof..." << std::endl;
        auto start_proof = std::chrono::high_resolution_clock::now();
        // 5. Generate Issue Material (Blinded signature/Commitments)
        BICYCL::CL_HSMqk::CipherText C_sid = GenIssueMaterial(pp, messages, sig, st, wt, rho);
        
        mcl_COM com;
        mcl_RESP resp;

        
        
        
        auto C_prime = CL_BBS_proof(pp, st, C_sid, wt, rho, com, resp);
        
        auto end_proof = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> proof_duration = end_proof - start_proof;
        // --- End Timing: Proof Generation ---

        // --- Start Timing: Verification ---
        std::cout << "[*] Verifying CL-BBS Proof..." << std::endl;
        auto start_verify = std::chrono::high_resolution_clock::now();
        
        bool proof_valid = CL_BBS_verify(pp, st, C_sid, com, C_prime, resp);
        
        auto end_verify = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> verify_duration = end_verify - start_verify;
        // --- End Timing: Verification ---

        // 6. Output Statistics
        std::cout << "------------------------------------------" << std::endl;
        if (proof_valid) {
            std::cout << ">>> FINAL RESULT: PROOF VALID <<<" << std::endl;
        } else {
            std::cerr << ">>> FINAL RESULT: PROOF INVALID <<<" << std::endl;
        }
        
        std::cout << "Proof Generation Time : " << proof_duration.count() << " ms" << std::endl;
        std::cout << "Proof Verification Time: " << verify_duration.count() << " ms" << std::endl;
        std::cout << "------------------------------------------" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "CRITICAL ERROR: " << e.what() << std::endl;
        return 1;
    }
}