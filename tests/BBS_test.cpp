#include <set>
#include <iostream>
#include <vector>
#include "BBS_signature/sign.h"
#include "BBS_signature/prove.h"

#include <chrono> // Added for timing

int main() {
try {
std::cout << "=== Test ===" << std::endl;

    size_t max_msg_count = 2;
    BBSParams params(mcl::BLS12_381, max_msg_count);

    BBS engine(params);
 

    BBS::PublicKey pk;
    BBS::SecretKey sk;
    size_t current_msg_count = 2;
    engine.generate_keys(pk, sk);

    std::vector<mcl::bn::Fr> messages(current_msg_count);
    messages[0].setStr("123456789");
    messages[1].setStr("987654321");
    

    auto sig = engine.sign(messages, pk, sk);
  
    bool is_valid = engine.verify(messages, sig, pk);
   

    if (!is_valid) {
        std::cerr << ">>> BBS Result: INVALID <<<" << std::endl;
        return 1;
    }

    const size_t N = 1000000;
    const std::string filename = "g1_arr_1000000.bin";
    std::vector<size_t> Ns = {200, 500, 800, 1000};


    for (size_t i = 0; i < Ns.size(); i++) {

        // Initialize ACProof components
        ACProof::statement st;
        ACProof::witness wt;
        std::cout << "\n--- Testing for usage limit N=" << Ns[i] << " ---" << std::endl;
     
        publicParams pp(params, pk, Ns[i]);
      

        ACProof ac_engine(pp);
      

        // --- Measure Prove Algorithm ---
        auto start_prove = std::chrono::steady_clock::now();

  
        ac_engine.GenMatrials(messages, sig, st, wt);
   


        ACProof::proof pf = ac_engine.Prove(st, wt);
  

        auto end_prove = std::chrono::steady_clock::now();
        auto prove_duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(end_prove - start_prove).count();
        std::cout << "Prove Execution Time: " << prove_duration << " ms" << std::endl;

        // --- Measure Verify Algorithm ---
        auto start_verify = std::chrono::steady_clock::now();

      
        bool ac_valid = ac_engine.Verify(st, pf, filename, N);
        
        auto end_verify = std::chrono::steady_clock::now();
        auto verify_duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(end_verify - start_verify).count();
        std::cout << "Verify Execution Time: " << verify_duration << " ms" << std::endl;

        if (!ac_valid) {
            std::cerr << ">>> ACProof Result: INVALID <<<" << std::endl;
        } else {
            std::cout << ">>> ACProof Result: VALID <<<" << std::endl;
        }
    }

   
    return 0;

} catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
}
}