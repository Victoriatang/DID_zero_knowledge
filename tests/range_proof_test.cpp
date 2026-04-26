#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <chrono>
#include <mcl/bls12_381.hpp>
#include <mcl/gmp_util.hpp>
#include <random>
#include "BBS_signature/range_proof.h"

// 计时工具
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


// 辅助函数：根据索引生成确定的伪随机点 (对应 Rust 的 generate_random_point)
G1 generate_test_point(size_t index, const std::string& label) {
    std::string seed = label + std::to_string(index);
    G1 P;
    mcl::hashAndMapToG1(P, seed);
    return P;
}
// 用标签字符串确定性地生成一个 G1 点
static G1 hashToG1(const std::string& label) {
    G1 out;
    mcl::hashAndMapToG1(out, label.data(), label.size());
    return out;
}
static G1 hashToG1WithIndex(uint32_t index, uint32_t seed) {
    // label = seed + index，转成字符串作为哈希输入
    std::string label = std::to_string((uint64_t)seed + index);
    return hashToG1(label);
}
void test_inner_product_argument(size_t n) {
    std::cout << "Running IPA test for n = " << n << "..." << std::endl;

    // 1. 初始化基点向量 g_vec 和 h_vec
    std::vector<G1> g_vec(n);
    std::vector<G1> h_vec(n);
    std::string label = "innerproduct";

    for (size_t i = 0; i < n; ++i) {
        g_vec[i] = generate_test_point(i, label);
        h_vec[i] = generate_test_point(i + n, label); // 偏移 n 以确保不同
    }

    // 2. 生成辅助基点 Gx (ux)
    G1 Gx = generate_test_point(999, "Gx_label");

    // 3. 生成随机秘密向量 a 和 b
    std::vector<Fr> a(n), b(n);
    for (size_t i = 0; i < n; ++i) {
        a[i].setByCSPRNG();
        b[i].setByCSPRNG();
    }

    // 4. 计算内积 c = <a, b>
    Fr c = inner_product(a, b);

    // 5. 模拟 hi_tag 计算 (对应 Rust 的 y 盲化逻辑)
    std::vector<G1> hi_tag = h_vec; 
    
    Fr y; y.setByCSPRNG();
    Fr y_inv; Fr::inv(y_inv, y);
    Fr y_inv_pow(1);
    for(size_t i=0; i<n; ++i) {
        G1::mul(hi_tag[i], h_vec[i], y_inv_pow);
        Fr::mul(y_inv_pow, y_inv_pow, y_inv);
    }
    

    // 6. 构造承诺 P = MSM(g_vec, a) + MSM(hi_tag, b) + Gx * c
    G1 P, tmp;
    G1::mulVec(P, g_vec.data(), a.data(), n);   // P = sum(a_i * g_i)
    G1::mulVec(tmp, hi_tag.data(), b.data(), n); // tmp = sum(b_i * h_i)
    G1::add(P, P, tmp);
    
    G1 ux_c;
    G1::mul(ux_c, Gx, c);
    G1::add(P, P, ux_c); // P = P + Gx * c

    // 7. 执行证明 (Prover)
    std::vector<G1> L_vec, R_vec;
    InnerProductArg ipp = prove(g_vec, hi_tag, Gx, P, a, b, L_vec, R_vec);

    // 8. 执行验证 (Verifier)
    bool is_valid = fast_verify(ipp, g_vec, hi_tag, Gx, P);

    // 9. 断言结果
    if (is_valid) {
        std::cout << "Test PASSED for n = " << n << std::endl;
    } else {
        std::cerr << "Test FAILED for n = " << n << std::endl;
        assert(false);
    }
}

void test_helper_range_proof(uint32_t seed, size_t n, size_t m) {
    const size_t nm = n * m;

    // G = 固定生成元（对应 Rust 的 Point::generator()）
    G1 G;
    mcl::hashAndMapToG1(G, "generator_G", 11);

    // H = hash("1") 映射到的点（对应 Rust 的 label=1）
    G1 H = hashToG1("1");

    // g_vec[i] 对应 Rust 的 hash(i + seed)
    std::vector<G1> g_vec(nm);
    for (size_t i = 0; i < nm; i++) {
        g_vec[i] = hashToG1WithIndex(static_cast<uint32_t>(i), seed);
    }

    // h_vec[i] 对应 Rust 的 hash(n + i + seed)
    std::vector<G1> h_vec(nm);
    for (size_t i = 0; i < nm; i++) {
        h_vec[i] = hashToG1WithIndex(static_cast<uint32_t>(n + i), seed);
    }

    // range = 2^n，v_vec[i] ∈ [0, 2^n)
    // 用 C++ 随机数生成，模拟 BigInt::sample_below(&range)
    std::mt19937 rng(seed);  // 用 seed 初始化，保证可复现
    uint64_t range = (uint64_t)1 << n;  // 假设 n <= 63
    auto sampleBelow = [&]() -> Fr {
        uint64_t val = rng() % range;
        Fr r;
        r = (int64_t)val;  // mcl Fr 支持从整数赋值
        return r;
    };

    std::vector<Fr> v_vec(m), r_vec(m);
    for (size_t i = 0; i < m; i++) {
        v_vec[i] = sampleBelow();
        r_vec[i].setByCSPRNG();
    }
    // ped_com[i] = G * v_vec[i] + H * r_vec[i]
    std::vector<G1> ped_com(m);
    for (size_t i = 0; i < m; i++) {
        G1 Gv, Hr;
        G1::mul(Gv, G, v_vec[i]);
        G1::mul(Hr, H, r_vec[i]);
        G1::add(ped_com[i], Gv, Hr);
    }
    RangeProof proof;

    printf("\n=== Range Proof Proving ===\n");
    {
        Timer t_total("Total");
        proof = prove_range(g_vec, h_vec, G, H, v_vec, r_vec, n);

    }
    bool ok;
    printf("\n=== Range Proof Verifying ===\n");
    {
        Timer t_total("Total");
        ok = fast_verify_range(proof, g_vec, h_vec, G, H, ped_com, n);

    }
    
    // 9. 断言结果
    if (ok) {
        std::cout << "Range Proof Test PASSED for n = " << n << std::endl;
    } else {
        std::cerr << "Range Proof Test FAILED for n = " << n << std::endl;
        assert(false);
    }
}

int main() {
    initPairing(mcl::BLS12_381);
   /* printf("=== Inner Product Argument (n=16) ===\n");
    {
        Timer t_total("Total");
        test_inner_product_argument(16);
    }
        */
    printf("\n=== Range Proof ( n=8, m=2) ===\n");
    test_helper_range_proof(0, 8, 2);
    
    printf("\n=== Range Proof (n=16, m=2) ===\n");
    test_helper_range_proof(0, 16, 2);
    

    return 0;
}