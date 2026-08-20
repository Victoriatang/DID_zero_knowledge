#include "prove.h"
#include <stdexcept>
#include <string>
#include <utility>
#include <random>
#include <fstream>

// ── G1 Serialization ────────────────────
static constexpr size_t G1_BYTES = 48;
static void g1ToBytes(const G1& p, uint8_t out[G1_BYTES]) {
    size_t n = p.serialize(out, G1_BYTES);
    assert(n == G1_BYTES);
}


size_t ceil_log2(size_t N) {
    size_t k = 0;
    size_t v = 1;
    while (v < N) {
        v <<= 1;
        k++;
    }
    return k;
}

size_t next_pow2(size_t x) {
    size_t v = 1;
    while (v < x) {
        v <<= 1;
    }
    return v;
}
// File-based binary search (keep file open to avoid repeated open/close overhead)
int binarySearchBinary(const std::string& filename, size_t N, const G1& target) {
    uint8_t target_bytes[G1_BYTES];
    g1ToBytes(target, target_bytes);

    std::ifstream ifs(filename, std::ios::binary);
    int left = 0, right = static_cast<int>(N) - 1;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        ifs.seekg(static_cast<std::streamoff>(mid * G1_BYTES));
        uint8_t buf[G1_BYTES];
        ifs.read(reinterpret_cast<char*>(buf), G1_BYTES);

        int cmp = std::memcmp(buf, target_bytes, G1_BYTES);
        if (cmp == 0)      return mid;
        else if (cmp < 0)  left = mid + 1;
        else               right = mid - 1;
    }
    return -1;
}

ACProof::ACProof(publicParams& pp) : pp_(pp) {}
ACProof::~ACProof() = default;

ACProof::proof ACProof::Prove(const ACProof::statement& st, const ACProof::witness& wt){
    // 1. (a) Pick a random blinding factor r from Z*_q
    mcl::Fr r,r_bf;
    r.setByCSPRNG();
    r_bf.setByCSPRNG();

    // 2. (b) Compute blinded values
    com cm;
    
    // Compute A_bar = r * A
    mcl::G1::mul(cm.A_bar, wt.sig.A, r);
    // Compute B_bar = r*G1 + (r*sid)*H1 + (r*m_r)*H2 + (-e)*A_bar
    mcl::G1 term1, term2, term3, term4;
    mcl::G1::mul(term1, pp_.params.g1,r);           // term: r * G1
    
    mcl::Fr r_sid;
    mcl::Fr::mul(r_sid, r, wt.sid);
    mcl::G1::mul(term2, pp_.pk.pks_g1[0], r_sid);       // term: (r * sid) * H1
    
    mcl::Fr r_mr;
    mcl::Fr::mul(r_mr, r, wt.m_r);
    mcl::G1::mul(term3, pp_.pk.pks_g1[1], r_mr);        // term: (r * m_r) * H2
    
 
    mcl::G1::mul(term4, cm.A_bar, wt.sig.e);   // term: e * A_bar

    
    // B_bar = term1 + term2 + term3 + term4
    mcl::G1::add(cm.B_bar, term1, term2);
    mcl::G1::add(cm.B_bar, cm.B_bar, term3);
    mcl::G1::sub(cm.B_bar, cm.B_bar, term4);

 
    // Define vk_bar = r*vk + r_bf*L3
    mcl::G2 temp;
    mcl::G2::mul(cm.vk_bar,st.vk,r);
    mcl::G2::mul(temp,pp_.l3,r_bf);
    mcl::G2::add(cm.vk_bar,cm.vk_bar,temp);


    
    // 3. (c) Pick random commitment scalars (alpha values)
    mcl::Fr alpha_r, alpha_e, alpha_v, alpha_sid, alpha_Q, alpha_bf, alpha_mr, alpha_ctr, alpha_rT;
    alpha_r.setByCSPRNG();
    alpha_e.setByCSPRNG();
    alpha_v.setByCSPRNG();
    alpha_sid.setByCSPRNG();
    alpha_Q.setByCSPRNG();
    alpha_bf.setByCSPRNG();
    alpha_mr.setByCSPRNG();
    alpha_ctr.setByCSPRNG();
    alpha_rT.setByCSPRNG();

    // 4. (d) Compute Commitment vector 
    
    // U_1 = alpha_r * vk + alpha_bf*L3
    mcl::G2::mul(cm.U_1, st.vk, alpha_r);
    mcl::G2::mul(temp,pp_.l3,alpha_bf);
    mcl::G2::add(cm.U_1,cm.U_1,temp);

    // U_2 = alpha_sid * L1 + alpha_v * L2 + alpha_bf * L3
    mcl::G2 u2_t1, u2_t2, u2_t3;
    mcl::G2::mul(u2_t1, pp_.l1, alpha_sid);
    mcl::G2::mul(u2_t2, pp_.l2, alpha_v);
    mcl::G2::mul(u2_t3,pp_.l3,alpha_bf);
    
    mcl::G2::add(cm.U_2, u2_t1, u2_t2);
    mcl::G2::add(cm.U_2, cm.U_2, u2_t3);

   
    // U_B = alpha_r * G1 + alpha_e * A_bar + alpha_sid * H1 + alpha_mr * H2
    mcl::G1 ub_t1, ub_t2, ub_t3, ub_t4;
    mcl::G1::mul(ub_t1, pp_.params.g1, alpha_r);
    mcl::G1::mul(ub_t2, cm.A_bar, alpha_e);
    mcl::G1::mul(ub_t3, pp_.pk.pks_g1[0], alpha_sid);
    mcl::G1::mul(ub_t4, pp_.pk.pks_g1[1], alpha_mr);
    
    mcl::G1::add(cm.U_B, ub_t1, ub_t2);
    mcl::G1::add(cm.U_B, cm.U_B, ub_t3);
    mcl::G1::add(cm.U_B, cm.U_B, ub_t4);

    // U_Q = e(tag, L2)^alpha_Q
    mcl::GT e_tag_l2;
    mcl::pairing(e_tag_l2, st.tag, pp_.l2);
    mcl::GT::pow(cm.U_Q, e_tag_l2, alpha_Q);

    mcl::G2::mul(u2_t1,pp_.l1,alpha_ctr);
    mcl::G2::mul(cm.G2_ctr,pp_.l2,alpha_rT);
    mcl::G2::add(cm.G2_ctr,cm.G2_ctr,u2_t1);

    mcl::G1::mul(ub_t1,pp_.params.g1,alpha_ctr);
    mcl::G1::mul(cm.G1_ctr,pp_.pk.pks_g1[0],alpha_rT);
    mcl::G1::add(cm.G1_ctr,cm.G1_ctr,ub_t1);


    //For range proof
    size_t n = next_pow2(ceil_log2(pp_.N));
    std::vector<mcl::G1> g_vec(2*n);
    std::vector<mcl::G1> h_vec(2*n);
    for (size_t i = 0; i<2*n; i++)
    {
        std::string seed = "range proof for g" + std::to_string(i);
        mcl::hashAndMapToG1(g_vec[i], seed);
        seed = "range proof for h" + std::to_string(i);
        mcl::hashAndMapToG1(h_vec[i], seed);
    }
    RangeProof proof;
    std::vector<mcl::Fr> v_vec(2), r_vec(2);
    v_vec[0] = wt.ctr;
    mcl::Fr::sub(v_vec[1],pp_.N,wt.ctr);
    r_vec[0] = wt.r_T;
    mcl::Fr::neg(r_vec[1],wt.r_T);
    proof = prove_range(g_vec, h_vec, pp_.params.g1, pp_.pk.pks_g1[0],v_vec,r_vec,n);


    // 5. Generate challenge k
    mcl::Fr k = k_from_hash(st, cm);

    
    // 6. Compute responses (beta values)
    resp rp;
    // beta_r = alpha_r + r * k
    mcl::Fr r_k;
    mcl::Fr::mul(r_k, r, k);
    mcl::Fr::add(rp.beta_r, alpha_r, r_k);
    //beta_bf = alpha_bf+r_bf * k
    mcl::Fr::mul(r_k,r_bf,k);
    mcl::Fr::add(rp.beta_bf,r_k,alpha_bf);

     // beta_sid = alpha_sid + (r * sid) * k
    mcl::Fr r_sid_k;
    mcl::Fr::mul(r_sid_k, r_sid, k); // r_sid was calculated in Step 2b
    mcl::Fr::add(rp.beta_sid, alpha_sid, r_sid_k);
    // beta_vk = alpha_vk + (r * r_vk) * k
    mcl::Fr r_mr_k;
    mcl::Fr::mul(r_mr_k, r, wt.r_vk); 
    mcl::Fr::mul(r_mr_k, r_mr_k, k);
    mcl::Fr::add(rp.beta_v, alpha_v, r_mr_k);
    // beta_e = alpha_e - e * k
    mcl::Fr e_k;
    mcl::Fr::mul(e_k, wt.sig.e, k);
    mcl::Fr::sub(rp.beta_e, alpha_e, e_k);
    // beta_mr = alpha_mr + r * m_r * ch
    mcl::Fr::mul(r_mr_k,r,wt.m_r);
    mcl::Fr::mul(r_mr_k,r_mr_k,k);
    mcl::Fr::add(rp.beta_mr,alpha_mr,r_mr_k);

    // beta_Q = alpha_Q + (r_vk + r_T) * k
    mcl::Fr sum_rvk_rt, r_sum_k;
    mcl::Fr::add(sum_rvk_rt, wt.r_vk, wt.r_T);
    mcl::Fr::mul(r_sum_k, sum_rvk_rt, k);
    mcl::Fr::add(rp.beta_Q, alpha_Q, r_sum_k);

    //beta_ctr = alpha_ctr + k * ctr, beta_rT = alpha_rT + k * r_T
    mcl::Fr::mul(r_k,wt.ctr,k);
    mcl::Fr::add(rp.beta_ctr,alpha_ctr,r_k);
    mcl::Fr::mul(r_k,wt.r_T,k);
    mcl::Fr::add(rp.beta_rT,alpha_rT,r_k);
   
    return {cm, rp, proof};
}
bool ACProof::Verify(const ACProof::statement& st, const ACProof::proof& pf, const std::string& filename, size_t N) {
    // 1. Initial check
    if (pf.cm.A_bar.isZero()) {
        std::cout << "[Verify] Failed: A_bar is zero." << std::endl;
        return false;
    }
  
   
    // 2. Recompute challenge k
    mcl::Fr k = k_from_hash(st, pf.cm);

    // Equation 1: e(A_bar, X) = e(B_bar, G2)
    mcl::GT e_A_X, e_B_G2;
    mcl::pairing(e_A_X, pf.cm.A_bar, pp_.pk.pk_g2); 

    mcl::pairing(e_B_G2, pf.cm.B_bar, pp_.params.g2);
    if (e_A_X != e_B_G2) {
        std::cout << "[Verify] Failed at Equation 1 (Pairing check A_bar/B_bar)" << std::endl;
        return false;
    }
    // 5. Equation 2: U_vk1 + ch * vk_bar = beta_r * vk + beta_bf * L3
    mcl::G2 beta_L3, k_vk_bar, lhs2, rhs2;
    mcl::G2::mul(beta_L3, pp_.l3, pf.rp.beta_bf); 
    mcl::G2::mul(rhs2, st.vk, pf.rp.beta_r);
    mcl::G2::add(rhs2, rhs2, beta_L3);
    mcl::G2::mul(k_vk_bar, pf.cm.vk_bar, k); 
    mcl::G2::add(lhs2, pf.cm.U_1, k_vk_bar);
    if (lhs2 != rhs2) {
        std::cout << "[Verify] Failed at Equation 2 (U_1 check)" << std::endl;
        return false;
    }

    // Equation 3: U_vk2 + ch * vk_bar = beta_sid * L1 + beta_v * L2 + beta_bf * L3
    mcl::G2 lhs3, rhs3, term1, term2;
    mcl::G2::add(lhs3, pf.cm.U_2, k_vk_bar);
    mcl::G2::mul(term1, pp_.l1, pf.rp.beta_sid);
    mcl::G2::mul(term2, pp_.l2, pf.rp.beta_v);
    mcl::G2::mul(rhs3, pp_.l3, pf.rp.beta_bf);
    mcl::G2::add(rhs3, rhs3, term2);
    mcl::G2::add(rhs3, rhs3, term1);
    if (lhs3 != rhs3) {
        std::cout << "[Verify] Failed at Equation 3 (U_2 check)" << std::endl;
        return false;
    }

    // 6. Equation 4: U_B + ch * B_bar = beta_r * G1 + beta_e * A_bar + beta_sid * H1 + beta_mr * H2
    mcl::G1 lhs4, rhs4, ch_B, t1, t2, t3, t4;
    mcl::G1::mul(ch_B, pf.cm.B_bar, k);
    mcl::G1::add(lhs4, pf.cm.U_B, ch_B);

    mcl::G1::mul(t1, pp_.params.g1, pf.rp.beta_r);
    mcl::G1::mul(t2, pf.cm.A_bar, pf.rp.beta_e);
    mcl::G1::mul(t3, pp_.pk.pks_g1[0], pf.rp.beta_sid);
    mcl::G1::mul(t4, pp_.pk.pks_g1[1], pf.rp.beta_mr);
    
    mcl::G1::add(rhs4, t1, t2);
    mcl::G1::add(rhs4, rhs4, t3);
    mcl::G1::add(rhs4, rhs4, t4);
    if (lhs4 != rhs4) {
        std::cout << "[Verify] Failed at Equation 4 (U_B check)" << std::endl;
        return false;
    }

    // 7. Equation 5: U_Q + ch * Q = beta_Q * e(tag, L2)
    // 3. Compute Q = e(tag, vk + T) / e(G_epoch, L1)
    mcl::G2 vk_T;
    mcl::G2::add(vk_T, st.vk, st.T);
     mcl::GT pair_lhs, pair_rhs, Q;
    mcl::pairing(pair_lhs, st.tag, vk_T);
    mcl::pairing(pair_rhs, pp_.gepoch, pp_.l1);
    mcl::GT::div(Q, pair_lhs, pair_rhs); 
    mcl::GT lhs5, rhs5, Q_k, e_t_l2;
    mcl::GT::pow(Q_k, Q, k);
    mcl::GT::mul(lhs5, pf.cm.U_Q, Q_k);

    mcl::pairing(e_t_l2, st.tag, pp_.l2);
    mcl::GT::pow(rhs5, e_t_l2, pf.rp.beta_Q);
    if (lhs5 != rhs5) {
        std::cout << "[Verify] Failed at Equation 5 (U_Q check)" << std::endl;
        return false;
    }
    mcl::G2 rhs6, lhs6;
    mcl::G2::mul(vk_T,pp_.l1,pf.rp.beta_ctr);
    mcl::G2::mul(rhs6,pp_.l2,pf.rp.beta_rT);
    mcl::G2::add(rhs6,rhs6,vk_T);
    mcl::G2::mul(lhs6,st.T,k);
    mcl::G2::add(lhs6,lhs6,pf.cm.G2_ctr);
    if (lhs6 != rhs6)
    {
        std::cout << "[Verify] Failed at Equation 6" << std::endl;
        return false;

    }
    mcl::G1 rhs7, lhs7;
    mcl::G1::mul(t1,pp_.params.g1,pf.rp.beta_ctr);
    mcl::G1::mul(rhs7,pp_.pk.pks_g1[0],pf.rp.beta_rT);
    mcl::G1::add(rhs7,rhs7,t1);
    mcl::G1::mul(lhs7,st.T_prime,k);
    mcl::G1::add(lhs7,lhs7,pf.cm.G1_ctr);
    if (lhs7 != rhs7)
    {
        std::cout << "[Verify] Failed at Equation 7" << std::endl;
        return false;

    }


    size_t n = next_pow2(ceil_log2(pp_.N));
    
    std::vector<mcl::G1> g_vec(2*n);
    std::vector<mcl::G1> h_vec(2*n);
    for (size_t i = 0; i<2*n; i++)
    {
        std::string seed = "range proof for g" + std::to_string(i);
        mcl::hashAndMapToG1(g_vec[i], seed);
        seed = "range proof for h" + std::to_string(i);
        mcl::hashAndMapToG1(h_vec[i], seed);
    }
    std::vector<mcl::G1> ped_com(2);
    ped_com[0] = st.T_prime;
    mcl::G1::mul(ped_com[1],pp_.params.g1,pp_.N);
    mcl::G1::sub(ped_com[1],ped_com[1],ped_com[0]);
    bool ok =
    fast_verify_range(pf.proof,g_vec,h_vec,pp_.params.g1, pp_.pk.pks_g1[0],ped_com,n);
    if(!ok){
        std::cout << "[Verify] Failed at range proof" << std::endl;
        return false;
    }

    auto idx = binarySearchBinary(filename, N, st.tag);
    if (idx > -1){
        std::cout << "[Verify] Failed -- The tag has been revoked" << std::endl;
        return false;
    }
    


    std::cout << "[Verify] All checks passed!" << std::endl;
    return true;
}
 mcl::Fr ACProof::k_from_hash(const ACProof::statement& st, const ACProof::com &cm) {
        std::ostringstream output;
        output << st.T << st.tag << st.vk << cm.A_bar << cm.B_bar << cm.U_1 << cm.U_2 << cm.U_B << cm.U_Q << cm.vk_bar
        << cm.G1_ctr << cm.G2_ctr;
        mcl::Fr challenge;
        challenge.setHashOf(output.str());
        return challenge;
        
    }

void ACProof::GenMatrials(const std::vector<mcl::Fr>& messages,
    const BBS::Signature &sig,
    ACProof::statement &st, ACProof::witness &wt){
        wt.sid = messages[0];
        wt.m_r = messages[1];
        wt.sig = sig;

        std::mt19937 rng(0);  
        uint64_t val = rng() % pp_.N;
        wt.ctr = (int64_t)val;  
        wt.r_T.setByCSPRNG();
        wt.r_vk.setByCSPRNG();
        mcl::Fr r1;
        mcl::Fr::add(r1,wt.sid,wt.ctr);
        mcl::Fr::inv(r1,r1);
        mcl::G1::mul(st.tag,pp_.gepoch,r1);
        mcl::G2 R2;
        mcl::G2::mul(R2,pp_.l2,wt.r_T);
        mcl::G2::mul(st.T,pp_.l1,wt.ctr);
        mcl::G2::add(st.T,st.T,R2);
        mcl::G2::mul(R2,pp_.l2,wt.r_vk);
        mcl::G2::mul(st.vk,pp_.l1,wt.sid);
        mcl::G2::add(st.vk,st.vk,R2); 
        mcl::G1 R2_prime;
        mcl::G1::mul(R2_prime,pp_.pk.pks_g1[0],wt.r_T);
        mcl::G1::mul(st.T_prime,pp_.params.g1,wt.ctr);
        mcl::G1::add(st.T_prime,st.T_prime,R2_prime);
    }