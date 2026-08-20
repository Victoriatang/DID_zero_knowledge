#ifndef PROVECL_H
#define PROVECL_H

#include "sign.h"
#include "bicycl.hpp"

struct mcl_ST {
    mcl::G1 null, B_sid;
};

struct mcl_WT {
    mcl::Fr       usk, r_ob, sid, m_r;
    BBS::Signature sig;
};

struct mcl_COM {
    mcl::G1 A_bar, B_bar, N_bar;
    mcl::G1 U_N1, U_N2, U_barB, U_Bprime;
};

inline
static void Fr_to_Mpz(BICYCL::Mpz& out, const mcl::Fr& fr) {
    out = fr.getStr();
}

inline
static mcl::Fr computeChallenge(
    const mcl::G1& null_pt,
    const mcl::G1& B_sid,
    const BICYCL::CL_HSMqk::CipherText& C_sid,
    const mcl_COM& com,
    const BICYCL::CL_HSMqk::CipherText& C_prime)
{
    std::ostringstream oss;
    // statement
    oss << null_pt << B_sid;

    oss << C_sid.c1() << C_sid.c2();
    // commitment values
    oss << com.A_bar << com.B_bar << com.N_bar
        << com.U_N1  << com.U_N2
        << com.U_barB << com.U_Bprime;
    // C'
    oss << C_prime.c1() << C_prime.c2();

    mcl::Fr ch;
    const std::string s = oss.str();
    ch.setHashOf(s.data(), s.size());
    return ch;
}

struct IssueParams {
    const BBSParams&               params_;
    const BBS::PublicKey&          pk_;
    mcl::G1                        H_s;
    BICYCL::CL_HSMqk&              C_;
    BICYCL::RandGen&               randgen_;
    BICYCL::CL_HSMqk::PublicKey    cl_pk_;

    IssueParams(const BBSParams&                    params,
                const BBS::PublicKey&               pk,
                BICYCL::CL_HSMqk&                  C,
                BICYCL::RandGen&                   randgen,
                const BICYCL::CL_HSMqk::PublicKey& cl_pk)
        : params_(params), pk_(pk), C_(C), randgen_(randgen), cl_pk_(cl_pk)
    {
        mcl::hashAndMapToG1(H_s, "H_s");
    }
};

struct mcl_RESP {
    mcl::Fr beta_r, beta_usk, beta_ob, beta_e, beta_sid, beta_mr;
    BICYCL::Mpz beta_rho;
};

inline
BICYCL::CL_HSMqk::CipherText GenIssueMaterial(
    const IssueParams&    pp,
    const std::vector<mcl::Fr>& messages,
    const BBS::Signature&       sig,
    mcl_ST&               st,
    mcl_WT&               wt,
    BICYCL::Mpz&          rho)
{
    const BBSParams&                    bbs   = pp.params_;
    const BBS::PublicKey&               bpk   = pp.pk_;
    BICYCL::CL_HSMqk&                   CL    = pp.C_;
    BICYCL::RandGen&                    rng   = pp.randgen_;
    const BICYCL::CL_HSMqk::PublicKey&  cl_pk = pp.cl_pk_;
    const mcl::G1& G1  = bbs.g1;
    const mcl::G1& H1  = bpk.pks_g1[0];
    const mcl::G1& H2  = bpk.pks_g1[1];
    const mcl::G1& H_s = pp.H_s;

    // ----------------------------------------------------------
    // Get usk and r_ob from messages
    // ----------------------------------------------------------
    assert(messages.size() >= 2 && "messages must contain at least usk and r_ob");
    wt.usk  = messages[0];
    wt.r_ob = messages[1];
    wt.sig  = sig;

    // ----------------------------------------------------------
    // Step 1. Compute the nullifier
    //   null <- usk * H_s in G1
    // ----------------------------------------------------------
    mcl::G1::mul(st.null, H_s, wt.usk);

    // ----------------------------------------------------------
    // Step 2a. Pick sid <-$ Z_q
    //          Encrypt: C_sid = (g_q^rho, f^sid * pk_cl^rho)
    // ----------------------------------------------------------
    wt.sid.setByCSPRNG();

    // Convert sid (mcl::Fr) to BICYCL::Mpz
    BICYCL::Mpz sid_mpz;
    Fr_to_Mpz(sid_mpz, wt.sid);

    // Pick rho <-$ [0, bound]
    {
        BICYCL::Mpz bound(CL.encrypt_randomness_bound());
        rho = rng.random_mpz(bound);
    }

    // C_sid = CL.Enc(pk_cl, sid; rho)
    BICYCL::CL_HSMqk::CipherText C_sid =
        CL.encrypt(cl_pk, BICYCL::CL_HSMqk::ClearText(CL, sid_mpz), rho);

    // ----------------------------------------------------------
    // Step 2b. Compute B_sid
    //   m_r <-$ Z_q*
    //   B_sid = G1 + sid*H1 + m_r*H2 in G1
    // ----------------------------------------------------------
    wt.m_r.setByCSPRNG();

    {
        mcl::G1 t1, t2;
        mcl::G1::mul(t1, H1, wt.sid);
        mcl::G1::mul(t2, H2, wt.m_r);
        mcl::G1::add(st.B_sid, G1, t1);
        mcl::G1::add(st.B_sid, st.B_sid, t2);
    }

    return C_sid;
}

inline
BICYCL::CL_HSMqk::CipherText CL_BBS_proof(
    const IssueParams&                    pp,
    const mcl_ST&                         st,
    const BICYCL::CL_HSMqk::CipherText&  C_sid,
    const mcl_WT&                         wt,
    const BICYCL::Mpz&                    rho,
    mcl_COM & com_, mcl_RESP& resp_)
{
    const BBSParams&              bbs  = pp.params_;
    const BBS::PublicKey&         bpk  = pp.pk_;
    BICYCL::CL_HSMqk&             CL   = pp.C_;
    BICYCL::RandGen&              rng  = pp.randgen_;
    const BICYCL::CL_HSMqk::PublicKey& cl_pk = pp.cl_pk_;

    const mcl::G1& G1 = bbs.g1;
    const mcl::G1& H1 = bpk.pks_g1[0];
    const mcl::G1& H2 = bpk.pks_g1[1];
    const mcl::G1& H_s = pp.H_s;

    // Step 1a. Pick r <-$ Z_q*
    mcl::Fr r;
    r.setByCSPRNG();

    const mcl::G1& A_ob  = wt.sig.A;
    const mcl::Fr& e_ob  = wt.sig.e;
    const mcl::Fr& r_ob  = wt.r_ob;
    const mcl::Fr& usk   = wt.usk;
    const mcl::Fr& sid   = wt.sid;
    const mcl::Fr& m_r   = wt.m_r;

    // ---- A_bar = r * A_ob ----
    mcl::G1 A_bar, B_bar, N_bar;
    mcl::G1::mul(A_bar, A_ob, r);

    // ---- B_bar = r*G1 + (r*usk)*H1 + (r*r_ob)*H2 + (-e)*A_bar ----
    // Use independent temporaries (no in-place aliasing)
    {
        mcl::Fr r_usk, r_rob, neg_e;
        mcl::Fr::mul(r_usk, r, usk);
        mcl::Fr::mul(r_rob, r, r_ob);
        mcl::Fr::neg(neg_e, e_ob);

        mcl::G1 t1, t2, t3, t4, s1, s2;
        mcl::G1::mul(t1, G1,    r);
        mcl::G1::mul(t2, H1,    r_usk);
        mcl::G1::mul(t3, H2,    r_rob);
        mcl::G1::mul(t4, A_bar, neg_e);
        mcl::G1::add(s1, t1, t2);      // s1 = t1 + t2
        mcl::G1::add(s2, t3, t4);      // s2 = t3 + t4
        mcl::G1::add(B_bar, s1, s2);   // B_bar = s1 + s2 (all distinct)
    }

    mcl::G1::mul(N_bar, st.null, r);

    // Step 1c. Pick masking factors
    //   alpha_r, alpha_usk, alpha_ob, alpha_e, alpha_sid, alpha_mr <-$ Z_q
    //   alpha_rho <-$ [0, 2^{lambda+lambda_s} * D_q]
    mcl::Fr alpha_r, alpha_usk, alpha_ob, alpha_e, alpha_sid, alpha_mr;
    alpha_r.setByCSPRNG();
    alpha_usk.setByCSPRNG();
    alpha_ob.setByCSPRNG();
    alpha_e.setByCSPRNG();
    alpha_sid.setByCSPRNG();
    alpha_mr.setByCSPRNG();

    BICYCL::Mpz alpha_rho;
    {
        BICYCL::Mpz bound(CL.encrypt_randomness_bound());
        BICYCL::Mpz::mulby2k(bound, bound, 256);
        BICYCL::Mpz::mulby2k(bound, bound, 40);
        alpha_rho = rng.random_mpz(bound);
    }

    // Step 1d. Compute the commitment components
    //   U_N1    = alpha_r * null
    //   U_N2    = alpha_usk * H_s
    //   U_barB  = alpha_r*G1 + alpha_e*A_bar + alpha_usk*H1 + alpha_ob*H2
    //   U_Bprime = alpha_sid*H1 + alpha_mr*H2
    //   C'      = (g_q^{alpha_rho}, f^{alpha_sid} * pk_cl^{alpha_rho})
    mcl::G1::mul(com_.U_N1, st.null, alpha_r);
    mcl::G1::mul(com_.U_N2, H_s,    alpha_usk);

    // U_barB = alpha_r*G1 + alpha_e*A_bar + alpha_usk*H1 + alpha_ob*H2
    {
        mcl::G1 t1, t2, t3, t4, s1, s2;
        mcl::G1::mul(t1, G1,    alpha_r);
        mcl::G1::mul(t2, A_bar, alpha_e);
        mcl::G1::mul(t3, H1,    alpha_usk);
        mcl::G1::mul(t4, H2,    alpha_ob);
        mcl::G1::add(s1, t1, t2);
        mcl::G1::add(s2, t3, t4);
        mcl::G1::add(com_.U_barB, s1, s2);
    }

    // U_Bprime = alpha_sid*H1 + alpha_mr*H2
    {
        mcl::G1 t1, t2;
        mcl::G1::mul(t1, H1, alpha_sid);
        mcl::G1::mul(t2, H2, alpha_mr);
        mcl::G1::add(com_.U_Bprime, t1, t2);
    }

    // C': encrypt alpha_sid with randomness alpha_rho
    BICYCL::Mpz alpha_sid_mpz;
    Fr_to_Mpz(alpha_sid_mpz, alpha_sid);
    BICYCL::CL_HSMqk::CipherText C_prime(
        CL, cl_pk, BICYCL::CL_HSMqk::ClearText(CL, alpha_sid_mpz), alpha_rho);

    com_.A_bar = A_bar;
    com_.B_bar = B_bar;
    com_.N_bar = N_bar;

    // Step 2. Fiat-Shamir challenge value ch
    mcl::Fr ch = computeChallenge(st.null, st.B_sid, C_sid, com_, C_prime);
    BICYCL::Mpz ch_mpz;
    Fr_to_Mpz(ch_mpz, ch);

    // Step 3. Compute the responses
    //   beta_r    = alpha_r    + r * ch
    //   beta_usk  = alpha_usk  + r * usk * ch
    //   beta_ob   = alpha_ob   + r * r_ob * ch
    //   beta_e    = alpha_e    - e_ob * ch
    //   beta_sid  = alpha_sid  + sid * ch
    //   beta_mr   = alpha_mr   + m_r * ch
    //   beta_rho  = alpha_rho  + rho * ch   (Mpz operation)
    mcl::Fr ch_r, ch_rusk, ch_usk;
    mcl::Fr::mul(ch_r,    r,   ch);
    mcl::Fr::mul(ch_usk,  usk, ch);
    mcl::Fr::mul(ch_rusk, r,   ch_usk);   // r * usk * ch

    // beta_r = alpha_r + r*ch
    mcl::Fr::add(resp_.beta_r,   alpha_r,   ch_r);
    // beta_usk = alpha_usk + r*usk*ch
    mcl::Fr::add(resp_.beta_usk, alpha_usk, ch_rusk);

    // beta_ob = alpha_ob + r*r_ob*ch
    {
        mcl::Fr tmp;
        mcl::Fr::mul(tmp, r_ob, ch_r);
        mcl::Fr::add(resp_.beta_ob, alpha_ob, tmp);
    }
    // beta_e = alpha_e - e_ob*ch
    {
        mcl::Fr tmp;
        mcl::Fr::mul(tmp, e_ob, ch);
        mcl::Fr::sub(resp_.beta_e, alpha_e, tmp);
    }
    // beta_sid = alpha_sid + sid*ch
    {
        mcl::Fr tmp;
        mcl::Fr::mul(tmp, sid, ch);
        mcl::Fr::add(resp_.beta_sid, alpha_sid, tmp);
    }
    // beta_mr = alpha_mr + m_r*ch
    {
        mcl::Fr tmp;
        mcl::Fr::mul(tmp, m_r, ch);
        mcl::Fr::add(resp_.beta_mr, alpha_mr, tmp);
    }
    // beta_rho = alpha_rho + rho*ch  (Mpz)
    {
        BICYCL::Mpz rho_ch;
        BICYCL::Mpz::mul(rho_ch, rho, ch_mpz);
        BICYCL::Mpz::add(resp_.beta_rho, alpha_rho, rho_ch);
    }
    return C_prime;
}

inline
bool CL_BBS_verify(
    const IssueParams&                   pp,
    const mcl_ST&                        st,
    const BICYCL::CL_HSMqk::CipherText&  C_sid,
    const mcl_COM&                       com_,
    const BICYCL::CL_HSMqk::CipherText&  C_prime_,
    const mcl_RESP&                      resp_)
{
    const BBSParams&              bbs   = pp.params_;
    const BBS::PublicKey&         bpk   = pp.pk_;
    BICYCL::CL_HSMqk&             CL    = pp.C_;
    const BICYCL::CL_HSMqk::PublicKey& cl_pk = pp.cl_pk_;

    const mcl::G1& G1  = bbs.g1;
    const mcl::G2& G2  = bbs.g2;
    const mcl::G1& H1  = bpk.pks_g1[0];
    const mcl::G1& H2  = bpk.pks_g1[1];
    const mcl::G1& H_s = pp.H_s;

    const mcl::G1& A_bar = com_.A_bar;
    const mcl::G1& B_bar = com_.B_bar;
    const mcl::G1& N_bar = com_.N_bar;

    // --------------------------------------------------------
    // a. Check that A_bar is not the identity element (A_bar != 1_{G1})
    // --------------------------------------------------------
    if (A_bar.isZero()) {
        std::cerr << "[Verify Failed] Step (a): A_bar is zero point." << std::endl;
        return false;
    }

    // --------------------------------------------------------
    // Recompute the challenge (Fiat-Shamir heuristic)
    // --------------------------------------------------------
    mcl::Fr ch = computeChallenge(st.null, st.B_sid, C_sid, com_, C_prime_);
    BICYCL::Mpz ch_mpz;
    Fr_to_Mpz(ch_mpz, ch);

    // --------------------------------------------------------
    // b. Range check: beta_rho in [0, 2^{lambda + lambda_s} * D_q]
    // --------------------------------------------------------
    const BICYCL::Mpz& beta_rho_val = resp_.beta_rho;
    {
        BICYCL::Mpz bound(CL.encrypt_randomness_bound());
        BICYCL::Mpz::mulby2k(bound, bound, 256); // security parameter lambda
        BICYCL::Mpz::mulby2k(bound, bound, 40);  // statistical parameter lambda_s
        if (beta_rho_val < BICYCL::Mpz(0L) || beta_rho_val >= bound) {
            std::cerr << "[Verify Failed] Step (b): beta_rho is out of range." << std::endl;
            return false;
        }
    }

    // --------------------------------------------------------
    // c. Pairing check: e(A_bar, X) = e(B_bar, G2)
    //    Verifies the structure of the BBS signature
    // --------------------------------------------------------
    {
        mcl::GT lhs, rhs;
        mcl::pairing(lhs, A_bar, bpk.pk_g2);   // e(A_bar, pk_g2)
        mcl::pairing(rhs, B_bar, G2);          // e(B_bar, g2)
        if (lhs != rhs) {
            std::cerr << "[Verify Failed] Step (c): Pairing e(A_bar, X) == e(B_bar, G2) failed." << std::endl;
            return false;
        }
    }

    // --------------------------------------------------------
    // d. Linear relation: U_N1 + ch * N_bar = beta_r * null
    // --------------------------------------------------------
    {
        mcl::G1 lhs, rhs, ch_Nbar;
        mcl::G1::mul(ch_Nbar, N_bar, ch);
        mcl::G1::add(lhs, com_.U_N1, ch_Nbar);
        mcl::G1::mul(rhs, st.null, resp_.beta_r);
        if (lhs != rhs) {
            std::cerr << "[Verify Failed] Step (d): Equation U_N1 + ch*N_bar == beta_r*null failed." << std::endl;
            return false;
        }
    }

    // --------------------------------------------------------
    // e. Linear relation: U_N2 + ch * N_bar = beta_usk * H_s
    // --------------------------------------------------------
    {
        mcl::G1 lhs, rhs, ch_Nbar;
        mcl::G1::mul(ch_Nbar, N_bar, ch);
        mcl::G1::add(lhs, com_.U_N2, ch_Nbar);
        mcl::G1::mul(rhs, H_s, resp_.beta_usk);
        if (lhs != rhs) {
            std::cerr << "[Verify Failed] Step (e): Equation U_N2 + ch*N_bar == beta_usk*H_s failed." << std::endl;
            return false;
        }
    }

    // --------------------------------------------------------
    // f. Commitment consistency: U_barB + ch * B_bar = ...
    // --------------------------------------------------------
    {
        mcl::G1 lhs, rhs;
        mcl::G1 ch_Bbar;
        mcl::G1::mul(ch_Bbar, B_bar, ch);
        mcl::G1::add(lhs, com_.U_barB, ch_Bbar);

        mcl::G1 t1, t2, t3, t4;
        mcl::G1::mul(t1, G1,    resp_.beta_r);
        mcl::G1::mul(t2, A_bar, resp_.beta_e);
        mcl::G1::mul(t3, H1,    resp_.beta_usk);
        mcl::G1::mul(t4, H2,    resp_.beta_ob);
        mcl::G1::add(rhs, t1, t2);
        mcl::G1::add(rhs, rhs, t3);
        mcl::G1::add(rhs, rhs, t4);

        if (lhs != rhs) {
            std::cerr << "[Verify Failed] Step (f): Linear relation for B_bar failed." << std::endl;
            return false;
        }
    }

    // --------------------------------------------------------
    // g. Commitment consistency: U_Bprime + ch * (B_sid - G1) = ...
    // --------------------------------------------------------
    {
        mcl::G1 lhs, rhs;
        mcl::G1 Bsid_minus_G1;
        mcl::G1 negG1;
        mcl::G1::neg(negG1, G1);
        mcl::G1::add(Bsid_minus_G1, st.B_sid, negG1);

        mcl::G1 ch_tmp;
        mcl::G1::mul(ch_tmp, Bsid_minus_G1, ch);
        mcl::G1::add(lhs, com_.U_Bprime, ch_tmp);

        mcl::G1 t1, t2;
        mcl::G1::mul(t1, H1, resp_.beta_sid);
        mcl::G1::mul(t2, H2, resp_.beta_mr);
        mcl::G1::add(rhs, t1, t2);

        if (lhs != rhs) {
            std::cerr << "[Verify Failed] Step (g): Linear relation for B_prime failed." << std::endl;
            return false;
        }
    }

    // --------------------------------------------------------
    // h. Homomorphic check: Enc(beta_sid; beta_rho) == C' * C_sid^ch
    // --------------------------------------------------------
    {
        BICYCL::Mpz beta_sid_mpz;
        Fr_to_Mpz(beta_sid_mpz, resp_.beta_sid);

        // Compute LHS: re-encrypt using the response values
        BICYCL::CL_HSMqk::CipherText lhs_ct = CL.encrypt(cl_pk,
                       BICYCL::CL_HSMqk::ClearText(CL, beta_sid_mpz),
                       beta_rho_val);

        // Compute RHS: C' * C_sid^ch using homomorphic properties
        BICYCL::CL_HSMqk::CipherText C_sid_ch = CL.scal_ciphertexts(cl_pk, C_sid, ch_mpz, BICYCL::Mpz("0"));
        BICYCL::CL_HSMqk::CipherText rhs_ct = CL.add_ciphertexts(cl_pk, C_prime_, C_sid_ch, BICYCL::Mpz("0"));

        // Compare ciphertext components c1 and c2
        if (!(lhs_ct.c1() == rhs_ct.c1()) || !(lhs_ct.c2() == rhs_ct.c2())) {
            std::cerr << "[Verify Failed] Step (h): CL homomorphic ciphertext equation failed." << std::endl;
            return false;
        }
    }

    // All checks passed: verification succeeds
    return true;
}

inline BICYCL::Mpz factorial(size_t value) {
    BICYCL::Mpz result("1");
    for (size_t current = 2; current <= value; ++current) {
        BICYCL::Mpz::mul(result, result, current);
    }
    return result;
}

inline BICYCL::Mpz cl_lagrange_at_zero(size_t bound,
                                       size_t party_id,
                                       const BICYCL::Mpz& delta) {
    BICYCL::Mpz numerator("1");
    BICYCL::Mpz denominator("1");
    BICYCL::Mpz result;

    for (size_t current_party_id = 1; current_party_id <= bound; current_party_id++) {
        if (current_party_id == party_id) {
            continue;
        }

        BICYCL::Mpz::mul(numerator, numerator, current_party_id);
        if (current_party_id > party_id) {
            BICYCL::Mpz::mul(denominator, denominator, current_party_id - party_id);
        } else {
            BICYCL::Mpz::mul(denominator, denominator, party_id - current_party_id);
            denominator.neg();
        }
    }

    BICYCL::Mpz::divexact(result, delta, denominator);
    BICYCL::Mpz::mul(result, result, numerator);
    return result;
}

inline
void partial_decrypt(const BICYCL::CL_HSMqk&             CL,
                     const BICYCL::CL_HSMqk::SecretKey& secret_key_share,
                     const BICYCL::CL_HSMqk::CipherText& ciphertext,
                     BICYCL::QFI& partial_decryption) {
    BICYCL::Mpz secret_key_mpz(secret_key_share);
    BICYCL::Mpz::mod(secret_key_mpz, secret_key_mpz, CL.secretkey_bound());

    CL.Cl_G().nupow(partial_decryption, ciphertext.c1(), secret_key_mpz);
    if (CL.compact_variant()) {
        CL.from_Cl_DeltaK_to_Cl_Delta(partial_decryption);
    }
}

inline
BICYCL::CL_HSMqk::ClearText aggregate_partial_ciphertext(const BICYCL::CL_HSMqk&             CL, size_t bound,
    const std::vector<BICYCL::QFI>& partial_decryptions,
    const BICYCL::CL_HSMqk::CipherText& ciphertext,
    const BICYCL::Mpz delta)  {
    BICYCL::QFI aggregated_ciphertext = ciphertext.c2();
    for (size_t party_id = 1; party_id <= bound; party_id++) {
        BICYCL::QFI lagrange_component;
        CL.Cl_G().nupow(
            lagrange_component,
            partial_decryptions[party_id - 1],
            cl_lagrange_at_zero(bound, party_id, delta));
        CL.Cl_Delta().nucompinv(
            aggregated_ciphertext, aggregated_ciphertext, lagrange_component);
    }

    return BICYCL::CL_HSMqk::ClearText(
        CL, CL.dlog_in_F(aggregated_ciphertext));
}

#endif // PROVECL_H
