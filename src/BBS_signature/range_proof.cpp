#include "range_proof.h"
#include <mcl/gmp_util.hpp>

Fr inner_product(const std::vector<Fr>& a, const std::vector<Fr>& b) {
    assert(a.size() == b.size());
    Fr res(0);
    for (size_t i = 0; i < a.size(); ++i) {
        Fr tmp;
        Fr::mul(tmp, a[i], b[i]);
        Fr::add(res, res, tmp);
    }
    return res;
}

Fr get_x_challenge(const G1& L, const G1& R, const G1& ux) {
    std::ostringstream output;
    output << L << R << ux;
    Fr challenge;
    challenge.setHashOf(output.str());
    return challenge;
}

InnerProductArg prove(
    const std::vector<G1>& G,
    const std::vector<G1>& H,
    const G1& ux,
    const G1& P,
    const std::vector<Fr>& a,
    const std::vector<Fr>& b,
    std::vector<G1>& L_vec,
    std::vector<G1>& R_vec
) {
    size_t n = G.size();
    if (n == 1) {
        return {L_vec, R_vec, a[0], b[0]};
    }

    size_t half = n / 2;

    // 1. Split vectors (simulating Rust's split_at)
    std::vector<Fr> a_L(a.begin(), a.begin() + half), a_R(a.begin() + half, a.end());
    std::vector<Fr> b_L(b.begin(), b.begin() + half), b_R(b.begin() + half, b.end());
    std::vector<G1> G_L(G.begin(), G.begin() + half), G_R(G.begin() + half, G.end());
    std::vector<G1> H_L(H.begin(), H.begin() + half), H_R(H.begin() + half, H.end());

    // 2. Compute L and R
    // L = MSM(G_R, a_L) + MSM(H_L, b_R) + ux * <a_L, b_R>
    Fr c_L = inner_product(a_L, b_R);
    G1 L, tmp1, tmp2;
    G1::mulVec(L, G_R.data(), a_L.data(), half); // Efficient MSM
    G1::mulVec(tmp1, H_L.data(), b_R.data(), half);
    G1::add(L, L, tmp1);
    G1::mul(tmp2, ux, c_L);
    G1::add(L, L, tmp2);

    // R = MSM(G_L, a_R) + MSM(H_R, b_L) + ux * <a_R, b_L>
    Fr c_R = inner_product(a_R, b_L);
    G1 R;
    G1::mulVec(R, G_L.data(), a_R.data(), half);
    G1::mulVec(tmp1, H_R.data(), b_L.data(), half);
    G1::add(R, R, tmp1);
    G1::mul(tmp2, ux, c_R);
    G1::add(R, R, tmp2);

    L_vec.push_back(L);
    R_vec.push_back(R);

    // 3. Compute challenge value x
    Fr x = get_x_challenge(L, R, ux);
    Fr x_inv;
    Fr::inv(x_inv, x);

    // 4. Update vectors a_new, b_new, G_new, H_new
    std::vector<Fr> a_new(half), b_new(half);
    std::vector<G1> G_new(half), H_new(half);

    for (size_t i = 0; i < half; ++i) {
        // a_i' = a_L_i * x + a_R_i * x_inv
        Fr t1, t2;
        Fr::mul(t1, a_L[i], x);
        Fr::mul(t2, a_R[i], x_inv);
        Fr::add(a_new[i], t1, t2);

        // b_i' = b_R_i * x + b_L_i * x_inv
        Fr::mul(t1, b_R[i], x);
        Fr::mul(t2, b_L[i], x_inv);
        Fr::add(b_new[i], t1, t2);

        // G_i' = G_L_i * x_inv + G_R_i * x
        G1 p1, p2;
        G1::mul(p1, G_L[i], x_inv);
        G1::mul(p2, G_R[i], x);
        G1::add(G_new[i], p1, p2);

        // H_i' = H_L_i * x + H_R_i * x_inv
        G1::mul(p1, H_L[i], x);
        G1::mul(p2, H_R[i], x_inv);
        G1::add(H_new[i], p1, p2);
    }

    return prove(G_new, H_new, ux, P, a_new, b_new, L_vec, R_vec);
}

bool fast_verify(
    const InnerProductArg& proof,
    const std::vector<G1>& g_vec,
    const std::vector<G1>& hi_tag,
    const G1& ux,
    const G1& P
) {
    size_t n = g_vec.size();
    size_t lg_n = proof.L.size();

    std::vector<Fr> x_sq_vec, x_inv_sq_vec;
    Fr all_inv(1);

    for (size_t i = 0; i < lg_n; ++i) {
        Fr x = get_x_challenge(proof.L[i], proof.R[i], ux);
        Fr x_inv; Fr::inv(x_inv, x);
        
        Fr x_sq, x_inv_sq;
        Fr::mul(x_sq, x, x);
        Fr::mul(x_inv_sq, x_inv, x_inv);

        x_sq_vec.push_back(x_sq);
        x_inv_sq_vec.push_back(x_inv_sq);
        Fr::mul(all_inv, all_inv, x_inv);
    }

    // Compute s vector
    std::vector<Fr> s(n);
    s[0] = all_inv;
    for (size_t i = 1; i < n; ++i) {
        size_t lg_i = 31 - __builtin_clz(i); // Get high bit index
        size_t k = 1 << lg_i;
        Fr x_sq = x_sq_vec[(lg_n - 1) - lg_i];
        Fr::mul(s[i], s[i - k], x_sq);
    }

    // Construct scalars and points for combined MSM
    std::vector<Fr> scalars;
    std::vector<G1> points;

    // 1. g_vec section: a_tag * s_i
    for (size_t i = 0; i < n; ++i) {
        Fr tmp; Fr::mul(tmp, proof.a_tag, s[i]);
        scalars.push_back(tmp);
        points.push_back(g_vec[i]);
    }
    // 2. hi_tag section: b_tag * s_i_inv
    for (size_t i = 0; i < n; ++i) {
        Fr s_inv, tmp;
        Fr::inv(s_inv, s[i]);
        Fr::mul(tmp, proof.b_tag, s_inv);
        scalars.push_back(tmp);
        points.push_back(hi_tag[i]);
    }
    // 3. L and R section: -x^2 and -x^-2
    for (size_t i = 0; i < lg_n; ++i) {
        Fr m_x_sq, m_x_inv_sq;
        Fr::neg(m_x_sq, x_sq_vec[i]);
        Fr::neg(m_x_inv_sq, x_inv_sq_vec[i]);
        scalars.push_back(m_x_sq);
        scalars.push_back(m_x_inv_sq);
        points.push_back(proof.L[i]);
        points.push_back(proof.R[i]);
    }

    // Verification: P_expect = MSM(points, scalars) + ux * (a_tag * b_tag)
    G1 expect_P, ux_part;
    G1::mulVec(expect_P, points.data(), scalars.data(), points.size());
    
    Fr c;
    Fr::mul(c, proof.a_tag, proof.b_tag);
    G1::mul(ux_part, ux, c);
    G1::add(expect_P, expect_P, ux_part);

    return expect_P == P;
}

std::vector<Fr> power_vector(Fr y, size_t n) {
    std::vector<Fr> res(n);
    res[0] = 1;
    for (size_t i = 1; i < n; ++i) {
        Fr::mul(res[i], res[i - 1], y);
    }
    return res;
}

// Hash a sequence of G1 points into a scalar
static Fr chainPoints(std::initializer_list<const G1*> pts) {
    std::vector<uint8_t> buf;
    for (auto* p : pts) {
        // G1 Serialization (Compressed format, 48 bytes)
        uint8_t tmp[48];
        size_t n = p->serialize(tmp, sizeof(tmp));
        buf.insert(buf.end(), tmp, tmp + n);
    }
    Fr result;
    result.setHashOf(buf.data(), buf.size());
    return result;
}

// Hash a sequence of Fr scalars into a scalar
static Fr chainScalars(std::initializer_list<const Fr*> scalars) {
    std::vector<uint8_t> buf;
    for (auto* f : scalars) {
        uint8_t tmp[32];
        size_t n = f->serialize(tmp, sizeof(tmp));
        buf.insert(buf.end(), tmp, tmp + n);
    }
    Fr result;
    result.setHashOf(buf.data(), buf.size());
    return result;
}

// Integer power for Fr: out = base^exp
static void frPow(Fr& out, const Fr& base, uint32_t exp) {
    out = 1;
    Fr b = base;
    uint32_t e = exp;
    while (e > 0) {
        if (e & 1) Fr::mul(out, out, b);
        Fr::mul(b, b, b);
        e >>= 1;
    }
}

// Random scalar generation
static Fr randomScalar() {
    Fr r;
    r.setByCSPRNG();
    return r;
}

RangeProof prove_range(
    const std::vector<G1>& g_vec,
    const std::vector<G1>& h_vec,
    const G1& G, 
    const G1& H, 
    std::vector<Fr> secret, // Array of original secret values
    const std::vector<Fr>& blinding, // Corresponding blinding factors
    size_t bit_length
) {
    size_t num_of_proofs = secret.size();
    size_t nm = num_of_proofs * bit_length;

    // ── Random blinding factors ──────────────────────────────────────────────────────
    Fr alpha = randomScalar();
    Fr rho   = randomScalar();

    // A = H * alpha,  S = H * rho
    G1 A, S;
    G1::mul(A, H, alpha);
    G1::mul(S, H, rho);

    // ── Aggregate secrets and expand each bit ─────────────────────────────────────────
    // aL[i] in {0,1}, aR[i] = aL[i] - 1
    std::vector<bool> secret_bits(nm, false);
    std::vector<Fr>   aL(nm), aR(nm);
    {
        Fr one = 1; 
        for (size_t i = 0; i < nm; i++) {
            size_t seg = i / bit_length;   // Corresponds to secret[seg]
            size_t bit = i % bit_length;   // Bit index within segment (LSB = 0)

            // Fr::serialize → little-endian byte array
            uint8_t bytes[32] = {};
            secret[seg].serialize(bytes, sizeof(bytes));
            bool b = (bytes[bit / 8] >> (bit % 8)) & 1;

            secret_bits[i] = b;
            aL[i] = (b ? 1 : 0);
            Fr::sub(aR[i], aL[i], one);   // aR = aL - 1 (mod order)
        }
    }

    // ── Calculate A ────────────────────────────────────────────────────────────
    // A = H*alpha + sum_{aL[i]=1} g_vec[i] - sum_{aL[i]=0} h_vec[i]
    for (size_t i = 0; i < nm; i++) {
        if (secret_bits[i]) {
            G1::add(A, A, g_vec[i]);
        } else {
            G1::sub(A, A, h_vec[i]);
        }
    }

    // ── Random vectors SL, SR, calculate S ───────────────────────────────────────────
    std::vector<Fr> SR(nm), SL(nm);
    for (size_t i = 0; i < nm; i++) {
        SL[i] = randomScalar();
        SR[i] = randomScalar();
    }


    // S = H*rho + sum_i (g_vec[i]*SL[i] + h_vec[i]*SR[i])
    for (size_t i = 0; i < nm; i++) {
        G1 tmp1, tmp2;
        G1::mul(tmp1, g_vec[i], SL[i]);
        G1::mul(tmp2, h_vec[i], SR[i]);
        G1::add(S, S, tmp1);
        G1::add(S, S, tmp2);
    }

    // ── Fiat-Shamir: y, z ─────────────────────────────────────────────────
    Fr y = chainPoints({&A, &S});

    // Rust: yG = base_point * y, then z = hash(yG)
    // Using the generator G passed by caller (consistent with Rust side)
    G1 yG;
    G1::mul(yG, G, y);
    Fr z = chainPoints({&yG});

    // ── yi = [1, y, y^2, ..., y^{nm-1}] ─────────────────────────────────
    std::vector<Fr> yi(nm);
    yi[0] = 1;
    for (size_t i = 1; i < nm; i++) Fr::mul(yi[i], yi[i-1], y);

    // ── t2 = < SL ○ yi, SR > ─────────────────────────────────────────────
    Fr t2; t2 = 0;
    for (size_t i = 0; i < nm; i++) {
        Fr tmp;
        Fr::mul(tmp, SR[i], yi[i]);
        Fr::mul(tmp, tmp, SL[i]);
        Fr::add(t2, t2, tmp);
    }

    // ── vec_2n = [1, 2, 4, ..., 2^{bit_length-1}] ────────────────────────
    std::vector<Fr> vec_2n(bit_length);
    vec_2n[0] = 1;
    {
        Fr two; two = 2;
        for (size_t i = 1; i < bit_length; i++) Fr::mul(vec_2n[i], vec_2n[i-1], two);
    }

    // ── t1 (Linear term coefficient) ──────────────────────────────────────────────────
    Fr t1; t1 = 0;
    for (size_t i = 0; i < nm; i++) {
        uint32_t j = static_cast<uint32_t>(i / bit_length) + 2;
        size_t   k = i % bit_length;

        Fr z_index; frPow(z_index, z, j);

        // t1_1 = aR[i] + z
        Fr t1_1; Fr::add(t1_1, aR[i], z);
        // t1_2 = t1_1 * yi[i]
        Fr t1_2; Fr::mul(t1_2, t1_1, yi[i]);
        // t1_3 = SL[i] * t1_2
        Fr t1_3; Fr::mul(t1_3, SL[i], t1_2);
        // t1_4 = aL[i] - z
        Fr t1_4; Fr::sub(t1_4, aL[i], z);
        // t1_5 = SR[i] * yi[i]
        Fr t1_5; Fr::mul(t1_5, SR[i], yi[i]);
        // t1_6 = t1_4 * t1_5
        Fr t1_6; Fr::mul(t1_6, t1_4, t1_5);
        // t1_7 = z^j * 2^k
        Fr t1_7; Fr::mul(t1_7, z_index, vec_2n[k]);
        // t1_8 = t1_7 * SL[i]
        Fr t1_8; Fr::mul(t1_8, t1_7, SL[i]);
        // t1_68 = t1_6 + t1_8
        Fr t1_68; Fr::add(t1_68, t1_6, t1_8);
        // row = t1_3 + t1_68
        Fr row; Fr::add(row, t1_3, t1_68);
        Fr::add(t1, t1, row);
    }

    // ── Pedersen Commitments T1, T2 ──────────────────────────────────────────────
    Fr tau1 = randomScalar();
    Fr tau2 = randomScalar();

    G1 T1, T2;
    {
        G1 tmp;
        G1::mul(T1, G, t1);
        G1::mul(tmp, H, tau1);
        G1::add(T1, T1, tmp);

        G1::mul(T2, G, t2);
        G1::mul(tmp, H, tau2);
        G1::add(T2, T2, tmp);
    }

    // ── Fiat-Shamir: fs_challenge (= x in the paper) ─────────────────────
    Fr fs_challenge = chainPoints({&T1, &T2, &G, &H});
    Fr fs_challenge_square;
    Fr::mul(fs_challenge_square, fs_challenge, fs_challenge);

    // ── tau_x = x*tau1 + x^2*tau2 + sum_i z^{i+2} * blinding[i] ─────────
    Fr taux_1; Fr::mul(taux_1, fs_challenge, tau1);
    Fr taux_2; Fr::mul(taux_2, fs_challenge_square, tau2);

    // taux_3 = taux_2 + sum_i z^{i+2} * blinding[i], accumulating from taux_2
    Fr taux_3 = taux_2;
    for (size_t i = 0; i < num_of_proofs; i++) {
        Fr z_j; frPow(z_j, z, static_cast<uint32_t>(i + 2));
        Fr tmp; Fr::mul(tmp, z_j, blinding[i]);
        Fr::add(taux_3, taux_3, tmp);
    }

    Fr tau_x; Fr::add(tau_x, taux_1, taux_3);

    // miu = rho * x + alpha
    Fr miu;
    Fr::mul(miu, rho, fs_challenge);
    Fr::add(miu, miu, alpha);

    // ── Lp, Rp vectors ───────────────────────────────────────────────────────
    std::vector<Fr> Lp(nm), Rp(nm);
    for (size_t i = 0; i < nm; i++) {
        // Lp[i] = SL[i]*x + (aL[i] - z)
        Fr Lp_1; Fr::mul(Lp_1, SL[i], fs_challenge);
        Fr Lp_2; Fr::sub(Lp_2, aL[i], z);
        Fr::add(Lp[i], Lp_1, Lp_2);

        // Rp[i] = yi[i] * (z + aR[i] + SR[i]*x)  +  z^j * 2^k
        uint32_t j = static_cast<uint32_t>(i / bit_length) + 2;
        size_t   k = i % bit_length;
        Fr Rp_1; Fr::mul(Rp_1, SR[i], fs_challenge);

        Fr z_index; frPow(z_index, z, j);
        Fr two_to_the_i = vec_2n[k];
        Fr Rp_2; Fr::mul(Rp_2, z_index, two_to_the_i);

        Fr Rp_3; Fr::add(Rp_3, z, aR[i]);
        Fr::add(Rp_3, Rp_3, Rp_1);

        Fr Rp_4; Fr::mul(Rp_4, yi[i], Rp_3);

        Fr::add(Rp[i], Rp_4, Rp_2);
    }
    // tx = < Lp, Rp >
    Fr tx; tx = 0;
    for (size_t i = 0; i < nm; i++) {
        Fr tmp; Fr::mul(tmp, Lp[i], Rp[i]);
        Fr::add(tx, tx, tmp);
    }
    Fr tx_fe = tx;

    // ── challenge_x, construct P ───────────────────────────────────────────────
    Fr challenge_x = chainScalars({&tau_x, &miu, &tx});
    G1 Gx;
    G1::mul(Gx, G, challenge_x);    // Gx = G * challenge_x (= u^{cx})

    G1 P;
    G1::mul(P, Gx, tx_fe);          // P = Gx * tx

    // ── hi_tag[i] = h_vec[i] * yi_inv[i] ────────────────────────────────
    std::vector<Fr> yi_inv(nm);
    for (size_t i = 0; i < nm; i++) Fr::inv(yi_inv[i], yi[i]);

    std::vector<G1> hi_tag(nm);
    for (size_t i = 0; i < nm; i++) G1::mul(hi_tag[i], h_vec[i], yi_inv[i]);

    // P += sum_i g_vec[i] * Lp[i]
    for (size_t i = 0; i < nm; i++) {
        G1 tmp; G1::mul(tmp, g_vec[i], Lp[i]);
        G1::add(P, P, tmp);
    }
    // P += sum_i hi_tag[i] * Rp[i]
    for (size_t i = 0; i < nm; i++) {
        G1 tmp; G1::mul(tmp, hi_tag[i], Rp[i]);
        G1::add(P, P, tmp);
    }

    // ── Call Inner Product Argument sub-algorithm ────────────────────────────────────────────────
    std::vector<G1> L_vec, R_vec;
    InnerProductArg inner_product_proof =
        prove(g_vec, hi_tag, Gx, P, Lp, Rp, L_vec, R_vec);

    return RangeProof { A, S, T1, T2, tau_x, miu, tx_fe, inner_product_proof };
}

bool fast_verify_range(
    const RangeProof& proof,
    const std::vector<G1>& g_vec,
    const std::vector<G1>& h_vec,
    const G1& G,
    const G1& H,
    const std::vector<G1>& ped_com,
    size_t bit_length)
{
    const size_t num_of_proofs = ped_com.size();
    const size_t nm = num_of_proofs * bit_length;

    // ── Fiat-Shamir: y, z ─────────────────────────────────────────────────
    Fr y = chainPoints({&proof.A, &proof.S});

    G1 yG;
    G1::mul(yG, G, y);
    Fr z = chainPoints({&yG});

    Fr z_bn = z;
    Fr z_minus; Fr::neg(z_minus, z);
    Fr z_minus_fe = z_minus;
    Fr z_squared; Fr::mul(z_squared, z, z);

    // ── yi = [1, y, y^2, ..., y^{nm-1}] ─────────────────────────────────
    std::vector<Fr> yi(nm);
    yi[0] = 1;
    for (size_t i = 1; i < nm; i++) Fr::mul(yi[i], yi[i-1], y);

    // scalar_mul_yn = sum(yi)
    Fr scalar_mul_yn = 0;
    for (size_t i = 0; i < nm; i++) Fr::add(scalar_mul_yn, scalar_mul_yn, yi[i]);

    // vec_2n = [1, 2, 4, ..., 2^{bit_length-1}]
    std::vector<Fr> vec_2n(bit_length);
    vec_2n[0] = 1;
    {
        Fr two = 2;
        for (size_t i = 1; i < bit_length; i++) Fr::mul(vec_2n[i], vec_2n[i-1], two);
    }

    // scalar_mul_2n = sum(vec_2n) = 2^bit_length - 1
    Fr scalar_mul_2n = 0;
    for (size_t i = 0; i < bit_length; i++) Fr::add(scalar_mul_2n, scalar_mul_2n, vec_2n[i]);

    // z_cubed_scalar_mul_2n = sum_{i=0}^{n-1} z^{i+3} * scalar_mul_2n
    Fr z_cubed_scalar_mul_2n = 0;
    for (size_t i = 0; i < num_of_proofs; i++) {
        Fr z_j; frPow(z_j, z_bn, static_cast<uint32_t>(i + 3));
        Fr tmp; Fr::mul(tmp, z_j, scalar_mul_2n);
        Fr::add(z_cubed_scalar_mul_2n, z_cubed_scalar_mul_2n, tmp);
    }

    // delta = (z - z^2) * sum(yi) - z_cubed_scalar_mul_2n
    Fr z_minus_zsq; Fr::sub(z_minus_zsq, z_bn, z_squared);
    Fr z_minus_zsq_scalar_mul_yn; Fr::mul(z_minus_zsq_scalar_mul_yn, z_minus_zsq, scalar_mul_yn);
    Fr delta; Fr::sub(delta, z_minus_zsq_scalar_mul_yn, z_cubed_scalar_mul_2n);

    // ── hi_tag[i] = h_vec[i] * yi_inv[i] ────────────────────────────────
    std::vector<Fr> yi_inv(nm);
    for (size_t i = 0; i < nm; i++) Fr::inv(yi_inv[i], yi[i]);

    std::vector<G1> hi_tag(nm);
    for (size_t i = 0; i < nm; i++) G1::mul(hi_tag[i], h_vec[i], yi_inv[i]);

    // ── Fiat-Shamir: fs_challenge ─────────────────────────────────────────
    Fr fs_challenge = chainPoints({&proof.T1, &proof.T2, &G, &H});
    Fr fs_challenge_square; Fr::mul(fs_challenge_square, fs_challenge, fs_challenge);

    // ── Verify Eq 65: G*tx + H*tau_x == G*delta + T1*x + T2*x^2 + sum(ped_com[i]*z^{i+2}) ──
    G1 left_side;
    {
        G1 Gtx, Htaux;
        G1::mul(Gtx,   G, proof.tx);
        G1::mul(Htaux, H, proof.tau_x);
        G1::add(left_side, Gtx, Htaux);
    }

    G1 right_side;
    {
        // ped_com_sum = sum_i ped_com[i] * z^{i+2}
        G1 ped_com_sum;
        bool first = true;
        for (size_t i = 0; i < num_of_proofs; i++) {
            Fr z_2_m; frPow(z_2_m, z_bn, static_cast<uint32_t>(i + 2));
            G1 tmp; G1::mul(tmp, ped_com[i], z_2_m);
            if (first) { ped_com_sum = tmp; first = false; }
            else        G1::add(ped_com_sum, ped_com_sum, tmp);
        }

        G1 Gdelta, Tx, Tx_sq;
        G1::mul(Gdelta, G,          delta);
        G1::mul(Tx,     proof.T1,   fs_challenge);
        G1::mul(Tx_sq,  proof.T2,   fs_challenge_square);

        G1::add(right_side, ped_com_sum, Gdelta);
        G1::add(right_side, right_side,  Tx);
        G1::add(right_side, right_side,  Tx_sq);
    }

    // ── Construct P ────────────────────────────────────────────────────────────
    Fr challenge_x = chainScalars({&proof.tau_x, &proof.miu, &proof.tx});

    G1 Gx;
    G1::mul(Gx, G, challenge_x);        // u^{cx}

    G1 P;
    G1::mul(P, Gx, proof.tx);           // P = Gx * tx

    // P = H*(-miu) + P + A + S*x
    {
        Fr minus_miu; Fr::neg(minus_miu, proof.miu);
        G1 Hmiu, Sx;
        G1::mul(Hmiu, H,       minus_miu);
        G1::mul(Sx,   proof.S, fs_challenge);
        G1::add(P, P,    Hmiu);
        G1::add(P, P,    proof.A);
        G1::add(P, P,    Sx);
    }

    // P1 = P + sum_i hi_tag[i] * (z*yi[i] + z^{j+2} * 2^k)
    //        where j = i/bit_length, k = i%bit_length
    {
        for (size_t i = 0; i < nm; i++) {
            size_t   j = i / bit_length;
            size_t   k = i % bit_length;

            Fr z_yn;  Fr::mul(z_yn, z_bn, yi[i]);

            Fr z_j;   frPow(z_j, z_bn, static_cast<uint32_t>(j + 2));
            Fr z_j_2n; Fr::mul(z_j_2n, z_j, vec_2n[k]);

            Fr zyn_zsq2n; Fr::add(zyn_zsq2n, z_yn, z_j_2n);

            G1 tmp; G1::mul(tmp, hi_tag[i], zyn_zsq2n);
            G1::add(P, P, tmp);
        }
    }

    // P = P + sum_i g_vec[i] * (-z)
    {
        for (size_t i = 0; i < nm; i++) {
            G1 tmp; G1::mul(tmp, g_vec[i], z_minus_fe);
            G1::add(P, P, tmp);
        }
    }

    // ── Call Inner Product Argument verification sub-algorithm ────────────────────────────────────────────
    bool inner_ok = fast_verify(proof.inner_product_proof, g_vec, hi_tag, Gx, P);

    return inner_ok && (left_side == right_side);
}