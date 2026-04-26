#ifndef THRESHOLD_CL_H
#define THRESHOLD_CL_H

#include <bicycl.hpp>

namespace BICYCL {


class CL_HSMqk_Part_Dec_ZKProof {
public:
    CL_HSMqk_Part_Dec_ZKProof(const CL_HSMqk& cryptosystem,
                             HashAlgo& hash_algo,
                             const CL_HSMqk::PublicKey& public_key,
                             const CL_HSMqk::CipherText& ciphertext,
                             const QFI& partial_decryption,
                             const CL_HSMqk::SecretKey& secret_key,
                             RandGen& random_generator) {
        const int soundness = hash_algo.digest_nbits();

        Mpz randomness_bound(cryptosystem.encrypt_randomness_bound());
        Mpz::mulby2k(randomness_bound, randomness_bound, static_cast<size_t>(soundness));
        Mpz::mulby2k(randomness_bound, randomness_bound, cryptosystem.lambda_distance());

        Mpz witness(random_generator.random_mpz(randomness_bound));
        cryptosystem.power_of_h(t1_, witness);
        cryptosystem.Cl_G().nupow(t2_, ciphertext.c1(), witness);

        k_ = challenge_from_hash(hash_algo, public_key, ciphertext, partial_decryption, t1_, t2_);
        Mpz::mul(z_, k_, secret_key);
        Mpz::add(z_, z_, witness);
    }

    bool verify(const CL_HSMqk& cryptosystem,
                HashAlgo&,
                const CL_HSMqk::PublicKey& public_key,
                const CL_HSMqk::CipherText& ciphertext,
                const QFI& partial_decryption)  {
        QFI lhs_1;
        QFI rhs_1;
        cryptosystem.power_of_h(lhs_1, z_);

        public_key.exponentiation(cryptosystem, rhs_1, k_);
        if (cryptosystem.compact_variant()) {
            cryptosystem.from_Cl_DeltaK_to_Cl_Delta(rhs_1);
        }
        cryptosystem.Cl_G().nucomp(rhs_1, rhs_1, t1_);
        if (!(lhs_1 == rhs_1)) {
            return false;
        }

        QFI lhs_2;
        QFI rhs_2;
        cryptosystem.Cl_G().nupow(lhs_2, ciphertext.c1(), z_);
        cryptosystem.Cl_G().nupow(rhs_2, partial_decryption, k_);
        cryptosystem.Cl_G().nucomp(rhs_2, rhs_2, t2_);
        return lhs_2 == rhs_2;
    }

private:
    static Mpz challenge_from_hash(HashAlgo& hash_algo,
                                   const CL_HSMqk::PublicKey& public_key,
                                   const CL_HSMqk::CipherText& ciphertext,
                                   const QFI& partial_decryption,
                                   const QFI& t1,
                                   const QFI& t2) {
        return Mpz(hash_algo(ciphertext, public_key, partial_decryption, t1, t2));
    }

    QFI t1_;
    QFI t2_;
    Mpz k_;
    Mpz z_;
};
} // namespace BICYCL

#endif //
