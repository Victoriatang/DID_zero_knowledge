#include <mcl/bls12_381.hpp>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <string>
#include <chrono>
#include <execution>

using namespace mcl;

// ── Time tool ──────────────────────────────────────────────────────────────────
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

// ── G1 Serialization ────────────────────
static constexpr size_t G1_BYTES = 48;

static void g1ToBytes(const G1& p, uint8_t out[G1_BYTES]) {
    size_t n = p.serialize(out, G1_BYTES);
    assert(n == G1_BYTES);
}

static G1 bytesToG1(const uint8_t in[G1_BYTES]) {
    G1 p;
    size_t n = p.deserialize(in, G1_BYTES);
    assert(n == G1_BYTES);
    return p;
}

static bool g1Less(const G1& a, const G1& b) {
    uint8_t ba[G1_BYTES], bb[G1_BYTES];
    g1ToBytes(a, ba);
    g1ToBytes(b, bb);
    return std::memcmp(ba, bb, G1_BYTES) < 0;
}

static bool g1Equal(const G1& a, const G1& b) {
    uint8_t ba[G1_BYTES], bb[G1_BYTES];
    g1ToBytes(a, ba);
    g1ToBytes(b, bb);
    return std::memcmp(ba, bb, G1_BYTES) == 0;
}


static void saveAllBinary(const std::string& filename, const std::vector<G1>& arr) {
    std::ofstream ofs(filename, std::ios::binary);
    for (const auto& p : arr) {
        uint8_t buf[G1_BYTES];
        g1ToBytes(p, buf);
        ofs.write(reinterpret_cast<const char*>(buf), G1_BYTES);
    }
}

// Directly jump to the idx-th record using seekg and read only 48 bytes
static G1 readG1AtBinary(std::ifstream& ifs, size_t idx) {
    ifs.seekg(static_cast<std::streamoff>(idx * G1_BYTES));
    uint8_t buf[G1_BYTES];
    ifs.read(reinterpret_cast<char*>(buf), G1_BYTES);
    return bytesToG1(buf);
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

// In-memory binary search (for performance comparison)
int binarySearchMem(const std::vector<G1>& arr, const G1& target) {
    int left = 0, right = static_cast<int>(arr.size()) - 1;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (g1Equal(arr[mid], target))     return mid;
        else if (g1Less(arr[mid], target)) left = mid + 1;
        else                               right = mid - 1;
    }
    return -1;
}

int main() {
    initPairing(mcl::BLS12_381);

    const size_t N = 1000000;
    const std::string filename = "g1_arr_1000000.bin";

    G1 base;
    mcl::hashAndMapToG1(base, "base", 4);
    
    /*// ── 1. Generate, Sort, and Save ──────────────────────────────────────────────────
    std::vector<G1> arr(N);
    {
        Timer t("generate");
        for (size_t i = 0; i < N; i++) {
            Fr r; r.setByCSPRNG();
            G1::mul(arr[i], base, r);
        }
    }
    {
        Timer t("sort");
        std::sort(arr.begin(), arr.end(), g1Less);
    }
    {
        Timer t("save binary");
        saveAllBinary(filename, arr);
    }
    printf("  file size: %.1f MB\n\n",
           static_cast<double>(N * G1_BYTES) / 1e6);
           
*/
    // ── 2. Pick test targets (already in memory, lookup time not counted) ────────────────────────
    std::vector<size_t> test_indices = {0, N/4, N/2, N*3/4, N-1};
    std::vector<G1> targets;
    std::ifstream ifs(filename);
    for (size_t i : test_indices) targets.push_back(readG1AtBinary(ifs,i));

    // Absent point
    Fr r_absent; r_absent.setByCSPRNG();
    G1 absent; G1::mul(absent, base, r_absent);

    // ── 3. File-based binary search timing ────────────────────────────────────────────────
    printf("=== binary search (file, fixed-width binary) ===\n");
    for (size_t k = 0; k < targets.size(); k++) {
        int idx;
        {
            Timer t("  seek index=" + std::to_string(test_indices[k]));
            idx = binarySearchBinary(filename, N, targets[k]);
        }
        assert(idx != -1);
        printf("  -> found at %d\n\n", idx);
    }
    {
        int idx;
        {
            Timer t("  seek (absent)");
            idx = binarySearchBinary(filename, N, absent);
        }
        printf("  -> %s\n\n", idx == -1 ? "not found (expected)" : "found (collision)");
    }
    /*

    // ── 4. In-memory binary search comparison ────────────────────────────────────────────────
    printf("=== binary search (in-memory, for comparison) ===\n");
    for (size_t k = 0; k < targets.size(); k++) {
        int idx;
        {
            Timer t("  mem index=" + std::to_string(test_indices[k]));
            idx = binarySearchMem(arr, targets[k]);
        }
        assert(idx != -1);
        printf("  -> found at %d\n\n", idx);
    }
    {
        int idx;
        {
            Timer t("  mem (absent)");
            idx = binarySearchMem(arr, absent);
        }
        printf("  -> %s\n\n", idx == -1 ? "not found (expected)" : "found (collision)");
    }
    */    

    return 0;
}