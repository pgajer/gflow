#include "kNN.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

enum class knn_cache_mode_t : int {
    none = 0,
    read = 1,
    write = 2,
    readwrite = 3
};

enum class knn_cache_load_status_t : int {
    loaded = 0,
    not_found = 1,
    invalid = 2,
    io_error = 3
};

struct knn_cache_header_t {
    char magic[8];
    std::uint32_t version;
    std::uint32_t endian_marker;
    std::uint64_t n_points;
    std::uint64_t n_features;
    std::uint64_t k;
};

constexpr std::array<char, 8> k_knn_cache_magic = {'G', 'F', 'L', 'K', 'N', 'N', '0', '1'};
constexpr std::uint32_t k_knn_cache_version = 1U;
constexpr std::uint32_t k_knn_cache_endian_marker = 0x01020304U;
std::atomic<std::uint64_t> g_knn_cache_tmp_counter(0U);

template <typename T>
bool write_binary(std::ofstream& out, const T& value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    return static_cast<bool>(out);
}

std::string make_knn_cache_tmp_path(const std::string& cache_path) {
    const std::uint64_t ts = static_cast<std::uint64_t>(
        std::chrono::steady_clock::now().time_since_epoch().count()
    );
    const std::uint64_t counter = g_knn_cache_tmp_counter.fetch_add(1U, std::memory_order_relaxed);
    std::ostringstream oss;
    oss << cache_path << ".tmp." << ts << "." << counter;
    return oss.str();
}

knn_result_t compute_knn_from_eigen_no_cache(
    const Eigen::SparseMatrix<double>& X,
    int k
) {
    const int n = static_cast<int>(X.rows());
    const int d = static_cast<int>(X.cols());

    if (n <= 0 || d <= 0) {
        throw std::invalid_argument("Matrix dimensions must be positive");
    }
    if (k <= 0 || k > n) {
        throw std::invalid_argument("k must be in range [1, n]");
    }

    // Convert Eigen sparse matrix to dense vector-of-vectors format
    // required by kNN() function.
    std::vector<std::vector<double>> X_dense(static_cast<size_t>(n), std::vector<double>(static_cast<size_t>(d)));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < d; ++j) {
            X_dense[static_cast<size_t>(i)][static_cast<size_t>(j)] = X.coeff(i, j);
        }
    }

    return kNN(X_dense, k);
}

knn_cache_load_status_t read_knn_cache_file_flat(
    const std::string& cache_path,
    int expected_n_points,
    int expected_n_features,
    int required_k,
    knn_result_t& out_knn_results,
    std::string& out_reason) {

    std::ifstream in(cache_path, std::ios::binary);
    if (!in.is_open()) {
        if (errno == ENOENT) {
            out_reason = "cache file not found";
            return knn_cache_load_status_t::not_found;
        }
        out_reason = std::string("cannot open cache file: ") + std::strerror(errno);
        return knn_cache_load_status_t::io_error;
    }

    knn_cache_header_t header{};
    if (!in.read(reinterpret_cast<char*>(&header), sizeof(knn_cache_header_t))) {
        out_reason = "cache header is truncated";
        return knn_cache_load_status_t::invalid;
    }

    if (!std::equal(k_knn_cache_magic.begin(), k_knn_cache_magic.end(), header.magic)) {
        out_reason = "cache magic mismatch";
        return knn_cache_load_status_t::invalid;
    }
    if (header.version != k_knn_cache_version) {
        out_reason = "cache version mismatch";
        return knn_cache_load_status_t::invalid;
    }
    if (header.endian_marker != k_knn_cache_endian_marker) {
        out_reason = "cache endianness mismatch";
        return knn_cache_load_status_t::invalid;
    }
    if (header.n_points != static_cast<std::uint64_t>(expected_n_points)) {
        out_reason = "cache nrow mismatch";
        return knn_cache_load_status_t::invalid;
    }
    if (header.n_features != static_cast<std::uint64_t>(expected_n_features)) {
        out_reason = "cache ncol mismatch";
        return knn_cache_load_status_t::invalid;
    }
    if (header.k < static_cast<std::uint64_t>(required_k)) {
        out_reason = "cache k is smaller than required k";
        return knn_cache_load_status_t::invalid;
    }
    if (header.k > static_cast<std::uint64_t>(INT_MAX)) {
        out_reason = "cache k exceeds INT_MAX";
        return knn_cache_load_status_t::invalid;
    }

    const int cached_k = static_cast<int>(header.k);
    const size_t expected_size = static_cast<size_t>(expected_n_points) * static_cast<size_t>(cached_k);

    knn_result_t loaded;
    loaded.n = expected_n_points;
    loaded.k = cached_k;
    loaded.indices.resize(expected_size);
    loaded.distances.resize(expected_size);

    in.read(reinterpret_cast<char*>(loaded.indices.data()),
            static_cast<std::streamsize>(expected_size * sizeof(int)));
    if (!in) {
        out_reason = "cache file is truncated (indices)";
        return knn_cache_load_status_t::invalid;
    }

    in.read(reinterpret_cast<char*>(loaded.distances.data()),
            static_cast<std::streamsize>(expected_size * sizeof(double)));
    if (!in) {
        out_reason = "cache file is truncated (distances)";
        return knn_cache_load_status_t::invalid;
    }

    char trailing = '\0';
    if (in.read(&trailing, 1)) {
        out_reason = "cache file has unexpected trailing bytes";
        return knn_cache_load_status_t::invalid;
    }

    for (size_t off = 0; off < expected_size; ++off) {
        const int idx = loaded.indices[off];
        const double dist = loaded.distances[off];
        if (idx < 0 || idx >= expected_n_points) {
            out_reason = "cache contains out-of-range neighbor index";
            return knn_cache_load_status_t::invalid;
        }
        if (!std::isfinite(dist) || dist < 0.0) {
            out_reason = "cache contains invalid distance";
            return knn_cache_load_status_t::invalid;
        }
    }

    out_knn_results = std::move(loaded);
    out_reason.clear();
    return knn_cache_load_status_t::loaded;
}

void write_knn_cache_file_atomic_flat(
    const std::string& cache_path,
    const knn_result_t& knn_results,
    int n_features) {

    if (cache_path.empty()) {
        throw std::runtime_error("knn.cache.path cannot be empty");
    }

    const size_t expected_size = static_cast<size_t>(knn_results.n) * static_cast<size_t>(knn_results.k);
    if (knn_results.indices.size() != expected_size || knn_results.distances.size() != expected_size) {
        throw std::runtime_error("kNN results have inconsistent dimensions for cache write");
    }

    const std::string tmp_path = make_knn_cache_tmp_path(cache_path);

    knn_cache_header_t header{};
    std::memcpy(header.magic, k_knn_cache_magic.data(), k_knn_cache_magic.size());
    header.version = k_knn_cache_version;
    header.endian_marker = k_knn_cache_endian_marker;
    header.n_points = static_cast<std::uint64_t>(knn_results.n);
    header.n_features = static_cast<std::uint64_t>(n_features);
    header.k = static_cast<std::uint64_t>(knn_results.k);

    {
        std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
        if (!out.is_open()) {
            throw std::runtime_error(
                std::string("Failed to create temp cache file '") + tmp_path + "': " + std::strerror(errno)
            );
        }

        if (!write_binary(out, header)) {
            std::remove(tmp_path.c_str());
            throw std::runtime_error("Failed to write cache header");
        }

        out.write(reinterpret_cast<const char*>(knn_results.indices.data()),
                  static_cast<std::streamsize>(expected_size * sizeof(int)));
        if (!out) {
            std::remove(tmp_path.c_str());
            throw std::runtime_error("Failed to write cache indices");
        }

        out.write(reinterpret_cast<const char*>(knn_results.distances.data()),
                  static_cast<std::streamsize>(expected_size * sizeof(double)));
        if (!out) {
            std::remove(tmp_path.c_str());
            throw std::runtime_error("Failed to write cache distances");
        }

        out.flush();
        if (!out) {
            std::remove(tmp_path.c_str());
            throw std::runtime_error("Failed to flush temp cache file");
        }
    }

    if (std::rename(tmp_path.c_str(), cache_path.c_str()) != 0) {
        const int rename_errno = errno;
#ifdef _WIN32
        std::remove(cache_path.c_str());
        if (std::rename(tmp_path.c_str(), cache_path.c_str()) == 0) {
            return;
        }
#endif
        std::remove(tmp_path.c_str());
        throw std::runtime_error(
            std::string("Failed to atomically move cache file to '") + cache_path + "': " +
            std::strerror(rename_errno)
        );
    }
}

} // namespace

knn_result_t compute_knn_from_eigen(
    const Eigen::SparseMatrix<double>& X,
    int k
) {
    return compute_knn_from_eigen_no_cache(X, k);
}

knn_result_t compute_knn_from_eigen(
    const Eigen::SparseMatrix<double>& X,
    int k,
    const std::string& knn_cache_path,
    int knn_cache_mode_raw,
    bool* cache_hit,
    bool* cache_written
) {
    if (cache_hit != nullptr) {
        *cache_hit = false;
    }
    if (cache_written != nullptr) {
        *cache_written = false;
    }

    if (knn_cache_mode_raw < static_cast<int>(knn_cache_mode_t::none) ||
        knn_cache_mode_raw > static_cast<int>(knn_cache_mode_t::readwrite)) {
        throw std::invalid_argument("knn.cache.mode must be in {0,1,2,3}");
    }
    const knn_cache_mode_t knn_cache_mode = static_cast<knn_cache_mode_t>(knn_cache_mode_raw);

    if (knn_cache_mode == knn_cache_mode_t::none) {
        return compute_knn_from_eigen_no_cache(X, k);
    }

    if (knn_cache_path.empty()) {
        throw std::invalid_argument("knn.cache.path must be provided when knn.cache.mode != 'none'");
    }

    const int n_points = static_cast<int>(X.rows());
    const int n_features = static_cast<int>(X.cols());

    knn_result_t knn_results;
    bool local_cache_hit = false;
    bool local_cache_written = false;

    if (knn_cache_mode == knn_cache_mode_t::read || knn_cache_mode == knn_cache_mode_t::readwrite) {
        std::string cache_reason;
        const knn_cache_load_status_t cache_status = read_knn_cache_file_flat(
            knn_cache_path,
            n_points,
            n_features,
            k,
            knn_results,
            cache_reason
        );

        if (cache_status == knn_cache_load_status_t::loaded) {
            local_cache_hit = true;
        } else if (cache_status == knn_cache_load_status_t::not_found &&
                   knn_cache_mode == knn_cache_mode_t::readwrite) {
            local_cache_hit = false;
        } else {
            throw std::runtime_error(
                std::string("Failed to read kNN cache '") + knn_cache_path + "': " + cache_reason
            );
        }
    }

    if (!local_cache_hit) {
        knn_results = compute_knn_from_eigen_no_cache(X, k);
        if (knn_cache_mode == knn_cache_mode_t::write || knn_cache_mode == knn_cache_mode_t::readwrite) {
            write_knn_cache_file_atomic_flat(knn_cache_path, knn_results, n_features);
            local_cache_written = true;
        }
    }

    if (cache_hit != nullptr) {
        *cache_hit = local_cache_hit;
    }
    if (cache_written != nullptr) {
        *cache_written = local_cache_written;
    }

    return knn_results;
}
