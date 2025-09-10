#pragma once

#include <algorithm>
#include <numeric>
#include <utility>

#if defined(__has_include)
#  if __has_include(<execution>)
#    include <execution>
#    define GFLOW_HAS_EXECUTION_HEADER 1
#  endif
#endif

// Conservative feature check for execution policies.
// Some libcs ship <execution> but not full parallel backends; we still compile,
// and "fast" just degrades to sequential.
#ifndef GFLOW_HAS_EXECUTION
#  if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
#    define GFLOW_HAS_EXECUTION 1
#  elif defined(GFLOW_HAS_EXECUTION_HEADER)
#    define GFLOW_HAS_EXECUTION 1
#  else
#    define GFLOW_HAS_EXECUTION 0
#  endif
#endif

namespace gflow {

// Call-site tags (intentionally NOT named `seq` in the global/Eigen namespaces)
struct seq_t  { explicit constexpr seq_t(int)  {} };
struct fast_t { explicit constexpr fast_t(int) {} };

// Inline tag objects
inline constexpr seq_t  seq{0};
inline constexpr fast_t fast{0};

// ---- algorithm wrappers ----
// for_each
template<class PolicyTag, class It, class Fn>
inline void for_each(PolicyTag, It first, It last, Fn&& fn) {
#if GFLOW_HAS_EXECUTION
  if constexpr (std::is_same_v<PolicyTag, seq_t>) {
    std::for_each(std::execution::seq, first, last, std::forward<Fn>(fn));
  } else { // fast
#   if defined(_WIN32)
    // Windows: prefer seq unless you've validated par_unseq with your toolchain
    std::for_each(std::execution::seq, first, last, std::forward<Fn>(fn));
#   else
    // Non-Windows: try par_unseq if present; otherwise seq (compiles even if no backend)
    std::for_each(std::execution::par_unseq, first, last, std::forward<Fn>(fn));
#   endif
  }
#else
  // No <execution>; use plain serial algorithm
  (void)sizeof(PolicyTag); // silence unused
  std::for_each(first, last, std::forward<Fn>(fn));
#endif
}

// transform
template<class PolicyTag, class InIt, class OutIt, class Fn>
inline OutIt transform(PolicyTag, InIt first, InIt last, OutIt out, Fn&& fn) {
#if GFLOW_HAS_EXECUTION
  if constexpr (std::is_same_v<PolicyTag, seq_t>) {
    return std::transform(std::execution::seq, first, last, out, std::forward<Fn>(fn));
  } else {
#   if defined(_WIN32)
    return std::transform(std::execution::seq, first, last, out, std::forward<Fn>(fn));
#   else
    return std::transform(std::execution::par_unseq, first, last, out, std::forward<Fn>(fn));
#   endif
  }
#else
  return std::transform(first, last, out, std::forward<Fn>(fn));
#endif
}

// reduce
template<class PolicyTag, class It, class T>
inline T reduce(PolicyTag, It first, It last, T init) {
#if GFLOW_HAS_EXECUTION
  if constexpr (std::is_same_v<PolicyTag, seq_t>) {
    return std::reduce(std::execution::seq, first, last, init);
  } else {
#   if defined(_WIN32)
    return std::reduce(std::execution::seq, first, last, init);
#   else
    return std::reduce(std::execution::par_unseq, first, last, init);
#   endif
  }
#else
  return std::accumulate(first, last, init);
#endif
}

// transform_reduce (binary)
template<class PolicyTag, class It1, class It2, class T, class BinOp1, class BinOp2>
inline T transform_reduce(PolicyTag, It1 f1, It1 l1, It2 f2, T init, BinOp1 reduce_op, BinOp2 xform_op) {
#if GFLOW_HAS_EXECUTION
  if constexpr (std::is_same_v<PolicyTag, seq_t>) {
    return std::transform_reduce(std::execution::seq, f1, l1, f2, init,
                                 std::forward<BinOp1>(reduce_op), std::forward<BinOp2>(xform_op));
  } else {
#   if defined(_WIN32)
    return std::transform_reduce(std::execution::seq, f1, l1, f2, init,
                                 std::forward<BinOp1>(reduce_op), std::forward<BinOp2>(xform_op));
#   else
    return std::transform_reduce(std::execution::par_unseq, f1, l1, f2, init,
                                 std::forward<BinOp1>(reduce_op), std::forward<BinOp2>(xform_op));
#   endif
  }
#else
  for (; f1 != l1; ++f1, ++f2) {
    init = reduce_op(init, xform_op(*f1, *f2));
  }
  return init;
#endif
}

} // namespace gflow

// Optional convenience macros if you prefer a single "default policy" name:
#if defined(_WIN32)
#  define GFLOW_EXEC_POLICY gflow::seq
#else
#  define GFLOW_EXEC_POLICY gflow::fast
#endif
