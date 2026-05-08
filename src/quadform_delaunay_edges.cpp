#include <Rcpp.h>

#include <algorithm>
#include <cstdio>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include "qhull/libqhull_r.h"
}

namespace {

FILE* qhull_tmpfile(const char* stream_name) {
    FILE* fp = std::tmpfile();
    if (fp == nullptr) {
        throw std::runtime_error(std::string("Failed to create temporary Qhull ") +
                                 stream_name + " stream.");
    }
    return fp;
}

std::vector<char> qhull_command(const std::string& options) {
    std::string command = "qhull d Qbb T0";
    if (!options.empty()) {
        command += " ";
        command += options;
    }
    std::vector<char> out(command.begin(), command.end());
    out.push_back('\0');
    return out;
}

Rcpp::IntegerMatrix edge_set_to_matrix(const std::set<std::pair<int, int> >& edges) {
    Rcpp::IntegerMatrix out(edges.size(), 2);
    int row = 0;
    for (const auto& edge : edges) {
        out(row, 0) = edge.first;
        out(row, 1) = edge.second;
        ++row;
    }
    return out;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List rcpp_quadform_delaunay_edges_3d(const Rcpp::NumericMatrix& X,
                                           std::string qhull_options = "Qt Qbb Qc") {
    if (X.ncol() != 3) {
        Rcpp::stop("'X' must have exactly three columns.");
    }
    if (X.nrow() < 4) {
        Rcpp::stop("'X' must contain at least four points for a 3D Delaunay tessellation.");
    }
    for (int i = 0; i < X.nrow(); ++i) {
        for (int j = 0; j < X.ncol(); ++j) {
            if (!R_finite(X(i, j))) {
                Rcpp::stop("'X' cannot contain NA, NaN, or Inf values.");
            }
        }
    }

    const int dim = 3;
    const int n = X.nrow();
    std::vector<coordT> points(static_cast<std::size_t>(n) * dim);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < dim; ++j) {
            points[static_cast<std::size_t>(i) * dim + j] = X(i, j);
        }
    }

    FILE* outfile = nullptr;
    FILE* errfile = nullptr;
    qhT qh_qh;
    qhT* qh = &qh_qh;
    int exitcode = 0;
    int curlong = 0;
    int totlong = 0;

    auto close_streams = [&]() {
        if (outfile != nullptr) {
            std::fclose(outfile);
            outfile = nullptr;
        }
        if (errfile != nullptr) {
            std::fclose(errfile);
            errfile = nullptr;
        }
    };

    try {
        outfile = qhull_tmpfile("stdout");
        errfile = qhull_tmpfile("stderr");
        QHULL_LIB_CHECK
        qh_zero(qh, errfile);
        std::vector<char> command = qhull_command(qhull_options);
        exitcode = qh_new_qhull(
            qh, dim, n, points.data(), False,
            command.data(), outfile, errfile
        );
        if (exitcode != 0) {
            qh_freeqhull(qh, !qh_ALL);
            qh_memfreeshort(qh, &curlong, &totlong);
            close_streams();
            Rcpp::stop("Qhull failed with exit code %d.", exitcode);
        }

        facetT* facet;
        vertexT* vertex;
        vertexT** vertexp;
        std::vector<std::vector<int> > tetrahedra;
        std::set<std::pair<int, int> > edge_set;

        FORALLfacets {
            if (facet->upperdelaunay) {
                continue;
            }
            if (!facet->simplicial) {
                qh_freeqhull(qh, !qh_ALL);
                qh_memfreeshort(qh, &curlong, &totlong);
                close_streams();
                Rcpp::stop("Qhull returned a non-simplicial Delaunay facet; try different options.");
            }
            std::vector<int> simplex;
            simplex.reserve(4);
            FOREACHvertex_(facet->vertices) {
                const int id = qh_pointid(qh, vertex->point) + 1;
                if (id < 1 || id > n) {
                    qh_freeqhull(qh, !qh_ALL);
                    qh_memfreeshort(qh, &curlong, &totlong);
                    close_streams();
                    Rcpp::stop("Qhull returned a vertex id outside the input point range.");
                }
                simplex.push_back(id);
            }
            if (simplex.size() != 4) {
                qh_freeqhull(qh, !qh_ALL);
                qh_memfreeshort(qh, &curlong, &totlong);
                close_streams();
                Rcpp::stop("Qhull returned a non-tetrahedral Delaunay simplex.");
            }
            std::sort(simplex.begin(), simplex.end());
            tetrahedra.push_back(simplex);
            for (std::size_t a = 0; a < simplex.size(); ++a) {
                for (std::size_t b = a + 1; b < simplex.size(); ++b) {
                    edge_set.insert(std::make_pair(simplex[a], simplex[b]));
                }
            }
        }

        std::sort(tetrahedra.begin(), tetrahedra.end());
        Rcpp::IntegerMatrix tetra_matrix(tetrahedra.size(), 4);
        for (std::size_t i = 0; i < tetrahedra.size(); ++i) {
            for (int j = 0; j < 4; ++j) {
                tetra_matrix(static_cast<int>(i), j) = tetrahedra[i][j];
            }
        }
        Rcpp::IntegerMatrix edge_matrix = edge_set_to_matrix(edge_set);

        qh_freeqhull(qh, !qh_ALL);
        qh_memfreeshort(qh, &curlong, &totlong);
        close_streams();

        return Rcpp::List::create(
            Rcpp::Named("edge_matrix") = edge_matrix,
            Rcpp::Named("tetrahedra") = tetra_matrix,
            Rcpp::Named("n_edges") = edge_matrix.nrow(),
            Rcpp::Named("n_tetrahedra") = tetra_matrix.nrow(),
            Rcpp::Named("qhull_options") = qhull_options,
            Rcpp::Named("qhull_command") = std::string("qhull d Qbb T0 ") + qhull_options,
            Rcpp::Named("exit_code") = exitcode,
            Rcpp::Named("curlong") = curlong,
            Rcpp::Named("totlong") = totlong
        );
    } catch (...) {
        close_streams();
        throw;
    }
}
