#include <gmpxx.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <utility>

int main() {
    size_t size;
    if (!(std::cin >> size) || size == 0) return 1;
    std::vector<std::vector<mpz_class>> matrix(size, std::vector<mpz_class>(size + 1));
    for (auto& row : matrix) {
        for (auto& entry : row) {
            if (!(std::cin >> entry)) return 1;
        }
    }

    mpz_class previous = 1;
    mpz_class buffer;
    for (size_t k = 0; k + 1 < size; ++k) {
        size_t pivot = k;
        while (pivot < size && matrix[pivot][k] == 0) ++pivot;
        if (pivot == size) throw std::runtime_error("singular block");
        std::swap(matrix[k], matrix[pivot]);
        for (size_t i = k + 1; i < size; ++i) {
            for (size_t j = k + 1; j <= size; ++j) {
                buffer = matrix[k][k] * matrix[i][j];
                mpz_submul(buffer.get_mpz_t(), matrix[i][k].get_mpz_t(), matrix[k][j].get_mpz_t());
                mpz_divexact(matrix[i][j].get_mpz_t(), buffer.get_mpz_t(), previous.get_mpz_t());
            }
        }
        previous = matrix[k][k];
    }
    std::vector<mpq_class> solution(size);
    for (size_t k = size; k-- > 0;) {
        if (matrix[k][k] == 0) throw std::runtime_error("singular block");
        solution[k] = matrix[k][size];
        for (size_t j = k + 1; j < size; ++j) solution[k] -= matrix[k][j] * solution[j];
        solution[k] /= matrix[k][k];
    }
    for (const auto& entry : solution) std::cout << entry << '\n';
}
