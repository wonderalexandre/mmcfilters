#pragma once

//#define NDEBUG  // Remove os asserts do código
#include <cassert>
#include <cstdint>
#include <list>
#include <unordered_set>
#include <unordered_map>

#include <memory>
#include <limits>
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <span>
#include <iostream>


#include "Image.hpp"
#include "PixelSetManager.hpp"
#include "../dataStructure/FastStack.hpp"
#include "../dataStructure/FastQueue.hpp"


namespace mmcfilters {


#define PRINT_LOG 0
#define PRINT_DEBUG 0

using NodeId = int;
constexpr NodeId InvalidNode = -1;  // ou std::numeric_limits<NodeId>::max()
inline bool isValidNode(NodeId id) noexcept {
    return id != InvalidNode;
}


/**
 * @brief Estrutura de dados para marcação eficiente (visited flags) usando carimbos de geração.
 *
 * A `GenerationStampSet` mantém um array de inteiros (stamps), cada posição
 * associada a um índice de elemento (ex.: nó de grafo). Em vez de limpar o
 * array inteiro a cada iteração, um contador de geração (`cur`) é incrementado
 * e usado como "marca lógica". 
 *
 *
 * @code
 * GenerationStampSet visited(numNodes);
 *
 * visited.mark(nodeIdx);
 *
 * if (!visited.isMarked(otherIdx)) {
 *     // processa nó não visitado
 * }
 *
 * visited.resetAll();  // O(1) para preparar nova iteração
 * @endcode
 */
struct GenerationStampSet {
    using gen_t = uint32_t;

    std::unique_ptr<gen_t[]> stamp; // array de carimbos
    size_t n{0};                    // tamanho
    gen_t cur{1};                   // geração atual (0 = “limpo”)

    GenerationStampSet() = default;
    explicit GenerationStampSet(size_t n) { resize(n); }

    void resize(size_t newN) {
        n = newN;
        stamp = std::make_unique<gen_t[]>(n);
        std::fill_n(stamp.get(), n, 0);
        cur = 1;
    }

    inline void mark(size_t idx) noexcept {
        stamp[idx] = cur;
    }

    inline bool isMarked(size_t idx) const noexcept {
        return stamp[idx] == cur;
    }

    // reset lógico em O(1)
    void resetAll() {
        if (++cur == 0) {
            std::fill_n(stamp.get(), n, 0);
            cur = 1;
        }
    }

    // limpeza forçada em O(N)
    void clearAll() {
        std::fill_n(stamp.get(), n, 0);
        cur = 1;
    }

    gen_t generation() const noexcept { return cur; }
};


} // namespace mmcfilters

