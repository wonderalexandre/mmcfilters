#pragma once


// ---------------------------------------------------------------------------
// Controle de assertivas
// ---------------------------------------------------------------------------
#if defined(MMCFILTERS_ENABLE_ASSERTS)
#  ifdef NDEBUG
#    undef NDEBUG
#  endif
#endif
#include <cassert>      // assert()

// ---------------------------------------------------------------------------
// Estruturas de dados (STL containers e algoritmos)
// ---------------------------------------------------------------------------
#include <list>          // Lista duplamente ligada
#include <vector>        // Vetor dinâmico redimensionável
#include <array>         // Array de tamanho fixo em tempo de compilação
#include <deque>         // Deque (fila dupla)
#include <stack>         // Pilha adaptada (baseada em deque por padrão)
#include <unordered_set> // Conjunto hash (não ordenado, busca O(1) médio)
#include <unordered_map> // Mapa hash (não ordenado, busca O(1) médio)
#include <typeindex>     // std::type_index, permite comparar/hashear typeid para usar em containers
#include <set>           // std::set, std::multiset (conjunto ordenado, árvore balanceada)
#include <map>           // std::map, std::multimap (dicionário ordenado, árvore balanceada)
#include <span>          // Visão não-dona de sequência contígua (C++20+)
#include <tuple>         // Estruturas heterogêneas fixas
#include <algorithm>     // Funções genéricas: sort, copy, fill, etc.
#include <iterator> //funções e utilitários para trabalhar com iteradores.
#include <utility>   // std::pair, std::move, std::swap

// ---------------------------------------------------------------------------
// Utilitários gerais
// ---------------------------------------------------------------------------
#include <cstdint>   // Tipos inteiros fixos (uint8_t, int32_t, etc.)
#include <limits>    // Limites numéricos: std::numeric_limits<T>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>     // Funções matemáticas em std:: (sqrt, sin, cos, pow, etc.)
#include <iostream>  // Entrada/saída padrão: std::cout, std::cin
#include <string>    // std::string
#include <iomanip>   // Manipuladores de IO (std::setw, std::setprecision, etc.)
#include <numeric>   // Algoritmos numéricos (std::accumulate, std::inner_product)
#include <stdexcept> // Exceções padrão (std::runtime_error, std::invalid_argument)
#include <sstream>  // std::istringstream, std::ostringstream, std::stringstream (streams baseados em string)
#include <numbers>  // Constantes matemáticas (C++20+): std::numbers::pi, e, phi, sqrt2, etc.


// ---------------------------------------------------------------------------
// Memória, funções e metaprogramação
// ---------------------------------------------------------------------------
#include <memory>       // Ponteiros inteligentes (shared_ptr, unique_ptr, weak_ptr)
#include <variant>      // std::variant (union segura com tipo discriminado)
#include <optional>     // std::optional (valor opcional)
#include <functional>   // std::function, std::bind, std::hash
#include <type_traits>  // Traits e SFINAE (std::is_same, std::enable_if, etc.)

// ---------------------------------------------------------------------------
// mmcfilters: includes comuns a todo o projeto
// ---------------------------------------------------------------------------
#include "Image.hpp"
#include "../dataStructure/FastStack.hpp"
#include "../dataStructure/FastQueue.hpp"


namespace mmcfilters {

//habilitar ou desabilitar logs/debugs.
constexpr bool PRINT_LOG   = false;
constexpr bool PRINT_DEBUG = false;

//tipo de dado NodeId
using NodeId = int; //não usar unsigned int
constexpr NodeId InvalidNode = -1; //-1 indica nó inválido
inline bool isValidNode(NodeId id) noexcept { return id != InvalidNode;}
inline bool isInvalid(NodeId id) noexcept { return id == InvalidNode; }


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
