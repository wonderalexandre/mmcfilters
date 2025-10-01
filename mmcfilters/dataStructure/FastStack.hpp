#pragma once

#include <vector>   // usado como armazenamento dinâmico
#include <cstddef>  // para size_t (tamanho de container)
#include <utility>  // para std::move nas operações push/pop

namespace mmcfilters {
/**
 * @brief Pilha (stack) simples e performática baseada em `std::vector`.
 *
 * `FastStack<T>` provê a interface essencial de uma pilha LIFO com
 * operações de custo amortizado O(1) e controle de capacidade via
 * `reserve`. É útil em rotinas de DFS, processamento de componentes,
 * e estruturas auxiliares onde o overhead de `std::stack` e alocações
 * frequentes deve ser evitado.
 *
 * ## Operações
 * - `push(const T&)`, `push(T&&)` — insere no topo (amortizado O(1)).
 * - `pop()` — remove e retorna o topo (amortizado O(1)).
 * - `top()` — acesso ao elemento do topo (O(1)).
 * - `empty()`, `size()` — consultas O(1).
 * - `reserve(n)`, `clear()` — gestão de capacidade e limpeza.
 *
 * ## Exemplo de uso
 * @code
 * FastStack<int> st;
 * st.reserve(1024);
 * st.push(3);
 * st.push(7);
 * int x = st.top();   // 7
 * x = st.pop();       // 7; agora o topo é 3
 * @endcode
 */
template <typename T>
struct FastStack {
private:
    std::vector<T> data_;

public:
    FastStack() = default;

    explicit FastStack(size_t n) {
        data_.reserve(n);
    }

    /// Reserva espaço inicial (opcional)
    void reserve(size_t n) { data_.reserve(n); }

    /// Remove todos os elementos
    void clear() { data_.clear(); }

    /// Retorna se a pilha está vazia
    bool empty() const { return data_.empty(); }

    /// Retorna o tamanho atual da pilha
    size_t size() const { return data_.size(); }

    /// Adiciona um elemento ao topo
    void push(const T& value) { data_.push_back(value); }

    void push(T&& value) { data_.push_back(std::move(value)); }

    /// Remove e retorna o elemento do topo
    T pop() {
        T value = std::move(data_.back());
        data_.pop_back();
        return value;
    }

    /// Acesso ao topo sem remover
    T& top() { return data_.back(); }
    const T& top() const { return data_.back(); }
};
}