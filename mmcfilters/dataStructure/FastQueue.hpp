#pragma once

#include <vector>   // container dinâmico usado como buffer
#include <cstddef>  // definição de size_t
#include <utility>  // std::move para otimizar push e pop

namespace mmcfilters {
/**
 * @brief Fila linear baseada em std::vector para alto desempenho em BFS.
 *
 * Essa estrutura encapsula um std::vector<T> e um índice de leitura (`head_`),
 * funcionando como uma fila FIFO.  
 * Ao contrário de std::queue, não tem overhead de alocações/push/pop, pois:
 *   - `push` adiciona ao fim do vetor.
 *   - `pop` apenas avança o índice `head_`.
 *   - `clear` reseta o índice e o tamanho para reutilização sem desalocar memória.
 *
 * Uso típico:
 * @code
 * FastQueue<int> q;
 * q.reserve(2048);
 * q.push(10);
 * q.push(20);
 * while (!q.empty()) {
 *     int x = q.pop();
 *     // processa x
 * }
 * @endcode
 */
template <typename T>
struct FastQueue {
private:
    std::vector<T> data_;
    size_t head_ = 0; // índice do próximo elemento a ser lido

public:
    FastQueue() = default;

    FastQueue(size_t n){
        data_.reserve(n); 
    } 

    /// Reserva espaço inicial (opcional, para evitar realocações)
    void reserve(size_t n) { data_.reserve(n); }

    /// Remove todos os elementos, reseta o índice de leitura
    void clear() { data_.clear(); head_ = 0; }

    /// Retorna se a fila está vazia
    bool empty() const { return head_ >= data_.size(); }

    /// Retorna o tamanho atual da fila
    size_t size() const { return data_.size() - head_; }

    /// Adiciona um elemento ao fim
    void push(const T& value) { data_.push_back(value); }

    void push(T&& value) { data_.push_back(std::move(value)); }

    /// Remove e retorna o próximo elemento
    T pop() { return std::move(data_[head_++]); }

    /// Acesso ao próximo elemento sem remover
    T& front() { return data_[head_]; }
    const T& front() const { return data_[head_]; }
};
}