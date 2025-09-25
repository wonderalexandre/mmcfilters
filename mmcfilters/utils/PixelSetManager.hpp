#pragma once

#include "Common.hpp"

namespace mmcfilters {


/**
 * @brief Gerenciador de conjuntos disjuntos de pixels (flat zones ou cnps) com listas circulares e mapeamentos O(1).
 *
 * O `PixelSetManager` mantém a relação entre pixels e seus conjuntos (flat zones ou cnps)
 * usando quatro vetores paralelos: `pixelToIndex`, `indexToPixel`, `sizeSets` e
 * `pixelsNext`. O desenho provê operações O(1) para consultas e splices
 * (concatenação de listas circulares) durante fusões de conjuntos, além de
 * *views* baseadas em `std::span` para iteração sem cópias.
 *
 * ## Estruturas internas
 * - `pixelToIndex[p]` → índice do conjunto ao qual o pixel representante `p` pertence.
 * - `indexToPixel[i]` → representante (pixel head) do conjunto de índice `i`.
 * - `sizeSets[i]` → tamanho (número de pixels) do conjunto `i`.
 * - `pixelsNext[p]` → próximo pixel na lista circular do conjunto ao qual `p` pertence.
 *
 * ## Operações principais
 * - `numSets()`, `numPixelsInSet(rep)`, `numPixelsInSets(reps)` — consultas O(1)/O(k).
 * - `mergeSetsByRep(repWinner, repLoser)` — fusão O(1) com splice de listas circulares.
 * - `shrinkToNumSets(n)` — reduz vetores de conjuntos ao número real de FZs ou |numNodes|
 * - *Views*: `view()`, `viewOf*()` expõem `std::span` para zero-cópia.
 * - Iteração de pixels por set: `getPixelsBySet(...)` retorna um range lazy.
 * - Iteração de representantes válidos: `getFlatzoneRepresentatives()` retorna um range.
 *
 * ## Complexidade
 * - Acesso por mapeamento: O(1).
 * - Fusão (`mergeSetsByRep`): O(1) (atualiza contadores e faz splice das listas).
 * - Iteração: O(|S|) proporcional ao número de pixels/sets percorridos.
 *
 * ## Invariantes
 * - Se `indexToPixel[i] == -1` então o slot do conjunto `i` é inválido (após fusões).
 * - Para qualquer representante `rep`: `pixelToIndex[rep]` aponta para um índice `i`
 *   tal que `indexToPixel[i] == headRep` e `sizeSets[i] > 0`.
 * - As listas de pixels de um mesmo conjunto formam um ciclo por `pixelsNext`.
 *
 * ## Exemplo mínimo
 * @code
 * PixelSetManager psm(numPixels);
 * // inicialização dos reps/estruturas omitida
 * int a = repA, b = repB;
 * psm.mergeSetsByRep(a, b);       // b fundido em a
 * for (int px : psm.getPixelsBySet(a)) {
 *   // processa pixels do conjunto resultante
 * }
 * for (int rep : psm.getFlatzoneRepresentatives()) {
 *   // percorre reps válidos
 * }
 * @endcode
 */
struct PixelSetManager{
    
    std::vector<int> pixelToIndex; //mapeamento do pixel representante para índice na lista de conjuntos disjuntos. Tamanho: numPixels
    std::vector<int> indexToPixel; // mapeamento de índice para pixel representante. Tamanho: numSets
    std::vector<int> sizeSets; // usada para armazenar o tamanho dos conjuntos disjuntos. Tamanho: numSets
    std::vector<int> pixelsNext; // mapa de pixels dos conjuntos disjuntos: Tamanho: numPixels

    PixelSetManager(int numPixels, int numSets)
        : pixelToIndex(numPixels, -1), indexToPixel(numSets, -1), sizeSets(numSets, 0), pixelsNext(numPixels, -1) { }
        
    PixelSetManager(int numPixels)
        : pixelToIndex(numPixels, -1), indexToPixel(numPixels, -1), sizeSets(numPixels, 0), pixelsNext(numPixels, -1) { }
    
    int numSets() const { return sizeSets.size(); }

    int numPixelsInSet(int rep){ return sizeSets[pixelToIndex[rep]]; }

    int numPixelsInSets(const std::vector<int>& reps){
        int sum = 0;
        for (int rep : reps) {
            sum += sizeSets[pixelToIndex[rep]];
        }
        return sum;
    }
    
    int indexOfPixel(int pixel) const {
        return pixelToIndex[pixel];
    }
    
    int pixelOfIndex(int idx) const {
        return indexToPixel[idx];
    }

    /**
     * @brief Redimensiona os vetores relacionados a conjuntos (flat zones)
     * para refletir o número real de FZs criadas.
     *
     * @param newNumSets Número real de conjuntos encontrados.
     */
    void shrinkToNumSets(int newNumSets) {
        indexToPixel.resize(newNumSets);
        sizeSets.resize(newNumSets);
    }
    
    void mergeSetsByRep(int repWinner, int repLoser) {
        
        // 1. Recupera índices dos representantes
        int idxRootWinner = pixelToIndex[repWinner];
        int idxRootLoser  = pixelToIndex[repLoser];
        sizeSets[idxRootWinner] += sizeSets[idxRootLoser];

        // 2. Splice O(1) das listas circulares (pixels)
        int nextWinner = pixelsNext[repWinner];
        int nextLoser  = pixelsNext[repLoser];
        pixelsNext[repWinner] = nextLoser;
        pixelsNext[repLoser]  = nextWinner;

        // 3. Invalida slot perdedor
        sizeSets[idxRootLoser]  = 0;
        indexToPixel[idxRootLoser] = -1;

        // 4. Redireciona lookups pelo antigo rep pixel
        pixelToIndex[repLoser] = idxRootWinner;
    }


    /**
     * @brief *View* de spans que expõe os vetores internos do gerenciador.
     */
    struct View {
        std::span<int> pixelToIndex;
        std::span<int> indexToPixel;
        std::span<int> sizeSets;
        std::span<int> pixelsNext;
    };

    View view() noexcept {
        return {std::span<int>(pixelToIndex), std::span<int>(indexToPixel), std::span<int>(sizeSets), std::span<int>(pixelsNext) };
    }

    
    std::span<int> viewOfPixelToIndex(){ return std::span<int>(pixelToIndex); }
    std::span<int> viewOfIndexToPixel(){ return std::span<int>(indexToPixel); }
    std::span<int> viewOfSizeSets(){ return std::span<int>(sizeSets); }
    std::span<int> viewOfPixelsNext(){ return std::span<int>(pixelsNext); }
    
    
    /**
     * @brief Faixa iterável sobre os pixels de uma ou mais sets.
     *
     * Esta classe encapsula um range de representantes de sets (ex.: `int`, 
     * `std::vector<int>`, `std::span<const int>`, `RepsOfCCRange` etc.) e fornece 
     * iteradores (`PixelsBySetIterator`) que percorrem todos os pixels desses sets. 
     *
     * O funcionamento é baseado nos seguintes princípios:
     *  - Cada representante identifica um set.
     *  - O iterador avança sobre todos os pixels de cada set usando a lista circular 
     *    interna (`pixelsNext`) e o tamanho registrado em `sizeSets`.
     *  - Ao término de um set, o iterador passa automaticamente para o próximo 
     *    set indicada no range de representantes.
     *  - O range é apenas uma "view" sobre os reps. Ele não copia os pixels, apenas 
     *    percorre dinamicamente as listas já armazenadas.
     *
     * Exemplos de uso:
     * @code
     *
     * // 1) Único representante
     * for (int px : v.getPixelsBySet(rep)) {
     *     // processa cada pixel do set de 'rep'
     * }
     *
     * // 2) Vários reps (std::vector<int>)
     * std::vector<int> reps = {rep1, rep2};
     * for (int px : v.getPixelsBySet(reps)) {
     *     // processa pixels dos sets de rep1 e rep2
     * }
     *
     * // 3) Usando span (sem cópia do vetor)
     * std::span<const int> s(reps.data(), reps.size());
     * for (int px : v.getPixelsBySet(s)) {
     *     // processa pixels, sem overhead de cópia
     * }
     *
     * // 4) Usando um range custom (ex.: RepsOfCCRange de um NodeMT)
     * for (int px : g.getPixelsBySet(node->getRepsOfCC())) {
     *     // processa pixels de todos os sets alcançadas na BFS
     * }
     * @endcode
     *
     * @tparam Range Tipo do container ou view que contém os reps.
     *               Pode ser `std::array<int,N>`, `std::vector<int>`, 
     *               `std::span<const int>`, ou um range custom compatível
     *               com `std::begin`/`std::end` que produza `int`.
     */
    template<class Range>
    class PixelsBySetRange {
    private:
        PixelSetManager::View   v_;  // guarda o View (spans) por valor
        Range reps_;   // range de representantes

        using RepIt = decltype(std::begin(std::declval<const Range&>()));

    public:
        // -------- Iterator ----------
        /**
         * @brief Iterador que percorre os pixels dos conjuntos fornecidos.
         */
        class PixelsBySetIterator {
            PixelSetManager::View   v_;
            RepIt it_, last_;
            int   cur_{-1};      // pixel atual
            int   remaining_{0}; // pixels restantes na FZ atual

            void startNextSegment() {
                cur_ = -1;
                remaining_ = 0;
                while (it_ != last_) {
                    const int rep = *it_;
                    const int idx = v_.pixelToIndex[rep];
                    if (idx >= 0) {
                        const int head = v_.indexToPixel[idx];
                        const int sz   = v_.sizeSets[idx];
                        if (head != -1 && sz > 0) {
                            cur_ = head;
                            remaining_ = sz;
                            return;
                        }
                    }
                    ++it_; // tenta próximo representante
                }
            }

        public:
            using iterator_category = std::forward_iterator_tag;
            using value_type        = int;
            using difference_type   = std::ptrdiff_t;
            using pointer           = const int*;
            using reference         = const int&;

            PixelsBySetIterator(PixelSetManager::View v, RepIt first, RepIt last) : v_(v), it_(first), last_(last) {
                startNextSegment();
            }

            reference operator*()  const { return cur_; }
            pointer   operator->() const { return &cur_; }

            PixelsBySetIterator& operator++() {
                if (remaining_ > 0) {
                    --remaining_;
                    if (remaining_ == 0) {
                        ++it_;
                        startNextSegment();        // próxima FZ
                    } else {
                        cur_ = v_.pixelsNext[cur_]; // próximo pixel na FZ
                    }
                }
                return *this;
            }

            PixelsBySetIterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }

            friend bool operator==(const PixelsBySetIterator& a, const PixelsBySetIterator& b) {
                return a.cur_ == b.cur_ && a.it_ == b.it_ && a.last_ == b.last_;
            }
            friend bool operator!=(const PixelsBySetIterator& a, const PixelsBySetIterator& b) {
                return !(a == b);
            }

        };

        // -------- Range ----------
        PixelsBySetRange(PixelSetManager::View v, Range reps)
            : v_(v), reps_(std::move(reps)) {}

        PixelsBySetIterator begin() const { return PixelsBySetIterator(v_, std::begin(reps_), std::end(reps_)); }
        PixelsBySetIterator end()   const { return PixelsBySetIterator(v_, std::end(reps_),  std::end(reps_));  }
    };

    //Um único representante
    auto getPixelsBySet(int rep) {
        return PixelsBySetRange<std::array<int,1>>(this->view(), std::array<int,1>{rep});
    }

    // Qualquer range (vector<int>, span<const int>, RepsOfCCRange, …)
    template<class Range>
    auto getPixelsBySet(Range reps) {
        return PixelsBySetRange<Range>(this->view(), std::move(reps));
    }




    /**
     * @brief Iterador para percorrer todos os representantes de sets válidos.
     *
     * Este iterador percorre o array `indexToPixel`, saltando entradas inválidas
     * (marcadas com -1). Cada elemento retornado é o pixel representante da set.
     *
     * Uso típico:
     * @code
     * for (int rep : getFlatzoneRepresentatives()) {
     *     // rep é um representante válido
     * }
     * @endcode
     */
    class RepresentativeIterator {
    private:
        std::span<int> indexToPixel_;
        size_t size_;
        size_t idx_;

        void skipInvalid() {
            while (idx_ < size_ && indexToPixel_[idx_] == -1) {
                ++idx_;
            }
        }

    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type        = int;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const int*;
        using reference         = const int&;

        RepresentativeIterator(std::span<int> data, size_t size, size_t startIdx)
            : indexToPixel_(data), size_(size), idx_(startIdx) {
            skipInvalid();
        }

        reference operator*()  const { return indexToPixel_[idx_]; }
        pointer   operator->() const { return &indexToPixel_[idx_]; }

        RepresentativeIterator& operator++() {
            ++idx_;
            skipInvalid();
            return *this;
        }

        RepresentativeIterator operator++(int) {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        friend bool operator==(const RepresentativeIterator& a, const RepresentativeIterator& b) {
            return a.idx_ == b.idx_;
        }
        friend bool operator!=(const RepresentativeIterator& a, const RepresentativeIterator& b) {
            return !(a == b);
        }
    };


    /**
     * @brief Faixa iterável de representantes de flat zones válidos.
     *
     * Retorna um objeto que pode ser usado em range-based for loops
     * para iterar apenas sobre os representantes ativos.
     */
    class RepresentativeRange {
    private:
        std::span<int> indexToPixel_;
        size_t size_;

    public:
        explicit RepresentativeRange(std::span<int> data, size_t size)
            : indexToPixel_(data), size_(size) {}

        RepresentativeIterator begin() const { return RepresentativeIterator(indexToPixel_, size_, 0); }
        RepresentativeIterator end()   const { return RepresentativeIterator(indexToPixel_, size_, size_); }
    };

    /**
     * @brief Obtém um range iterável de representantes de flat zones válidos.
     *
     * Uso:
     * @code
     * for (int rep : getFlatzoneRepresentatives()) {
     *     // processar rep
     * }
     * @endcode
     */
    RepresentativeRange getFlatzoneRepresentatives()  {
        return RepresentativeRange(viewOfIndexToPixel(), indexToPixel.size());
    }


};



}