#ifndef COMMONS_HPP  
#define COMMONS_HPP  


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


#define PRINT_LOG 0
#define PRINT_DEBUG 0

using NodeId = int;
constexpr NodeId InvalidNode = -1;  // ou std::numeric_limits<NodeId>::max()
inline bool isValidNode(NodeId id) noexcept {
    return id != InvalidNode;
}


/**
 * @brief Classe de imagem genérica 2D com armazenamento contíguo e controle de vida via std::shared_ptr.
 *
 * A `Image<PixelType>` representa uma imagem 2D em ordem row-major, encapsulando
 * largura, altura e um buffer contíguo gerenciado por `std::shared_ptr<PixelType[]>`.
 * Fornece utilitários para criação/copiar/preencher e acesso indexado em 1D.
 *
 * ## Semântica de propriedade do buffer
 * - `create(rows, cols)`: aloca um novo buffer e o gerencia (deleter padrão).
 * - `create(rows, cols, initValue)`: idem, preenchendo com o valor inicial.
 * - `fromExternal(rawPtr, rows, cols)`: **não** assume a propriedade
 *   (o deleter é vazio). Útil quando o ciclo de vida do ponteiro bruto é externo.
 * - `fromRaw(rawPtr, rows, cols)`: **assume** a propriedade do array
 *   (deleter padrão de array). Use quando a instância deve gerenciar a memória.
 *
 * ## Layout de memória
 * - Acesso linear (row-major): índice `i = row * numCols + col`.
 * - Operador `operator[](int)` fornece acesso por índice linear.
 *
 *
 * ## Exemplo de uso
 * @code
 * using ImgU8 = Image<uint8_t>;
 * auto img = ImgU8::create(480, 640, 0);     // aloca e zera
 * img->fill(255);                            // preenche com 255
 * int idx = ImageUtils::to1D(10, 20, img->getNumCols());
 * (*img)[idx] = 128;                         // acesso linear
 * auto clone = img->clone();                  // deep copy
 * bool eq = img->isEqual(clone);              // true
 * @endcode
 *
 * @tparam PixelType tipo do pixel armazenado (ex.: uint8_t, int32_t, float).
 */
template <typename PixelType>
class Image {
    private:
        int numRows;
        int numCols;
        std::shared_ptr<PixelType[]> data;
        using Ptr = std::shared_ptr<Image<PixelType>>;
        
    public:
    using Type = PixelType;

    Image(int rows, int cols): numRows(rows), numCols(cols), data(new PixelType[rows * cols], std::default_delete<PixelType[]>()) {}

    static Ptr create(int rows, int cols) {
        return std::make_shared<Image>(rows, cols);
    }

    static Ptr create(int rows, int cols, PixelType initValue) {
        auto img = create(rows, cols);
        img->fill(initValue);
        return img;
    }

    static Ptr fromExternal(PixelType* rawPtr, int rows, int cols) {
        auto img = create(rows, cols);
        img->data = std::shared_ptr<PixelType[]>(rawPtr, [](PixelType*) {
            // deleter vazio: não libera o ponteiro
        });
        return img;
    }

    static Ptr fromRaw(PixelType* rawPtr, int rows, int cols) {
        auto img = create(rows, cols);
        img->data = std::shared_ptr<PixelType[]>(rawPtr, std::default_delete<PixelType[]>());
        return img;
    }

    
    void fill(PixelType value) {
        std::fill_n(data.get(), numRows * numCols, value);
    }

    bool isEqual(const Ptr& other) const {
        if (numRows != other->numRows || numCols != other->numCols)
            return false;
        int n = numRows * numCols;
        for (int i = 0; i < n; ++i) {
            if (data[i] != (*other)[i])
                return false;
        }
        return true;
    }
    
    Ptr clone() const {
        auto newImg = create(numRows, numCols);
        std::copy(data.get(), data.get() + (numRows * numCols), newImg->data.get());
        return newImg;
    }

    std::shared_ptr<PixelType[]> rawDataPtr(){ return data; }
    PixelType* rawData() { return data.get(); }
    int getNumRows() const { return numRows; }
    int getNumCols() const { return numCols; }
    int getSize() const { return numRows * numCols; }
    PixelType& operator[](int index) { return data[index]; }
    const PixelType& operator[](int index) const { return data[index]; }


};

// Aliases
using ImageUInt8 = Image<uint8_t>;
using ImageInt32 = Image<int32_t>;
using ImageFloat = Image<float>;

using ImageUInt8Ptr = std::shared_ptr<ImageUInt8>;
using ImageInt32Ptr = std::shared_ptr<ImageInt32>;
using ImageFloatPtr = std::shared_ptr<ImageFloat>;

template <typename PixelType>
using ImagePtr = std::shared_ptr<Image<PixelType>>;



class ImageUtils{
public:
    // Converte (row, col) para índice 1D (row-major)
    inline static int to1D(int row, int col, int numCols) noexcept{
        return row * numCols + col;
    }

    // Converte índice 1D para (row, col) (row-major)
    inline static std::pair<int, int> to2D(int index, int numCols) noexcept {
        int row = index / numCols;
        int col = index - row * numCols;  // evita operador % => int col = index % numCols;
        return {row, col};
    }

    // Cria uma imagem colorida aleatória a partir de uma imagem em escala de cinza: [(R,G,B), (R,G,B), ...]
    static ImageUInt8Ptr createRandomColor(int* img, int numRowsOfImage, int numColsOfImage){
        int max = 0;
        int sizeImage = numColsOfImage * numRowsOfImage;
        for (int i = 0; i < sizeImage; i++){
            if (img[i] > max)
                max = img[i];
        }

        std::unique_ptr<int[]> r(new int[max + 1]);
        std::unique_ptr<int[]> g(new int[max + 1]);
        std::unique_ptr<int[]> b(new int[max + 1]);
        r[0] = 0;
        g[0] = 0;
        r[0] = 0;
        for (int i = 1; i <= max; i++){
            r[i] = rand() % 256;
            g[i] = rand() % 256;
            b[i] = rand() % 256;
        }
        
        int sizeOutput = sizeImage * 3; // [(R,G,B), (R,G,B), ...]
        ImageUInt8Ptr outImage = ImageUInt8::create(numRowsOfImage, numColsOfImage * 3);
        
        auto output = outImage->rawData();
            // Inicializa com zero
        std::fill_n(output, sizeOutput, 0);

        for (int pidx = 0; pidx < sizeImage; pidx++){
            int cpidx = pidx * 3; // (coloured) for 3 channels
            output[cpidx]     = r[img[pidx]];
            output[cpidx + 1] = g[img[pidx]];
            output[cpidx + 2] = b[img[pidx]];
        }
        return outImage;
    }


};



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
     * // 4) Usando um range custom (ex.: RepsOfCCRange de um NodeCT)
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



#endif 
