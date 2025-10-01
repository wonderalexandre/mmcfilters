#pragma once

#include "../utils/Common.hpp"


namespace mmcfilters {

class AdjacencyRelation;  // forward declaration
using AdjacencyRelationPtr = std::shared_ptr<AdjacencyRelation>;

/**
 * @brief Relação de adjacência em grade 2D com raio arbitrário e iteração eficiente.
 *
 * Define offsets de vizinhança para uma janela de raio real (ex.: 1.0 → 4-adj, 1.5 → 8-adj),
 * provendo utilitários para listar vizinhos de um pixel e iterar por eles com um iterador leve.
 * Também oferece uma variação "forward" que emite somente metade dos vizinhos (sem simetria),
 * útil para varreduras assimétricas e construção de arestas únicas.
 */
class AdjacencyRelation {
private:
    int id;
    
    int row;
    int col;    
    int numCols;
    int numRows;
    double radius;
    double radius2;
    int n;
    bool forwardOnly = false; //meia 

    std::vector<int> offsetRow;
    std::vector<int> offsetCol;
    std::vector<uint8_t> forwardMask;// máscara “forward” por offset i: true se (dy>0) || (dy==0 && dx>0)


public:
    /**
     * @brief Constrói uma relação de adjacência para imagem `numRows`×`numCols`.
     * @param numRows Número de linhas da imagem.
     * @param numCols Número de colunas da imagem.
     * @param radius Raio da vizinhança (1.0 ≈ 4-conexão, 1.5 ≈ 8-conexão).
     */    
    AdjacencyRelation(int numRows, int numCols, double radius){
        this->numRows = numRows;
        this->numCols = numCols;
        this->radius = radius;
        this->radius2 = radius * radius;

        int i, j, k, dx, dy, r0, r2, i0 = 0;
        this->n = 0;
        r0 = (int) radius;
        r2 = (int) radius2;
        for (dy = -r0; dy <= r0; dy++)
            for (dx = -r0; dx <= r0; dx++)
                if (((dx * dx) + (dy * dy)) <= r2)
                    this->n++;
        
        i = 0;
        this->offsetCol.resize(this->n);
        this->offsetRow.resize(this->n);
        
        for (dy = -r0; dy <= r0; dy++) {
            for (dx = -r0; dx <= r0; dx++) {
                if (((dx * dx) + (dy * dy)) <= r2) {
                    this->offsetCol[i] =dx;
                    this->offsetRow[i] =dy;
                    if ((dx == 0) && (dy == 0))
                        i0 = i;
                    i++;
                }
            }
        }
            
        float aux;
        std::vector<float> da(n);
        std::vector<float> dr(n);
        
        /* Set clockwise */
        for (i = 0; i < n; i++) {
            dx = this->offsetCol[i];
            dy = this->offsetRow[i];
            dr[i] = std::sqrt((dx * dx) + (dy * dy));
            if (i != i0) {
                da[i] = (std::atan2(-dy, -dx) * 180.0 / std::numbers::pi);
                if (da[i] < 0.0)
                    da[i] += 360.0;
            }
        }
        da[i0] = 0.0;
        dr[i0] = 0.0;

        /* place central pixel at first */
        aux = da[i0];
        da[i0] = da[0];
        da[0] = aux;

        aux = dr[i0];
        dr[i0] = dr[0];
        dr[0] = aux;

        int auxX, auxY;
        auxX = this->offsetCol[i0];
        auxY = this->offsetRow[i0];
        this->offsetCol[i0] = this->offsetCol[0];
        this->offsetRow[i0] = this->offsetRow[0];
            
        this->offsetCol[0] = auxX;
        this->offsetRow[0] = auxY;
            

        /* sort by angle */
        for (i = 1; i < n - 1; i++) {
            k = i;
            for (j = i + 1; j < n; j++)
                if (da[j] < da[k]) {
                    k = j;
                }
            aux = da[i];
            da[i] = da[k];
            da[k] = aux;
            aux = dr[i];
            dr[i] = dr[k];
            dr[k] = aux;

            auxX = this->offsetCol[i];
            auxY = this->offsetRow[i];
            this->offsetCol[i] = this->offsetCol[k];
            this->offsetRow[i] = this->offsetRow[k];
                
            this->offsetCol[k] = auxX;
            this->offsetRow[k] = auxY;
        }

        /* sort by radius for each angle */
        for (i = 1; i < n - 1; i++) {
            k = i;
            for (j = i + 1; j < n; j++)
                if ((dr[j] < dr[k]) && (da[j] == da[k])) {
                    k = j;
                }
            aux = dr[i];
            dr[i] = dr[k];
            dr[k] = aux;

            auxX = this->offsetCol[i];
            auxY = this->offsetRow[i];
            this->offsetCol[i] = this->offsetCol[k];
            this->offsetRow[i] = this->offsetRow[k];
                
            this->offsetCol[k] = auxX;
            this->offsetRow[k] = auxY;
                
        }

        // máscara forward: 1, se (dy>0) || (dy==0 && dx>0); caso contrario é 0
        forwardMask.resize(n, 0);
        for (int k = 1; k < n; ++k) {
            int dx = offsetCol[k], dy = offsetRow[k];
            forwardMask[k] = (dy > 0 || (dy == 0 && dx > 0)) ? 1 : 0;
        }
        forwardMask[0] = 0;
    }


    /**
     * @brief Avança para o próximo offset válido conforme os limites e máscara.
     * @return Índice do próximo offset válido, ou tamanho para fim.
     */
    int nextValid() {
        id += 1;
        while (id < n) {

            // checa "forward" se necessário
            if (forwardOnly && !forwardMask[id]) { id += 1; continue; }

            // coordenadas do vizinho
            const int newRow = row + offsetRow[id];
            const int newCol = col + offsetCol[id];

            if (newRow >= 0 && newRow < numRows && newCol >= 0 && newCol < numCols) {
                return id;
            }
            id += 1;
        }
        return n;
    }

    /**
     * @brief Retorna a quantidade de offsets no stencil atual.
     */
    int getSize(){
        return this->n;
    }

    /**
     * @brief Configura (row,col) e prepara iteração de adjacentes sem filtro forward. Esse método incluí a origem.
     */
    AdjacencyRelation& getAdjPixels(int row, int col){
        if (row < 0 || row >= this->numRows || col < 0 || col >= this->numCols) {
            throw std::out_of_range("Índice fora dos limites.");
        }
        this->row = row;
        this->col = col;
        this->id = -1;
        this->forwardOnly = false;

        return *this;
    }

    /**
     * @brief Configura por índice linear e prepara iteração de adjacentes sem filtro. Esse método incluí a origem.
     */
    AdjacencyRelation& getAdjPixels(int indexVector){
        return getAdjPixels(indexVector / this->numCols, indexVector % this->numCols);
    }
    

    /**
     * @brief Configura (row,col) e prepara iteração de vizinhos dentro dos limites. Esse método NÃO incluí a origem.
     */
        AdjacencyRelation& getNeighborPixels(int row, int col){
        if (row < 0 || row >= this->numRows || col < 0 || col >= this->numCols) {
            throw std::out_of_range("Índice fora dos limites.");
        }
        this->row = row;
        this->col = col;
        this->id = 0;
        this->forwardOnly = false;
        return *this;
    }

    /**
     * @brief Configura por índice linear e prepara iteração de vizinhos dentro dos limites. Esse método NÃO incluí a origem.
     */
    AdjacencyRelation& getNeighborPixels(int indexVector){
        return getNeighborPixels(indexVector / this->numCols, indexVector % this->numCols);
    }

    /**
     * @brief Prepara iteração apenas sobre metade dos vizinhos (forward-only) em (row,col).
     * @note Útil para gerar pares (p,q) sem duplicação (p→q, nunca q→p).
     */
    AdjacencyRelation& getNeighborPixelsForward(int row, int col){
        if (row < 0 || row >= this->numRows || col < 0 || col >= this->numCols) {
            throw std::out_of_range("Índice fora dos limites.");
        }
        this->row = row;
        this->col = col;
        this->id = 0;
        this->forwardOnly = true;
        return *this;
    }
    
    /**
     * @brief Versão por índice linear de `getNeighborPixelsForward`.
     */
    AdjacencyRelation& getNeighborPixelsForward(int indexVector){
        return getNeighborPixelsForward(indexVector / this->numCols, indexVector % this->numCols);
    }
    

    /**
     * @brief Verifica adjacência por índices lineares (p,q).
     */
    inline bool isAdjacent(int p, int q) const noexcept {
        int py = p / numCols, px = p % numCols;
        int qy = q / numCols, qx = q % numCols;

        return isAdjacent(px, py, qx, qy);
    }

    /**
     * @brief Verifica adjacência por coordenadas (px,py) e (qx,qy).
     */
    inline bool isAdjacent(int px, int py, int qx, int qy) const noexcept {
        int dx = px - qx;
        int dy = py - qy;
        return double(dx)*dx + double(dy)*dy <= radius2;
    }

    /**
     * @brief Retorna o raio em uso.
     */
    double getRadius(){
        return this->radius;
    }

    bool is4connectivity(){ return this->radius == 1;}
    bool is8connectivity(){return this->radius == 1.5;}
        
    bool isBorderDomainImage(int index){
        auto[row, col] = ImageUtils::to2D(index, this->numCols);
        return isBorderDomainImage(row, col);
    }
    bool isBorderDomainImage(int row, int col){
        return row == 0 || col == 0 || row == this->numRows - 1 || col == this->numCols - 1;
    }

    int getOffsetRow(int index){
        return offsetRow[index];
    }
    int getOffsetCol(int index){
        return offsetCol[index];
    }


    /**
     * @brief Iterador leve para percorrer vizinhos já configurados via `get*`.
     *
     * Produz índices lineares de pixels vizinhos válidos, respeitando os limites e,
     * quando configurado, a máscara forward-only.
     */
    class IteratorAdjacency { 
    private:
        int index;
        AdjacencyRelation* instance; 

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = int;

        IteratorAdjacency(AdjacencyRelation* obj, int id) :  index(id), instance(obj) { }

        AdjacencyRelation* getInstance() { return instance; } 

        IteratorAdjacency& operator++() { 
            this->index = instance->nextValid();  
            return *this; 
        }

        bool operator==(const IteratorAdjacency& other) const { 
            return index == other.index; 
        }
        bool operator!=(const IteratorAdjacency& other) const { 
            return !(*this == other);
        }

        int operator*() const { 
            return (instance->row + instance->offsetRow[index]) * instance->numCols + (instance->col + instance->offsetCol[index]); 
        }
    };
    /**
     * @brief Início da iteração de vizinhos conforme configuração atual.
     */
    IteratorAdjacency begin() { 
        return IteratorAdjacency(this, nextValid()); 
    }

    /**
     * @brief Marcador de fim da iteração de vizinhos.
     */
    IteratorAdjacency end() { 
        return IteratorAdjacency(this, this->n); 
    }
};

} // namespace mmcfilters



