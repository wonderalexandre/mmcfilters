#pragma once

#include <list>
#include <vector>
#include "../include/Common.hpp"

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
    AdjacencyRelation(int numRows, int numCols, double radius);
    /**
     * @brief Avança para o próximo offset válido conforme os limites e máscara.
     * @return Índice do próximo offset válido, ou tamanho para fim.
     */
    int nextValid();
    /**
     * @brief Retorna a quantidade de offsets no stencil atual.
     */
    int getSize();
    /**
     * @brief Configura (row,col) e prepara iteração de adjacentes sem filtro forward. Esse método incluí a origem.
     */
    AdjacencyRelation& getAdjPixels(int row, int col);
    /**
     * @brief Configura por índice linear e prepara iteração de adjacentes sem filtro. Esse método incluí a origem.
     */
    AdjacencyRelation& getAdjPixels(int index);
    /**
     * @brief Configura (row,col) e prepara iteração de vizinhos dentro dos limites. Esse método NÃO incluí a origem.
     */
    AdjacencyRelation& getNeighborPixels(int row, int col);
    /**
     * @brief Configura por índice linear e prepara iteração de vizinhos dentro dos limites. Esse método NÃO incluí a origem.
     */
    AdjacencyRelation& getNeighborPixels(int index);
    /**
     * @brief Prepara iteração apenas sobre metade dos vizinhos (forward-only) em (row,col).
     * @note Útil para gerar pares (p,q) sem duplicação (p→q, nunca q→p).
     */
    AdjacencyRelation& getNeighborPixelsForward(int row, int col);
    /**
     * @brief Versão por índice linear de `getNeighborPixelsForward`.
     */
    AdjacencyRelation& getNeighborPixelsForward(int index);

    /**
     * @brief Verifica adjacência por índices lineares (p,q).
     */
    inline bool isAdjacent(int p, int q)const noexcept;
    /**
     * @brief Verifica adjacência por coordenadas (px,py) e (qx,qy).
     */
    inline bool isAdjacent(int px, int py, int qx, int qy) const noexcept;
    /**
     * @brief Retorna o raio em uso.
     */
    double getRadius();

    int getOffsetRow(int index);
    int getOffsetCol(int index);
    
    bool is4connectivity();
    bool is8connectivity();
    
    bool isBorderDomainImage(int index);
    bool isBorderDomainImage(int row, int col);


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
    IteratorAdjacency begin();
    /**
     * @brief Marcador de fim da iteração de vizinhos.
     */
   IteratorAdjacency end();	 
};

