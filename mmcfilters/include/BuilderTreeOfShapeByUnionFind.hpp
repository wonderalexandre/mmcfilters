#include <algorithm>
#include <climits>
#include <vector>
#include <utility>
#include <array>
#include <list>

#include "../include/AdjacencyRelation.hpp"


#ifndef BUILDER_TREE_OF_SHAPE_BY_UNION_FIND_H
#define BUILDER_TREE_OF_SHAPE_BY_UNION_FIND_H

class BuilderTreeOfShapeByUnionFind {
private:
    int interpNumRows;
    int interpNumCols;
    int* interpolationMin;
    int* interpolationMax;
    int* parent;
    int* imgR; 
    int* imgU;
    AdjacencyRelation* adj;


    class PriorityQueueToS {
    private:
        std::array<std::vector<int>, 256> buckets; 
        int currentPriority;
        int numElements;

    public:
        PriorityQueueToS() : currentPriority(0), numElements(0) {
            for (int i = 0; i < 256; ++i) {
                buckets[i] = std::vector<int>(); // Inicializar cada bucket como um vetor vazio
            }
        }

        void initial(int element, int priority) {
            currentPriority = priority;
            buckets[priority].push_back(element);
            numElements++;
        }
        int getCurrentPriority()  {return currentPriority;}
        bool isEmpty()  {return numElements == 0;}

        void priorityPush(int element, int lower, int upper) {
            int priority;
            if (lower > currentPriority) {
                priority = lower;
            } else if (upper < currentPriority) {
                priority = upper;
            } else {
                priority = currentPriority;
            }
            numElements++;
            buckets[priority].push_back(element);
        }

        int priorityPop() {
            // Se o bucket atual estiver vazio, precisamos ajustar a prioridade
            if (buckets[currentPriority].empty()) {
                int i = currentPriority;
                int j = currentPriority;
                while (true) {
                    // Tentar aumentar a prioridade
                    if (i < 256 && buckets[i].empty()) {
                        i++;
                    }
                    if (i < 256 && !buckets[i].empty()) { // Encontrou o próximo bucket não vazio aumentando a prioridade
                        currentPriority = i;
                        break;
                    }
                    // Tentar diminuir a prioridade
                    if (j > 0 && buckets[j].empty()) {
                        j--;
                    }
                    if (!buckets[j].empty()) { // Encontrou o próximo bucket não vazio diminuindo a prioridade
                        currentPriority = j;
                        break;
                    }
                }
            }

            int element = buckets[currentPriority].at(buckets[currentPriority].size() - 1); // Mudança aqui!
            buckets[currentPriority].pop_back();
            numElements--;  
            return element;
        }
    };
    

public:

    int getInterpNumRows();
    int getInterpNumCols();
    int* getInterpolationMin();
    int* getInterpolationMax();
    int* getImgR();
    int* getImgU();
    int* getParent();

    BuilderTreeOfShapeByUnionFind();
    ~BuilderTreeOfShapeByUnionFind();
    void interpolateImage(int* img, int num_rows, int num_cols);
    void sort();
    int findRoot(int zPar[], int x);
    void createTreeByUnionFind();

};

#endif