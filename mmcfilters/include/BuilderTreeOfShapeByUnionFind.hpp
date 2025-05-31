#include <algorithm>
#include <climits>
#include <vector>
#include <utility>
#include <list>
#include <deque>

#include "../include/AdjacencyUC.hpp"
#include "../include/Common.hpp"


#ifndef BUILDER_TREE_OF_SHAPE_BY_UNION_FIND_H
#define BUILDER_TREE_OF_SHAPE_BY_UNION_FIND_H

class BuilderTreeOfShapeByUnionFind {
private:
    int interpNumRows;
    int interpNumCols;
    std::unique_ptr<uint8_t[]> interpolationMin;
    std::unique_ptr<uint8_t[]> interpolationMax;
    std::unique_ptr<uint8_t[]> imgU;
    std::unique_ptr<int[]> parent;
    std::unique_ptr<int[]> imgR; 
    
    std::unique_ptr<AdjacencyUC> adj;
    bool is4c8cConnectivity;

    class PriorityQueueToS {
    private:
        std::vector<std::deque<int>> buckets;
        int currentPriority;
        int numElements;
        int maxPriorityLevels;
        

    public:
        PriorityQueueToS(int depthOfImage=8) : currentPriority(0), numElements(0), maxPriorityLevels(1 << depthOfImage){
            buckets.resize(maxPriorityLevels);
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

                    // Tentar diminuir a prioridade
                    if (j > 0 && buckets[j].empty()) {
                        j--;
                    }
                    if (!buckets[j].empty()) { // Encontrou o próximo bucket não vazio diminuindo a prioridade
                        currentPriority = j;
                        break;
                    }

                    // Tentar aumentar a prioridade
                    if (i < maxPriorityLevels && buckets[i].empty()) {
                        i++;
                    }
                    if (i < maxPriorityLevels && !buckets[i].empty()) { // Encontrou o próximo bucket não vazio aumentando a prioridade
                        currentPriority = i;
                        break;
                    }
                }
            }

            int element = buckets[currentPriority].front(); 
            buckets[currentPriority].pop_front();           

            numElements--;  
            return element;
        }
    };
    

public:

    int getInterpNumRows();
    int getInterpNumCols();
    
    uint8_t* getImgU();
    int* getImgR();
    int* getParent();
    
    BuilderTreeOfShapeByUnionFind();
    ~BuilderTreeOfShapeByUnionFind();
    void interpolateImage(ImageUInt8Ptr img);
    void interpolateImage4c8c(ImageUInt8Ptr img);
    void sort();
    int findRoot(int zPar[], int x);
    void createTreeByUnionFind();


};

#endif