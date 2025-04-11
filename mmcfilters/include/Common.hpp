#ifndef COMMONS_HPP  
#define COMMONS_HPP  


#define NDEBUG  // Remove os asserts do código
#include <cassert>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <memory>
#include <limits>

#define PRINT_LOG 1 
#define PRINT_DEBUG 0 

// Forward declaration 
class MorphologicalTree;
class NodeMT;


using NodeMTPtr = std::shared_ptr<NodeMT>;
using MorphologicalTreePtr = std::shared_ptr<MorphologicalTree>; 




#endif 
