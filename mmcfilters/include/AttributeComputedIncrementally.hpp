

#ifndef ATTRIBUTE_COMPUTED_INCREMENTALLY_H
#define ATTRIBUTE_COMPUTED_INCREMENTALLY_H

#include "../include/NodeMT.hpp"
#include "../include/ImageUtils.hpp"
#include "../include/MorphologicalTree.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits> // Para usar std::numeric_limits<float>::epsilon()
#include <unordered_map>

#define PI 3.14159265358979323846


class AttributeNames {
public:
    
	std::unordered_map<std::string, int> mapIndexes;
	const int NUM_ATTRIBUTES;

    AttributeNames(int n) : NUM_ATTRIBUTES(29) {
        mapIndexes = {
            {"AREA", 0 * n},
            {"VOLUME", 1 * n},
            {"LEVEL", 2 * n},
            {"MEAN_LEVEL", 3 * n},
            {"VARIANCE_LEVEL", 4 * n},
            {"STANDARD_DEVIATION", 5 * n},
            {"BOX_WIDTH", 6 * n},
            {"BOX_HEIGHT", 7 * n},
            {"RECTANGULARITY", 8 * n},
            {"RATIO_WH", 9 * n},
            {"CENTRAL_MOMENT_20", 10 * n},
            {"CENTRAL_MOMENT_02", 11 * n},
            {"CENTRAL_MOMENT_11", 12 * n},
            {"CENTRAL_MOMENT_30", 13 * n},
            {"CENTRAL_MOMENT_03", 14 * n},
            {"CENTRAL_MOMENT_21", 15 * n},
            {"CENTRAL_MOMENT_12", 16 * n},
            {"ORIENTATION", 17 * n},
            {"LENGTH_MAJOR_AXIS", 18 * n},
            {"LENGTH_MINOR_AXIS", 19 * n},
            {"ECCENTRICITY", 20 * n},
            {"COMPACTNESS", 21 * n},
            {"HU_MOMENT_1_INERTIA", 22 * n},
            {"HU_MOMENT_2", 23 * n},
            {"HU_MOMENT_3", 24 * n},
            {"HU_MOMENT_4", 25 * n},
            {"HU_MOMENT_5", 26 * n},
            {"HU_MOMENT_6", 27 * n},
            {"HU_MOMENT_7", 28 * n}
        };
    }

};



class AttributeComputedIncrementally{

public:

 
    virtual void preProcessing(NodeMTPtr v);

    virtual void mergeChildren(NodeMTPtr parent, NodeMTPtr child);

    virtual void postProcessing(NodeMTPtr parent);

    void computerAttribute(NodeMTPtr root);

	static void computerAttribute(NodeMTPtr root, 
										std::function<void(NodeMTPtr)> preProcessing,
										std::function<void(NodeMTPtr, NodeMTPtr)> mergeChildren,
										std::function<void(NodeMTPtr)> postProcessing ){
		
		preProcessing(root);
			
		for(NodeMTPtr child: root->getChildren()){
			AttributeComputedIncrementally::computerAttribute(child, preProcessing, mergeChildren, postProcessing);
			mergeChildren(root, child);
		}

		postProcessing(root);
	}


	static std::vector<std::unordered_set<int>> extractCountors(MorphologicalTreePtr tree){
		std::vector<std::unordered_set<int>> contours(tree->getNumNodes());
		std::vector<std::list<int>> contoursToRemove(tree->getNumNodes());
		std::vector<int> ncount(tree->getNumRowsOfImage() * tree->getNumColsOfImage(), 0);
		AdjacencyRelationPtr adj4 = std::make_shared<AdjacencyRelation>(tree->getNumRowsOfImage(), tree->getNumColsOfImage(), 1);
		LCAEulerRMQ lca(tree);	

		AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
			[](NodeMTPtr node) -> void { // pre-processing

			},
			[&contours, &ncount, tree, adj4](NodeMTPtr parent, NodeMTPtr child) -> void { // merge-processing
				std::unordered_set<int> &Ncontour = contours[parent->getIndex()];
				for (int p : contours[child->getIndex()]){
					Ncontour.insert(p);
				}
			},
			[&contours, &contoursToRemove, &lca, &ncount, tree, adj4](NodeMTPtr node) -> void { // post-processing
				// Initialise contours of node "N"
				std::unordered_set<int> &Ncontour = contours[node->getIndex()];
				std::list<int> &NcontourToRemove = contoursToRemove[node->getIndex()];
				for(int p: NcontourToRemove){ //pixels que sao contornos de nodes descendentes ao NodeAtual
					bool isPixelToBeRemoved = true;
					for (int q : adj4->getNeighboringPixels(p)) { //Existe um nodeQ ascendente de NodeAtual contendo p como contorno? (p, q) in A
						NodeMTPtr nodeQ = tree->getSC(q); 
						if (tree->isStrictAncestor(nodeQ, node) ) { 
							contoursToRemove[nodeQ->getIndex()].push_back(p); 
							isPixelToBeRemoved = false;	
					  	}
					}
					if(!adj4->isBorderDomainImage(p) && isPixelToBeRemoved){
						Ncontour.erase(p);
					}
				}
			
				for (int p : node->getCNPs()) {
					if (adj4->isBorderDomainImage(p)){
						ncount[p]++;
					}
					for (int q : adj4->getNeighboringPixels(p)) {
						NodeMTPtr nodeQ = tree->getSC(q); 
						if(!tree->isComparable(node, tree->getSC(q))){ //se os nodeP e nodeQ não sao comparaveis, então p pode ser removido pelo LCA de nodeP e nodeQ 
							NodeMTPtr nodeLCA = lca.findLowestCommonAncestor(node, nodeQ);
							std::list<int> &NcontourToRemove = contoursToRemove[nodeLCA->getIndex()];
							NcontourToRemove.push_back(p);
							ncount[p]++;
						}
						else if(tree->isStrictDescendant(node, nodeQ)){  //maxtree:  SC(p) \subset SC(q) <=> f(p) > f(q)
					  		ncount[p]++;
						}else if (tree->isStrictAncestor(node, nodeQ)) { ////maxtree:  SC(q) \subset SC(p) <=> f(p) < f(q)
					  		ncount[q]--;
							if (ncount[q] == 0) {
								Ncontour.erase(q);
							}
						}
				  	}

				  	if (ncount[p] > 0){
						Ncontour.insert(p);
					}
				}

			}
		);
				  
		return contours;
	  }


	static float* computerAttribute(MorphologicalTreePtr tree, std::string attrName){
		const int n = tree->getNumNodes();
		float *attr = new float[n];
		auto [attributeNames, ptrValues] = AttributeComputedIncrementally::computerBasicAttributes(tree);
		int index = attributeNames.mapIndexes[attrName];
		for(int i = 0; i < n; i++){
			attr[i] = ptrValues[i + index];
		}
		delete[] ptrValues;
		return attr;
	}

	static float* computerStructTreeAttributes(MorphologicalTreePtr tree){
		const int numAttribute = 10;
		const int n = tree->getNumNodes();
		float *attrs = new float[n * numAttribute];
		/*
		0 - Altura do node
		1 - Profundidade do node
		2 - é node folha
		3 - é node root
		4 - Número de filhos do node
		5 - Número de irmãos do node
		6 - Número de node folha do node
		7 - Número de descendentes do node
		8 - Número de antecessores do node
		9 - Número de descendentes folha do node
		
		*/
		return attrs;
	}

	static std::pair<AttributeNames, float*> computerBasicAttributes(MorphologicalTreePtr tree){
	    
		/*
		0 - area
		1 - volume
		2 - level
		3 - mean level
		4 - variance level
		5 - standard deviation
		6 - Box width
		7 - Box height
		8 - rectangularity
		9 - ratio (Box width, Box height)
		
		10 - momentos centrais 20
		11 - momentos centrais 02
		12 - momentos centrais 11
		13 - momentos centrais 30
		14 - momentos centrais 03
		15 - momentos centrais 21
		16 - momentos centrais 12
		17 - orientation
		18 - lenght major axis
		19 - lenght minor axis
		20 - eccentricity = alongation
		21 - compactness = circularity
		22 - momentos de Hu 1 => inertia
		23 - momentos de Hu 2
		24 - momentos de Hu 3
		25 - momentos de Hu 4
		26 - momentos de Hu 5
		27 - momentos de Hu 6
		28 - momentos de Hu 7
		*/
		int n = tree->getNumNodes();
		AttributeNames attributeNames(n);
		
		float *attrs = new float[n * attributeNames.NUM_ATTRIBUTES];
		std::unordered_map<std::string, int> ATTR = attributeNames.mapIndexes;
		

		std::unique_ptr<int[]> xmax(new int[n]);
		std::unique_ptr<int[]> ymax(new int[n]);
		std::unique_ptr<int[]> xmin(new int[n]);
		std::unique_ptr<int[]> ymin(new int[n]);
		
		//momentos geometricos para calcular o centroide
		std::unique_ptr<long int[]> sumX(new long int[n]);//sum x
		std::unique_ptr<long int[]> sumY(new long int[n]);//sum y
		

		std::unique_ptr<long[]> sumGrayLevelSquare(new long[n]);
		int numCols = tree->getNumColsOfImage();
		int numRows = tree->getNumRowsOfImage();
		

		//computação dos atributos: area, volume, gray level, mean of gray level, variance of gray level, standard deviation gray level, Box width, Box height, rectangularity, ratio (Box width, Box height) e momentos geometricos 
	    AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
						[&ATTR, &attrs, n,  &xmax, &ymax, &xmin, &ymin, &sumX, &sumY, &sumGrayLevelSquare, numCols, numRows](NodeMTPtr node) -> void {
							attrs[node->getIndex() + ATTR["AREA"]  ] = node->getCNPs().size(); //area
							attrs[node->getIndex() + ATTR["VOLUME"]] = node->getCNPs().size() * node->getLevel(); //volume =>  \sum{ f }
							attrs[node->getIndex() + ATTR["LEVEL"] ] = node->getLevel(); //level

							xmax[node->getIndex()] = 0;
							ymax[node->getIndex()] = 0;
							xmin[node->getIndex()] = numCols;
							ymin[node->getIndex()] = numRows;
							sumX[node->getIndex()] = 0;
							sumY[node->getIndex()] = 0;
							sumGrayLevelSquare[node->getIndex()] = std::pow(node->getLevel(), 2) * node->getCNPs().size(); //computando: \sum{ f^2 }
							for(int p: node->getCNPs()) {
								int x = p % numCols;
								int y = p / numCols;
								xmin[node->getIndex()] = std::min(xmin[node->getIndex()], x);
								ymin[node->getIndex()] = std::min(ymin[node->getIndex()], y);
								xmax[node->getIndex()] = std::max(xmax[node->getIndex()], x);
								ymax[node->getIndex()] = std::max(ymax[node->getIndex()], y);

								sumX[node->getIndex()] += x;
								sumY[node->getIndex()] += y;
							}
						},
						[&ATTR, &attrs, n, &xmax, &ymax, &xmin, &ymin, &sumX, &sumY, &sumGrayLevelSquare](NodeMTPtr parent, NodeMTPtr child) -> void {
							attrs[parent->getIndex() + ATTR["AREA"]  ] += attrs[child->getIndex()]; //area
							attrs[parent->getIndex() + ATTR["VOLUME"]] += attrs[child->getIndex() + n]; //volume
							
							sumGrayLevelSquare[parent->getIndex()] += sumGrayLevelSquare[child->getIndex()]; //computando: \sum{ f^2 }

							ymax[parent->getIndex()] = std::max(ymax[parent->getIndex()], ymax[child->getIndex()]);
							xmax[parent->getIndex()] = std::max(xmax[parent->getIndex()], xmax[child->getIndex()]);
							ymin[parent->getIndex()] = std::min(ymin[parent->getIndex()], ymin[child->getIndex()]);
							xmin[parent->getIndex()] = std::min(xmin[parent->getIndex()], xmin[child->getIndex()]);
		
							sumX[parent->getIndex()] += sumX[child->getIndex()];
							sumY[parent->getIndex()] += sumY[child->getIndex()];
							
						},
						[&ATTR, &attrs, n, &xmax, &ymax, &xmin, &ymin, &sumGrayLevelSquare](NodeMTPtr node) -> void {
							
							float area = attrs[node->getIndex() + ATTR["AREA"]];
							float volume = attrs[node->getIndex() + ATTR["VOLUME"]]; 
							float width = xmax[node->getIndex()] - xmin[node->getIndex()] + 1;	
							float height = ymax[node->getIndex()] - ymin[node->getIndex()] + 1;	
							
							float meanGrayLevel = volume / area; //mean graylevel - // E(f)
							double meanGrayLevelSquare = sumGrayLevelSquare[node->getIndex()] / area; // E(f^2)
							float var = meanGrayLevelSquare - (meanGrayLevel * meanGrayLevel); //variance: E(f^2) - E(f)^2
							attrs[node->getIndex() + ATTR["VARIANCE_LEVEL"] ] = var > 0? var : 0; //variance
							
							if (attrs[node->getIndex() + ATTR["VARIANCE_LEVEL"]] >= 0) {
								attrs[node->getIndex() + ATTR["STANDARD_DEVIATION"]] = std::sqrt(attrs[node->getIndex() + ATTR["VARIANCE_LEVEL"]]); // desvio padrão do graylevel
							} else {
								attrs[node->getIndex() + ATTR["STANDARD_DEVIATION"]] = 0.0; // Se a variância for negativa, definir desvio padrão como 0
							}
							
							attrs[node->getIndex() + ATTR["MEAN_LEVEL"] ] = meanGrayLevel;
							attrs[node->getIndex() + ATTR["BOX_WIDTH"] ] = width;
							attrs[node->getIndex() + ATTR["BOX_HEIGHT"] ] = height;
							attrs[node->getIndex() + ATTR["RECTANGULARITY"] ] = area / (width * height);
							attrs[node->getIndex() + ATTR["RATIO_WH"] ] = std::max(width, height) / std::min(width, height);
		});

		

		//Computação dos momentos centrais e momentos de Hu
		AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
			[&ATTR, &attrs, n,  &sumX, &sumY, numCols](NodeMTPtr node) -> void {				
				attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_20"]] = 0; // momento central 20
				attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_02"]] = 0; // momento central 02
				attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_11"]] = 0; // momento central 11
				attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_30"]] = 0; // momento central 30
				attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_03"]] = 0; // momento central 03
				attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_21"]] = 0; // momento central 21
				attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_12"]] = 0; // momento central 12

				float xCentroid = sumX[node->getIndex()] / attrs[node->getIndex() + ATTR["AREA"]];
				float yCentroid = sumY[node->getIndex()] / attrs[node->getIndex() + ATTR["AREA"]];		
				for(int p: node->getCNPs()) {
					int x = p % numCols;
					int y = p / numCols;
					float dx = x - xCentroid;
            		float dy = y - yCentroid;
					
					// Momentos centrais de segunda ordem
					attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_20"] ] += std::pow(dx, 2);
					attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_02"] ] += std::pow(dy, 2);
					attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_11"] ] += dx * dy;
	
					// Momentos centrais de terceira ordem
					attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_30"]] += std::pow(dx, 3);
					attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_03"]] += std::pow(dy, 3);
					attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_21"]] += std::pow(dx, 2) * dy;
					attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_12"]] += dx * std::pow(dy, 2);
				}

			},
			[&ATTR, &attrs, n](NodeMTPtr parent, NodeMTPtr child) -> void {
				attrs[parent->getIndex() + ATTR["CENTRAL_MOMENT_20"]] += attrs[child->getIndex() + ATTR["CENTRAL_MOMENT_20"]];
				attrs[parent->getIndex() + ATTR["CENTRAL_MOMENT_02"]] += attrs[child->getIndex() + ATTR["CENTRAL_MOMENT_02"]];
				attrs[parent->getIndex() + ATTR["CENTRAL_MOMENT_11"]] += attrs[child->getIndex() + ATTR["CENTRAL_MOMENT_11"]];
				attrs[parent->getIndex() + ATTR["CENTRAL_MOMENT_30"]] += attrs[child->getIndex() + ATTR["CENTRAL_MOMENT_30"]];
				attrs[parent->getIndex() + ATTR["CENTRAL_MOMENT_03"]] += attrs[child->getIndex() + ATTR["CENTRAL_MOMENT_03"]];
				attrs[parent->getIndex() + ATTR["CENTRAL_MOMENT_21"]] += attrs[child->getIndex() + ATTR["CENTRAL_MOMENT_21"]];
				attrs[parent->getIndex() + ATTR["CENTRAL_MOMENT_12"]] += attrs[child->getIndex() + ATTR["CENTRAL_MOMENT_12"]];
			
			},
			[&ATTR, &attrs, n](NodeMTPtr node) -> void {
				
				float area = attrs[node->getIndex() + ATTR["AREA"]]; // area
				auto normMoment = [area](float moment, int p, int q){ 
					return moment / std::pow( area, (p + q + 2.0) / 2.0); 
				}; //função para normalizacao dos momentos				
				

				//Momentos centrais
				float mu20 = attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_20"]];
				float mu02 = attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_02"]];
				float mu11 = attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_11"]];
				float mu30 = attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_30"]];
				float mu03 = attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_03"]];
				float mu21 = attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_21"]];
				float mu12 = attrs[node->getIndex() + ATTR["CENTRAL_MOMENT_12"]];
					
				float discriminant = std::pow(mu20 - mu02, 2) + 4 * std::pow(mu11, 2);
					
				// Verificar se o denominador é zero antes de calcular atan2 para evitar divisão por zero
				if (mu20 != mu02 || mu11 != 0) {
					float radians = 0.5 * std::atan2(2 * mu11, mu20 - mu02);// orientação em radianos
					float degrees = radians * (180.0 / M_PI); // Converter para graus
					if (degrees < 0) { // Ajustar para o intervalo [0, 360] graus
						degrees += 360.0;
					}
					attrs[node->getIndex() + ATTR["ORIENTATION"]] = degrees; // Armazenar a orientação em graus no intervalo [0, 360]
				} else {
					attrs[node->getIndex() + ATTR["ORIENTATION"]] = 0.0; // Se não for possível calcular a orientação, definir um valor padrão
				}

				// Verificar se o discriminante é positivo para evitar raiz quadrada de números negativos
				if (discriminant < 0) {
					std::cerr << "Erro: Discriminante negativo, ajustando para zero." << std::endl;
					discriminant = 0;
				}	
				float a1 = mu20 + mu02 + std::sqrt(discriminant); // autovalores (correspondente ao eixo maior)
				float a2 = mu20 + mu02 - std::sqrt(discriminant); // autovalores (correspondente ao eixo menor)

				// Verificar se a1 e a2 são positivos antes de calcular sqrt para evitar NaN
				if (a1 > 0) {
					attrs[node->getIndex() + ATTR["LENGTH_MAJOR_AXIS"]] = std::sqrt((2 * a1) / area); // length major axis
				} else {
					attrs[node->getIndex() + ATTR["LENGTH_MAJOR_AXIS"]] = 0.0; // Definir valor padrão
				}

				if (a2 > 0) {
					attrs[node->getIndex() + ATTR["LENGTH_MINOR_AXIS"]] = std::sqrt((2 * a2) / area); // length minor axis
				} else {
					attrs[node->getIndex() + ATTR["LENGTH_MINOR_AXIS"]] = 0.0; // Definir valor padrão
				}

				// Verificar se a2 é diferente de zero antes de calcular a excentricidade
				attrs[node->getIndex() + ATTR["ECCENTRICITY"]] = (std::abs(a2) > std::numeric_limits<float>::epsilon()) ? a1 / a2 : a1 / 0.1; // eccentricity

				// Verificar se moment20 + mu02 é diferente de zero antes de calcular a compacidade
				if ((mu20 + mu02) > std::numeric_limits<float>::epsilon()) {
					attrs[node->getIndex() + ATTR["COMPACTNESS"]] = (1.0 / (2 * PI)) * (area / (mu20 + mu02)); // compactness
				} else {
					attrs[node->getIndex() + ATTR["COMPACTNESS"]] = 0.0; // Definir valor padrão
				}


				// Calcular os momentos normalizados
				float eta20 = normMoment(mu20, 2, 0);
				float eta02 = normMoment(mu02, 0, 2);
				float eta11 = normMoment(mu11, 1, 1);
				float eta30 = normMoment(mu30, 3, 0);
				float eta03 = normMoment(mu03, 0, 3);
				float eta21 = normMoment(mu21, 2, 1);
				float eta12 = normMoment(mu12, 1, 2);

				// Cálculo dos momentos de Hu
				attrs[node->getIndex() + ATTR["HU_MOMENT_1_INERTIA"]] = eta20 + eta02; // primeiro momento de Hu => inertia
				attrs[node->getIndex() + ATTR["HU_MOMENT_2"]] = std::pow(eta20 - eta02, 2) + 4 * std::pow(eta11, 2);
				attrs[node->getIndex() + ATTR["HU_MOMENT_3"]] = std::pow(eta30 - 3 * eta12, 2) + std::pow(3 * eta21 - eta03, 2);
				attrs[node->getIndex() + ATTR["HU_MOMENT_4"]] = std::pow(eta30 + eta12, 2) + std::pow(eta21 + eta03, 2);
				
				attrs[node->getIndex() + ATTR["HU_MOMENT_5"]] = 
					(eta30 - 3 * eta12) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - 3 * std::pow(eta21 + eta03, 2)) +
					(3 * eta21 - eta03) * (eta21 + eta03) * (3 * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));
				
				attrs[node->getIndex() + ATTR["HU_MOMENT_6"]] = 
					(eta20 - eta02) * (std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2)) + 
					4 * eta11 * (eta30 + eta12) * (eta21 + eta03);
				
				attrs[node->getIndex() + ATTR["HU_MOMENT_7"]] = 
					(3 * eta21 - eta03) * (eta30 + eta12) * (std::pow(eta30 + eta12, 2) - 3 * std::pow(eta21 + eta03, 2)) -
					(eta30 - 3 * eta12) * (eta21 + eta03) * (3 * std::pow(eta30 + eta12, 2) - std::pow(eta21 + eta03, 2));

				
		});
		return std::make_pair(attributeNames, attrs);
    }



};

#endif 