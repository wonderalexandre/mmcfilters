

#ifndef ATTRIBUTE_COMPUTED_INCREMENTALLY_H
#define ATTRIBUTE_COMPUTED_INCREMENTALLY_H

#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include <algorithm> 
#include <cmath>
#include <iostream>
#include <limits> // Para usar std::numeric_limits<float>::epsilon()

#define PI 3.14159265358979323846

class AttributeComputedIncrementally{

public:
    virtual void preProcessing(NodeCT *v);

    virtual void mergeChildren(NodeCT *parent, NodeCT *child);

    virtual void postProcessing(NodeCT *parent);

    void computerAttribute(NodeCT *root);

	static void computerAttribute(NodeCT* root, 
										std::function<void(NodeCT*)> preProcessing,
										std::function<void(NodeCT*, NodeCT*)> mergeChildren,
										std::function<void(NodeCT*)> postProcessing ){
		
		preProcessing(root);
			
		for(NodeCT *child: root->getChildren()){
			AttributeComputedIncrementally::computerAttribute(child, preProcessing, mergeChildren, postProcessing);
			mergeChildren(root, child);
		}

		postProcessing(root);
	}

	static float* computerAttribute(ComponentTree *tree, int indexAttr){
		const int n = tree->getNumNodes();
		float *attr = new float[n];
		float* attributes = AttributeComputedIncrementally::computerBasicAttributes(tree);
		for(int i = 0; i < n; i++){
			attr[i] = attributes[i + (n * indexAttr)];
		}
		delete[] attributes;
		return attr;
	}

	static float* computerStructTreeAttributes(ComponentTree *tree){
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

    static float* computerBasicAttributes(ComponentTree *tree){
	    const int numAttribute = 19;
		const int n = tree->getNumNodes();
		float *attrs = new float[n * numAttribute];
		/*
		0 - area
		1 - volume
		2 - level
		3 - mean level
		4 - variance level
		5 - Box width
		6 - Box height
		7 - rectangularity
		8 - ratio (Box width, Box height)
		9 - momentos centrais 20
		10 - momentos centrais 02
		11 - momentos centrais 11
		12 - momentos de Hu 1 => inertia
		13 - orientation
		14 - lenght major axis
		15 - lenght minor axis
		16 - eccentricity = alongation
		17 - compactness = circularity
		18 - standard deviation
		*/
		
		int xmax[n]; //min value
		int ymax[n]; //min value
		int xmin[n]; //max value
		int ymin[n]; //max value
		long int m10[n];
		long int m01[n];
		long int m20[n];
		long int m02[n];
		long int m11[n];
		long sumGrayLevelSquare[n];

		int numCols = tree->getNumColsOfImage();
		int numRows = tree->getNumRowsOfImage();
		

	    AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
						[&attrs, n,  &xmax, &ymax, &xmin, &ymin, &m10, &m01,&m11, &m20, &m02, &sumGrayLevelSquare, numCols, numRows](NodeCT* node) -> void {
							attrs[node->getIndex()         ] = node->getCNPs().size(); //area
							attrs[node->getIndex() + n     ] = node->getCNPs().size() * node->getLevel(); //volume =>  \sum{ f }
							attrs[node->getIndex() + (n*2) ] = node->getLevel(); //level

							xmax[node->getIndex()] = 0;
							ymax[node->getIndex()] = 0;
							xmin[node->getIndex()] = numCols;
							ymin[node->getIndex()] = numRows;
							m10[node->getIndex()] = 0;
							m11[node->getIndex()] = 0;
							m01[node->getIndex()] = 0;
							m20[node->getIndex()] = 0;
							m02[node->getIndex()] = 0;
							sumGrayLevelSquare[node->getIndex()] = std::pow(node->getLevel(), 2) * node->getCNPs().size(); //computando: \sum{ f^2 }
							for(int p: node->getCNPs()) {
								int x = p % numCols;
								int y = p / numCols;
								xmin[node->getIndex()] = std::min(xmin[node->getIndex()], x);
								ymin[node->getIndex()] = std::min(ymin[node->getIndex()], y);
								xmax[node->getIndex()] = std::max(xmax[node->getIndex()], x);
								ymax[node->getIndex()] = std::max(ymax[node->getIndex()], y);
								m10[node->getIndex()] += x;
								m01[node->getIndex()] += y;
								m11[node->getIndex()] += x*y;
								m20[node->getIndex()] += x*x;
								m02[node->getIndex()] += y*y;								
							}


						},
						[&attrs, n, &xmax, &ymax, &xmin, &ymin, &m10, &m01,&m11, &m20, &m02, &sumGrayLevelSquare](NodeCT* parent, NodeCT* child) -> void {
							attrs[parent->getIndex()       ] += attrs[child->getIndex()]; //area
							attrs[parent->getIndex() + n   ] += attrs[child->getIndex() + n]; //volume
							
							sumGrayLevelSquare[parent->getIndex()] += sumGrayLevelSquare[child->getIndex()]; //computando: \sum{ f^2 }

							ymax[parent->getIndex()] = std::max(ymax[parent->getIndex()], ymax[child->getIndex()]);
							xmax[parent->getIndex()] = std::max(xmax[parent->getIndex()], xmax[child->getIndex()]);
							ymin[parent->getIndex()] = std::min(ymin[parent->getIndex()], ymin[child->getIndex()]);
							xmin[parent->getIndex()] = std::min(xmin[parent->getIndex()], xmin[child->getIndex()]);
		
							m10[parent->getIndex()] += m10[child->getIndex()];
							m01[parent->getIndex()] += m01[child->getIndex()];
							m11[parent->getIndex()] += m11[child->getIndex()];
							m20[parent->getIndex()] += m20[child->getIndex()];
							m02[parent->getIndex()] += m02[child->getIndex()];

						},
						[&attrs, n, &xmax, &ymax, &xmin, &ymin, &m10, &m01, &m11, &m20, &m02, &sumGrayLevelSquare](NodeCT* node) -> void {
							
							float area = attrs[node->getIndex()];
							float volume = attrs[node->getIndex() + n]; 
							float width = xmax[node->getIndex()] - xmin[node->getIndex()] + 1;	
							float height = ymax[node->getIndex()] - ymin[node->getIndex()] + 1;	
							
							attrs[node->getIndex() + (n*3) ] = volume / area; //mean graylevel
							float meanGrayLevel = attrs[node->getIndex() + (n*3) ]; // E(f) 
							double meanGrayLevelSquare = sumGrayLevelSquare[node->getIndex()] / area; // E(f^2)
							attrs[node->getIndex() + (n*4) ] = meanGrayLevelSquare - (meanGrayLevel * meanGrayLevel); //variance: E(f^2) - E(f)^2
							if (attrs[node->getIndex() + (n*4)] >= 0) {
								attrs[node->getIndex() + (n*18)] = std::sqrt(attrs[node->getIndex() + (n*4)]); // desvio padrão do graylevel
							} else {
								attrs[node->getIndex() + (n*18)] = 0.0; // Se a variância for negativa, definir desvio padrão como 0
							}

							attrs[node->getIndex() + (n*5) ] = width;
							attrs[node->getIndex() + (n*6) ] = height;
							attrs[node->getIndex() + (n*7) ] = area / (width * height);
							attrs[node->getIndex() + (n*8) ] = std::max(width, height) / std::min(width, height);
		});

		
		AttributeComputedIncrementally::computerAttribute(tree->getRoot(),
			[&attrs, n,  &m10, &m01, &m11, &m20, &m02, numCols](NodeCT* node) -> void {
				
				attrs[node->getIndex() + (n*9) ] = 0; //central_moments20
				attrs[node->getIndex() + (n*10) ] = 0; //central_moment02
				attrs[node->getIndex() + (n*11) ] = 0; //central_moment11
				
				float xCentroid = m10[node->getIndex()] / attrs[node->getIndex()];
				float yCentroid = m01[node->getIndex()] / attrs[node->getIndex()];		
				for(int p: node->getCNPs()) {
					int x = p % numCols;
					int y = p / numCols;
					attrs[node->getIndex() + (n*9) ] += std::pow(x - xCentroid, 2);
					attrs[node->getIndex() + (n*10) ] += std::pow(y - yCentroid, 2);
					attrs[node->getIndex() + (n*11) ] += (x - xCentroid) * (y - yCentroid);
				}
			},
			[&attrs, n, &m10, &m01, &m11, &m20, &m02](NodeCT* parent, NodeCT* child) -> void {
				attrs[parent->getIndex() + (n*9) ] += attrs[child->getIndex() + (n*9) ];
				attrs[parent->getIndex() + (n*10) ] += attrs[child->getIndex() + (n*10) ];
				attrs[parent->getIndex() + (n*4) ] += attrs[child->getIndex() + (n*4) ];
			},
			[&attrs, n, &m10, &m01, &m11, &m20, &m02](NodeCT* node) -> void {
				
				
				//Momentos centrais
				float moment20 = attrs[node->getIndex() + (n*9)];  // moment20
				float moment02 = attrs[node->getIndex() + (n*10)]; // moment02
				float moment11 = attrs[node->getIndex() + (n*11)]; // moment11
				float area = attrs[node->getIndex()]; // area
				
				auto normMoment = [area](float moment, int p, int q){ 
					return moment / std::pow( area, (p + q + 2.0) / 2.0); 
				}; //função para normalizacao dos momentos


				attrs[node->getIndex() + (n*12)] = normMoment(moment02, 0, 2) + normMoment(moment20, 2, 0); // primeiro momento de Hu => inertia
					
				float discriminant = std::pow(moment20 - moment02, 2) + 4 * std::pow(moment11, 2);
					
				// Verificar se o denominador é zero antes de calcular atan2 para evitar divisão por zero
				if (moment20 != moment02 || moment11 != 0) {
					
					float radians = 0.5 * std::atan2(2 * moment11, moment20 - moment02);// orientação em radianos
					float degrees = radians * (180.0 / M_PI); // Converter para graus
					if (degrees < 0) { // Ajustar para o intervalo [0, 360] graus
						degrees += 360.0;
					}
					attrs[node->getIndex() + (n*13)] = degrees; // Armazenar a orientação em graus no intervalo [0, 360]
				} else {
					attrs[node->getIndex() + (n*13)] = 0.0; // Se não for possível calcular a orientação, definir um valor padrão
				}

				// Verificar se o discriminante é positivo para evitar raiz quadrada de números negativos
				if (discriminant < 0) {
					std::cerr << "Erro: Discriminante negativo, ajustando para zero." << std::endl;
					discriminant = 0;
				}	
				float a1 = moment20 + moment02 + std::sqrt(discriminant); // autovalores (correspondente ao eixo maior)
				float a2 = moment20 + moment02 - std::sqrt(discriminant); // autovalores (correspondente ao eixo menor)

				// Verificar se a1 e a2 são positivos antes de calcular sqrt para evitar NaN
				if (a1 > 0) {
					attrs[node->getIndex() + (n*14)] = std::sqrt((2 * a1) / area); // length major axis
				} else {
					attrs[node->getIndex() + (n*14)] = 0.0; // Definir valor padrão
				}

				if (a2 > 0) {
					attrs[node->getIndex() + (n*15)] = std::sqrt((2 * a2) / area); // length minor axis
				} else {
					attrs[node->getIndex() + (n*15)] = 0.0; // Definir valor padrão
				}

				// Verificar se a2 é diferente de zero antes de calcular a excentricidade
				attrs[node->getIndex() + (n*16)] = (std::abs(a2) > std::numeric_limits<float>::epsilon()) ? a1 / a2 : a1 / 0.1; // eccentricity

				// Verificar se moment20 + moment02 é diferente de zero antes de calcular a compacidade
				if ((moment20 + moment02) > std::numeric_limits<float>::epsilon()) {
					attrs[node->getIndex() + (n*17)] = (1.0 / (2 * PI)) * (area / (moment20 + moment02)); // compactness
				} else {
					attrs[node->getIndex() + (n*17)] = 0.0; // Definir valor padrão
				}
				
		});
		return attrs;
    }



};

#endif 