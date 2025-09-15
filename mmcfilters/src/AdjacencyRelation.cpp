#include "../include/AdjacencyRelation.hpp"
#include <math.h>
#include <cmath> 
#include <stdexcept>
#define PI 3.14159265358979323846


AdjacencyRelation::AdjacencyRelation(int numRows, int numCols, double radius){
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
		dr[i] = sqrt((dx * dx) + (dy * dy));
		if (i != i0) {
			da[i] = (atan2(-dy, -dx) * 180.0 / PI);
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

int AdjacencyRelation::getSize(){
	return this->n;
}

double AdjacencyRelation::getRadius(){
	return this->radius;
}

/*
bool AdjacencyRelation::isAdjacent(int px, int py, int qx, int qy) {
	double distance = std::sqrt(std::pow(px - qx, 2) + std::pow(py - qy, 2));
    return (distance <= radius);
}
*/
inline bool AdjacencyRelation::isAdjacent(int px, int py, int qx, int qy) const noexcept {
    int dx = px - qx;
    int dy = py - qy;
    return double(dx)*dx + double(dy)*dy <= radius2;
}

inline bool AdjacencyRelation::isAdjacent(int p, int q) const noexcept {
    int py = p / numCols, px = p % numCols;
    int qy = q / numCols, qx = q % numCols;

    return isAdjacent(px, py, qx, qy);
}


int AdjacencyRelation::nextValid() {
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

AdjacencyRelation::IteratorAdjacency AdjacencyRelation::begin() { 
    return IteratorAdjacency(this, nextValid()); 
}

AdjacencyRelation::IteratorAdjacency AdjacencyRelation::end() { 
    return IteratorAdjacency(this, this->n); 
}

AdjacencyRelation& AdjacencyRelation::getAdjPixels(int row, int col){
	if (row < 0 || row >= this->numRows || col < 0 || col >= this->numCols) {
        throw std::out_of_range("Índice fora dos limites.");
    }
    this->row = row;
    this->col = col;
    this->id = -1;
	this->forwardOnly = false;

    return *this;
}

AdjacencyRelation& AdjacencyRelation::getAdjPixels(int indexVector){
    return getAdjPixels(indexVector / this->numCols, indexVector % this->numCols);
}

AdjacencyRelation& AdjacencyRelation::getNeighborPixels(int row, int col){
	if (row < 0 || row >= this->numRows || col < 0 || col >= this->numCols) {
        throw std::out_of_range("Índice fora dos limites.");
    }
    this->row = row;
    this->col = col;
    this->id = 0;
	this->forwardOnly = false;
    return *this;
}

AdjacencyRelation& AdjacencyRelation::getNeighborPixels(int indexVector){
    return getNeighborPixels(indexVector / this->numCols, indexVector % this->numCols);
}

AdjacencyRelation& AdjacencyRelation::getNeighborPixelsForward(int row, int col){
	if (row < 0 || row >= this->numRows || col < 0 || col >= this->numCols) {
        throw std::out_of_range("Índice fora dos limites.");
    }
    this->row = row;
    this->col = col;
    this->id = 0;
	this->forwardOnly = true;
    return *this;

}

bool AdjacencyRelation::is4connectivity(){ return this->radius == 1;}
bool AdjacencyRelation::is8connectivity(){return this->radius == 1.5;}
    
bool AdjacencyRelation::isBorderDomainImage(int index){
	auto[row, col] = ImageUtils::to2D(index, this->numCols);
	return isBorderDomainImage(row, col);
}
bool AdjacencyRelation::isBorderDomainImage(int row, int col){
	return row == 0 || col == 0 || row == this->numRows - 1 || col == this->numCols - 1;
}

AdjacencyRelation& AdjacencyRelation::getNeighborPixelsForward(int indexVector){
	return getNeighborPixelsForward(indexVector / this->numCols, indexVector % this->numCols);
}

int AdjacencyRelation::getOffsetRow(int index){
	return offsetRow[index];
}
int AdjacencyRelation::getOffsetCol(int index){
	return offsetCol[index];
}