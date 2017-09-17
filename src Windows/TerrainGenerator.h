#ifndef _TERRAINGENERATOR_
#define _TERRAINGENERATOR_

#include "Globals.h"
#include "EasyBMP/EasyBMP.h"
#include <string>

class TerrainGenerator
{
public:
	// Functions
	TerrainGenerator(int size);
	~TerrainGenerator();

	void writeToFile(std::string filename);
	void writeToFile(std::string filename, float *data);

	// Fractal Algorithms
	void makeTriangleDivision(float roughness, int seed, float firValue);
	void makeDiamondSquare(float roughness, int seed, float firValue);
	void makeFourierTransformation(float roughness);

	// Procedural Algorithms
	void makeFaultFormation(int iterations, int filterIterations, float firValue);
	void makeMidpointDisplacement(float roughness, int seed, float firValue);
	void makeParticleDeposition(int nMountain, int moveDrop, int particle, float caldera, float firValue);
	void makePerlinNoise(float persistence, float frequency, float amplitude, int octaves, int seed, float firValue);
	void makeVoronoiDiagram(int points, int seed, float firValue);

	// Erosion Algorithms
	void makeThermalErosion(float talus, int iterations);
	void makeHydraulicErosion(float water, float sediment, float evaporation, float capacity, int iterations);

	// Game Maps
	void makeAccessabilityMap(std::string filename, float threshold);
	void makeFlatnessMap(std::string filename, float threshold);

	//Texture
	void wrapAround(int &i, int &j, int max_i, int max_j);
	void makeProceduralTexture(std::string filename);

private:
	// Functions
	inline int computeIndex(int a, int b);
	float random(float a, float b);
	int irand( int a, int b);
	void normalize();
	void filterTerrain(float filter);
	void filterFIR(float *value, int ofs, float filter);

	void FFT2D(int n, int m, bool inverse, double *gRe, double *gIm, double *GRe, double *GIm);
	void createCalderas(int x, int y, float height);						// Particle Deposition

	void createGameMap(float threshold, float *data);

	//Variables
	int m_iSize;				// size of the terrain in pixels
	int m_iSize2;				// pow of size
	int m_iSizeMask;		// Size-1
	float *m_fTerrain;	// array with heightmap of terrain
};

#endif