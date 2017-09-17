#include <iostream>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string>

#include "Globals.h"
#include "TerrainGenerator.h"

int main()
{
	int seed = 2;							// startseed for random numbers
	int algorithm;						// number of chosen algorithm
	float roughness = 0.1f;		// roughness of terrain
	float fir = 0.65f;				// value for the FIR filter to blur the Heightmap
	int size;									// size of Heightmap
	std::string filename;			// filename
	// Fault formation
	int iterations=16;
	int filterIterations=8;
	// Particle Deposition
	int mountain=20;
	int moveDrop=10;
	int particle=100000;
	float caldera=7.0f;
	// Perlin Noise
	double persistence;
	double frequency;
	double amplitude;
	int octaves;
	// Voronoi
	int points;

	// Erosion
	bool erosion = false;					// Erosion?
	int erosion_type;							// type of Erosion Algorithm
	bool erosion_debug = false;		// save Heightmap with and without Erosion
	int erosion_iterations;

	// tick counter
	clock_t start, end;

	std::cout << "Heightmap Generator - Werner Ziegelwanger 2011" << std::endl << std::endl;
	// choose algorithm
	do{
		std::cout << "Algorithms:" << std::endl;
		std::cout << "Fraktale Algorithmen" << std::endl;
		std::cout << "1. Triangle Division" << std::endl;
		std::cout << "2. Diamond-Square" << std::endl;
		std::cout << "3. Midpoint Displacement" << std::endl;
		std::cout << "Prozedurale Algorithmen" << std::endl;
		std::cout << "4. Fault Formation" << std::endl;
		std::cout << "5. Particle Deposition" << std::endl;
		std::cout << "6. Perlin Noise" << std::endl;
		std::cout << "7. Voronoi Diagram" << std::endl;
		std::cout << "- - - - - - - - - - - - -" << std::endl;
		std::cout << "0. Hilfe" << std::endl;
		std::cout << "--> ";
		std::cin >> algorithm;
	}while(algorithm<0 || algorithm>7);
	std::cout << std::endl;
	// set paramters
	switch(algorithm)
	{
	case 0:
		std::cout << "Hilfe" << std::endl;
		std::cout << "Nach der Auswahl eines Algorithmus werden einzelne Parameter für die Berechnung abgefragt. Werte in [] sind Defaultwerte." << std::endl;
		break;
	case 1: //Triangle Division
		std::cout << "Triangle Division" << std::endl;
		std::cout << "roughness [0.7]: --> ";
		std::cin >> roughness;
		std::cout << "FIR value [0.65]: --> ";
		std::cin >> fir;
		break;
	case 2: //Diamond Square
		std::cout << "Diamond Square" << std::endl;
		std::cout << "roughness [0.9]: --> ";
		std::cin >> roughness;
		std::cout << "FIR value [0.65]: --> ";
		std::cin >> fir;
		break;
	case 3: // Midpoint Displacement
		std::cout << "Midpoint Displacement" << std::endl;
		std::cout << "roughness [0.5]: --> ";
		std::cin >> roughness;
		std::cout << "FIR value [0.65]: --> ";
		std::cin >> fir;
		break;
	case 4: // Fault Formation
		std::cout << "Fault Formation" << std::endl;
		std::cout << "iterations [128]: --> ";
		std::cin >> iterations;
		std::cout << "filter iterations [8]: --> ";
		std::cin >> filterIterations;
		std::cout << "FIR value [0.4]: --> ";
		std::cin >> fir;
		break;
	case 5: // Particle Deposition
		std::cout << "Particle Deposition" << std::endl;
		std::cout << "mountains [20]: --> ";
		std::cin >> mountain;
		std::cout << "move drop position [10]: --> ";
		std::cin >> moveDrop;
		std::cout << "particles [100000]: --> ";
		std::cin >> particle;
		std::cout << "caldera height [7.0]: --> ";
		std::cin >> caldera;
		std::cout << "FIR value [0.65]: --> ";
		std::cin >> fir;
		break;
	case 6: // Perlin Noise
		std::cout << "Perlin Noise" << std::endl;
		std::cout << "octaves [5]: --> ";
		std::cin >> octaves;
		std::cout << "amplitude [0.5]: --> ";
		std::cin >> amplitude;
		std::cout << "frequency [0.02]: --> ";
		std::cin >> frequency;
		std::cout << "persistence [1.0]: --> ";
		std::cin >> persistence;
		std::cout << "FIR value [0.65]: --> ";
		std::cin >> fir;
		break;
	case 7: // Voronoi Diagram
		std::cout << "Voronoi Diagram" << std::endl;
		std::cout << "points [20]: --> ";
		std::cin >> points;
		std::cout << "FIR value [0.65]: --> ";
		std::cin >> fir;
		break;
	}

	// set seed
	std::cout << "seed [0]: --> ";
	std::cin >> seed;
	std::srand(seed);
	std::cout << "size [1024]: --> ";
	std::cin >> size;
	std::cout << "filename [heightmap]: --> ";
	std::cin >> filename;

	// erosion
	char c;
	std::cout << "erosion [j/n]: --> ";
	std::cin >> c;
	if(c == 'j' || c == 'J')
	{
		erosion = true;
		std::cout << "Erosion Types:" << std::endl;
		std::cout << "1. Thermal Erosion" << std::endl;
		std::cout << "2. Hydraulic Erosion" << std::endl;
		std::cout << "3. Both" << std::endl;
		std::cout << "--> ";
		std::cin >> erosion_type;
		if(erosion_type < 1 || erosion_type > 3)
			erosion_type = 3;
		std::cout << "iterations [100]: --> ";
		std::cin >> erosion_iterations;
		std::cout << "debug [j/n]: --> ";
		std::cin >> c;
		if(c == 'j' || c == 'J')
			erosion_debug = true;
	}

	// generate Terrain
	TerrainGenerator tg(size);
	switch(algorithm)
	{
	case 1: 
		start = clock();
		tg.makeTriangleDivision(roughness, seed, fir); 
		end = clock();
		std::cout << "Runtime (TriangleDivision): " << start << " to " << end << " --> " << end-start << " / " << CLOCKS_PER_SEC << " seconds" << std::endl;
		break;
	case 2: 
		start = clock();
		tg.makeDiamondSquare(roughness, seed, fir);
		end = clock();
		std::cout << "Runtime (DiamondSquare): " << start << " to " << end << " --> " << end-start << " / " << CLOCKS_PER_SEC << " seconds" << std::endl;
		break;
	case 3: 
		start = clock();
		tg.makeMidpointDisplacement(roughness, seed, fir);
		end = clock();
		std::cout << "Runtime (MidpointDisplacement): " << start << " to " << end << " --> " << end-start << " / " << CLOCKS_PER_SEC << " seconds" << std::endl;
		break;
	case 4: 
		start = clock();
		tg.makeFaultFormation(iterations, filterIterations, fir);
		end = clock();
		std::cout << "Runtime (FaultFormation): " << start << " to " << end << " --> " << end-start << " / " << CLOCKS_PER_SEC << " seconds" << std::endl;
		break;
	case 5: 
		start = clock();
		tg.makeParticleDeposition(mountain, moveDrop, particle, caldera, fir);
		end = clock();
		std::cout << "Runtime (ParticleDeposition): " << start << " to " << end << " --> " << end-start << " / " << CLOCKS_PER_SEC << " seconds" << std::endl;
		break;
	case 6: 
		start = clock();
		tg.makePerlinNoise(persistence, frequency, amplitude, octaves, seed, fir);
		end = clock();
		std::cout << "Runtime (PerlinNoise): " << start << " to " << end << " --> " << end-start << " / " << CLOCKS_PER_SEC << " seconds" << std::endl;
		break;
	case 7: 
		start = clock();
		tg.makeVoronoiDiagram(points, seed, fir);
		end = clock();
		std::cout << "Runtime (VoronoiDiagram): " << start << " to " << end << " --> " << end-start << " / " << CLOCKS_PER_SEC << " seconds" << std::endl;
		break;
	}
	
	// Erosion
	if(erosion)
	{
		// Save Heightmap without Erosion if Debug is on
		if(erosion_debug)
		{
			std::string f = filename;
			std::string file = f.append("_NoErosion");
			tg.writeToFile(file);
		}
		// Compute talus and compute Erosion
		float talus = static_cast<float>(4.0f/size);
		
		switch(erosion_type)
		{
		case 1: 
			start = clock();
			tg.makeThermalErosion(talus, erosion_iterations);
			end = clock();
			std::cout << "Runtime (ThermalErosion): " << start << " to " << end << " --> " << end-start << " / " << CLOCKS_PER_SEC << " seconds" << std::endl;
			break;
		case 2: 
			start = clock();
			tg.makeHydraulicErosion(0.1f, 0.1f, 0.5f, 0.6f, erosion_iterations);
			end = clock();
			std::cout << "Runtime (HydraulicErosion): " << start << " to " << end << " --> " << end-start << " / " << CLOCKS_PER_SEC << " seconds" << std::endl;
			break;
		case 3: 
			start = clock();
			tg.makeThermalErosion(talus, erosion_iterations*0.5);
			tg.makeHydraulicErosion(0.1f, 0.1f, 0.5f, 0.6f, erosion_iterations*0.5);
			end = clock();
			std::cout << "Runtime (Thermal and HydraulicErosion): " << start << " to " << end << " --> " << end-start << " / " << CLOCKS_PER_SEC << " seconds" << std::endl;
			break;
		}
	}

	// make and write Game Maps
	tg.makeAccessabilityMap(filename, 0.005f);
	tg.makeFlatnessMap(filename, 0.0025f);

	// generate a Texture
	std::string f = filename;
	std::string file = f.append("_Texture");
	tg.makeProceduralTexture(file);

	// write Heightmap to File
	start = clock();
	tg.writeToFile(filename);
	end = clock();
	std::cout << "Runtime (SaveHeightmap): " << start << " to " << end << " --> " << end-start << " / " << CLOCKS_PER_SEC << " seconds" << std::endl;

	std::cout << "Done." << std::endl;
}