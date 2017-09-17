#include "TerrainGenerator.h"
#include "Globals.h"
#include "PerlinNoise.h"
#include <time.h>
#include <vector>
#include <algorithm>

// point struct for voronoi diagram
struct vpoint
{
	int x;
	int y;
};

// --------------------------------------------------------------------------------------------------------------------------------
TerrainGenerator::TerrainGenerator(int size)
{
	m_iSize = size;																				// get size of terrain
	m_iSize2 = m_iSize * m_iSize;
	m_iSizeMask = m_iSize - 1;
	m_fTerrain = new float[ m_iSize2 ];										// create new array for heightmap
	memset( m_fTerrain, 0, m_iSize2 * sizeof( float ) );	// fill the array with 0
}

// --------------------------------------------------------------------------------------------------------------------------------
TerrainGenerator::~TerrainGenerator()
{
	delete[] m_fTerrain;      // safe delete heightmap
	m_fTerrain = 0;
}

// --------------------------------------------------------------------------------------------------------------------------------
void TerrainGenerator::writeToFile(std::string filename)
{
	// generate filename with filetype
	std::string file = filename.append(".bmp");

	//save BMP
	BMP oFile;
	oFile.SetBitDepth(8);
	oFile.SetSize(m_iSize,m_iSize);
	CreateGrayscaleColorTable(oFile);
	for(int i = 0; i < m_iSize; i++)
	{
		for(int j = 0; j < m_iSize; j++)
		{
			unsigned char o = (unsigned char)( m_fTerrain[ i + j * m_iSize ] * MAX_HEIGHT );
			oFile(j,i)->Red = o;
			oFile(j,i)->Green = o;
			oFile(j,i)->Blue = o;
		}
	}

	oFile.WriteToFile(file.c_str());

} // writeToFile

// --------------------------------------------------------------------------------------------------------------------------------
void TerrainGenerator::writeToFile(std::string filename, float *data)
{
	// generate filename with filetype
	std::string file = filename.append(".bmp");

	//save BMP
	BMP oFile;
	oFile.SetBitDepth(8);
	oFile.SetSize(m_iSize,m_iSize);
	CreateGrayscaleColorTable(oFile);
	for(int i = 0; i < m_iSize; i++)
	{
		for(int j = 0; j < m_iSize; j++)
		{
			unsigned char o = (unsigned char)( data[ i + j * m_iSize ] * MAX_HEIGHT );
			oFile(j,i)->Red = o;
			oFile(j,i)->Green = o;
			oFile(j,i)->Blue = o;
		}
	}

	oFile.WriteToFile(file.c_str());

} // writeToFile

// --------------------------------------------------------------------------------------------------------------------------------
void TerrainGenerator::makeTriangleDivision(float roughness, int seed, float firValue)
{
	int i, j;

	int	rectSize = m_iSize;
	float	dh = (float)m_iSize * 0.5f;
	
	srand( seed );

	// random values for the edge points
	m_fTerrain[computeIndex(0, 0)] = random( -dh, +dh );
	m_fTerrain[computeIndex(m_iSize-1, 0)] = random( -dh, +dh );
	m_fTerrain[computeIndex(0, m_iSize-1)] = random( -dh, +dh );
	m_fTerrain[computeIndex(m_iSize-1, m_iSize-1)] = random( -dh, +dh );

	while( rectSize > 0 )
	{
		dh  = pow(static_cast<float>(rectSize), roughness);	//k^r

		for ( i = 0; i < m_iSize; i += rectSize )
			for ( j = 0; j < m_iSize; j += rectSize )
			{
				int x = ( i + rectSize ) & m_iSizeMask;
				int y = ( j + rectSize ) & m_iSizeMask;

				int mx = i + rectSize / 2;
				int my = j + rectSize / 2;

				m_fTerrain[ computeIndex( mx, my ) ] =
					0.5f * (m_fTerrain[ computeIndex( x, j ) ] + m_fTerrain[ computeIndex( i, y ) ]) + random( -dh, +dh );

				// oben
				m_fTerrain[ computeIndex( mx, j ) ] =
					0.5f * ( m_fTerrain[ computeIndex( i, j ) ] + m_fTerrain[ computeIndex( x, j ) ]) + random( -dh, +dh );

				// links
				m_fTerrain[ computeIndex( i, my ) ] =
					0.5f * ( m_fTerrain[ computeIndex( i, j ) ] + m_fTerrain[ computeIndex( i, y ) ]) + random( -dh, +dh );
			} // for

			rectSize >>= 1;
	} // while

	if(firValue>0.0f)
		filterTerrain(firValue);
	normalize();
} // makeTriangleDivision

// --------------------------------------------------------------------------------------------------------------------------------
void TerrainGenerator::makeDiamondSquare(float roughness, int seed, float firValue)
{
	int i, j;

	int	rectSize = m_iSize;
	float	dh = (float)m_iSize * 0.25f;

	srand( seed );

	float it = 1.0f;

	while( rectSize > 0 )
	{
		dh  = pow(roughness, it);			//r^i

		// Diamond
		for ( i = 0; i < m_iSize; i += rectSize )
			for ( j = 0; j < m_iSize; j += rectSize )
			{
				int x = ( i + rectSize ) & m_iSizeMask;
				int y = ( j + rectSize ) & m_iSizeMask;

				int mx = i + rectSize / 2;
				int my = j + rectSize / 2;

				m_fTerrain[ computeIndex( mx, my ) ] =
					0.25f * ( m_fTerrain[ computeIndex( i, j ) ] + m_fTerrain[ computeIndex( x, j ) ] +
					m_fTerrain[ computeIndex( i, y ) ] + m_fTerrain[ computeIndex( x, y ) ] ) +
					random( -dh, +dh );
			} // for

		// Square
		for ( i = 0; i < m_iSize; i += rectSize )
			for ( j = 0; j < m_iSize; j += rectSize )
			{
				int x = ( i + rectSize ) & m_iSizeMask;
				int y = ( j + rectSize ) & m_iSizeMask;

				int mx = i + rectSize / 2;
				int my = j + rectSize / 2;

				int sx = ( i - rectSize / 2 + m_iSize ) & m_iSizeMask;
				int sy = ( j - rectSize / 2 + m_iSize ) & m_iSizeMask;

				// oben
				m_fTerrain[ computeIndex( mx, j ) ] =
					0.25f * ( m_fTerrain[ computeIndex( i, j ) ] + m_fTerrain[ computeIndex( x, j ) ] +
					m_fTerrain[ computeIndex( mx, sy ) ] + m_fTerrain[ computeIndex( mx, my ) ] ) +
					random( -dh, +dh );

				// links
				m_fTerrain[ computeIndex( i, my ) ] =
					0.25f * ( m_fTerrain[ computeIndex( i, j ) ] + m_fTerrain[ computeIndex( i, y ) ] +
					m_fTerrain[ computeIndex( sx, my ) ] + m_fTerrain[ computeIndex( mx, my ) ] ) +
					random( -dh, +dh );
			} // for

		rectSize >>= 1;
		it++;
	} // while

	if(firValue>0.0f)
		filterTerrain(firValue);
	normalize();
} // makeDiamondSquare

// --------------------------------------------------------------------------------------------------------------------------------
void TerrainGenerator::makeFourierTransformation(float roughness)
{
	double *fRe, *fIm; //the signal's real part, imaginary part, and amplitude
	double *FRe, *FIm; //the FT's real part, imaginary part and amplitude
	double *fRe2, *fIm2; //will become the signal again after IDFT of the spectrum
	double *FRe2, *FIm2; //filtered spectrum

	fRe = new double[m_iSize * m_iSize * 3];
	memset( fRe, 0, m_iSize2 * 3 * sizeof( double ) );	// fill the array with 0
	fIm = new double[m_iSize * m_iSize * 3];
	memset( fIm, 0, m_iSize2 * 3 * sizeof( double ) );	// fill the array with 0

	FRe = new double[m_iSize * m_iSize * 3];
	memset( FRe, 0, m_iSize2 * 3 * sizeof( double ) );	// fill the array with 0
	FIm = new double[m_iSize * m_iSize * 3];
	memset( FIm, 0, m_iSize2 * 3 * sizeof( double ) );	// fill the array with 0

	fRe2 = new double[m_iSize * m_iSize * 3];
	memset( fRe2, 0, m_iSize2 * 3 * sizeof( double ) );	// fill the array with 0
	fIm2 = new double[m_iSize * m_iSize * 3];
	memset( fIm2, 0, m_iSize2 * 3 * sizeof( double ) );	// fill the array with 0

	FRe2 = new double[m_iSize * m_iSize * 3];
	memset( FRe2, 0, m_iSize2 * 3 * sizeof( double ) );	// fill the array with 0
	FIm2 = new double[m_iSize * m_iSize * 3];
	memset( FIm2, 0, m_iSize2 * 3 * sizeof( double ) );	// fill the array with 0

	// create random points
	for(int x = 0; x < m_iSize2; ++x) 
	{
		double v = static_cast<double>(random(0.0f, 1.0f));
		fRe[x + 0 * 3] = v;
		fRe[x + 1 * 3] = v;
		fRe[x + 2 * 3] = v;
	}

	// 2D FFT
	FFT2D(m_iSize, m_iSize, false, fRe, fIm, FRe, FIm);

	// scaling
	for(int x = 0; x < m_iSize; ++x)
		for(int y = 0; y < m_iSize; ++y) 
		{
			double factor = std::max(1.0, 128.0 - (abs(x-(m_iSize/2)) + abs(y-(m_iSize/2)))) / (m_iSize/2);
			double v = FRe[x + y * m_iSize];

			double newv = v * (1 / (pow(factor,static_cast<double>(roughness))));;					// value of frequency domain * 1/f noise
			
			FRe[x + y * m_iSize + 0 * 3] = newv;
			FRe[x + y * m_iSize + 1 * 3] = newv;
			FRe[x + y * m_iSize + 2 * 3] = newv;

			double im_v = FIm[x + y * m_iSize];

			newv =  im_v * (1 / (pow(factor,static_cast<double>(roughness))));

			FIm[x + y * m_iSize + 0 * 3] = newv;
			FIm[x + y * m_iSize + 1 * 3] = newv;
			FIm[x + y * m_iSize + 2 * 3] = newv;
		}

	// inverse 2D FFT
	FFT2D(m_iSize, m_iSize, true, FRe, FIm, fRe2, fIm2);

	// get the Terrain
	for(int x = 0; x < m_iSize; ++x) 
		for(int y = 0; y < m_iSize; ++y) 
		{
			float v = static_cast<float>(fRe2[x + y * m_iSize]);
			m_fTerrain[computeIndex(x,y)] = v;
			//std::cout << v << std::endl;
		}

	normalize();

	delete[] fRe;
	delete[] fIm;
	delete[] FRe;
	delete[] FIm;
	delete[] fRe2;
	delete[] fIm2;
	delete[] FRe2;
	delete[] FIm2;
} // makeFourierTransformation

void TerrainGenerator::FFT2D(int n, int m, bool inverse, double *gRe, double *gIm, double *GRe, double *GIm)
{
	int l2n = 0, p = 1; //l2n will become log_2(n)
	while(p < n) {p *= 2; l2n++;}
	int l2m = 0; p = 1; //l2m will become log_2(m)
	while(p < m) {p *= 2; l2m++;}

	m = 1 << l2m; n = 1 << l2n; //Make sure m and n will be powers of 2, otherwise you'llget in an infinite loop

	//Erase all history of this array
	for(int x = 0; x < m; x++) //for each column
		for(int y = 0; y < m; y++) //for each row
			for(int c = 0; c < 3; c++) //for each color component
			{
				GRe[3 * m * x + 3 * y + c] = gRe[3 * m * x + 3 * y + c];
				GIm[3 * m * x + 3 * y + c] = gIm[3 * m * x + 3 * y + c];
			}

	//Bit reversal of each row
	int j;
	for(int y = 0; y < m; y++) //for each row
		for(int c = 0; c < 3; c++) //for each color component
		{
			j = 0;
			for(int i = 0; i < n - 1; i++)
			{
				GRe[3 * m * i + 3 * y + c] = gRe[3 * m * j + 3 * y + c];
				GIm[3 * m * i + 3 * y + c] = gIm[3 * m * j + 3 * y + c];
				int k = n / 2;
				while (k <= j) {j -= k; k/= 2;}
				j += k;
			}
		}

	//Bit reversal of each column
	double tx = 0, ty = 0;  
	for(int x = 0; x < n; x++) //for each column
		for(int c = 0; c < 3; c++) //for each color component
		{
			j = 0;
			for(int i = 0; i < m - 1; i++)
			{
				if(i < j)
				{
					tx = GRe[3 * m * x + 3 * i + c];
					ty = GIm[3 * m * x + 3 * i + c];
					GRe[3 * m * x + 3 * i + c] = GRe[3 * m * x + 3 * j + c];
					GIm[3 * m * x + 3 * i + c] = GIm[3 * m * x + 3 * j + c];
					GRe[3 * m * x + 3 * j + c] = tx;
					GIm[3 * m * x + 3 * j + c] = ty;
				}
				int k = m / 2;
				while (k <= j) {j -= k; k/= 2;}
				j += k;
			}
		}

	//Calculate the FFT of the columns
	for(int x = 0; x < n; x++) //for each column
		for(int c = 0; c < 3; c++) //for each color component
		{
			//This is the 1D FFT:
			double ca = -1.0;
			double sa = 0.0;
			int l1 = 1, l2 = 1;
			for(int l=0;l < l2n;l++)
			{
				l1 = l2;
				l2 *= 2;
				double u1 = 1.0;
				double u2 = 0.0;
				for(int j = 0; j < l1; j++)
				{
					for(int i = j; i < n; i += l2)
					{
						int i1 = i + l1;
						double t1 = u1 * GRe[3 * m * x + 3 * i1 + c] - u2 * GIm[3 * m * x + 3 * i1 + c];
						double t2 = u1 * GIm[3 * m * x + 3 * i1 + c] + u2 * GRe[3 * m * x + 3 * i1 + c];
						GRe[3 * m * x + 3 * i1 + c] = GRe[3 * m * x + 3 * i + c] - t1;
						GIm[3 * m * x + 3 * i1 + c] = GIm[3 * m * x + 3 * i + c] - t2;
						GRe[3 * m * x + 3 * i + c] += t1;
						GIm[3 * m * x + 3 * i + c] += t2;
					}
					double z =  u1 * ca - u2 * sa;
					u2 = u1 * sa + u2 * ca;
					u1 = z;
				}
				sa = sqrt((1.0 - ca) / 2.0);
				if(!inverse) sa = -sa;
				ca = sqrt((1.0 + ca) / 2.0);
			}
		}
	//Calculate the FFT of the rows
	for(int y = 0; y < m; y++) //for each row
		for(int c = 0; c < 3; c++) //for each color component
		{
			//This is the 1D FFT:
			double ca = -1.0;
			double sa = 0.0;
			int l1= 1, l2 = 1;
			for(int l = 0; l < l2m; l++)
			{
				l1 = l2;
				l2 *= 2;
				double u1 = 1.0;
				double u2 = 0.0;
				for(int j = 0; j < l1; j++)
				{
					for(int i = j; i < n; i += l2)
					{
						int i1 = i + l1;
						double t1 = u1 * GRe[3 * m * i1 + 3 * y + c] - u2 * GIm[3 * m * i1 + 3 * y + c];
						double t2 = u1 * GIm[3 * m * i1 + 3 * y + c] + u2 * GRe[3 * m * i1 + 3 * y + c];
						GRe[3 * m * i1 + 3 * y + c] = GRe[3 * m * i + 3 * y + c] - t1;
						GIm[3 * m * i1 + 3 * y + c] = GIm[3 * m * i + 3 * y + c] - t2;
						GRe[3 * m * i + 3 * y + c] += t1;
						GIm[3 * m * i + 3 * y + c] += t2;
					}
					double z =  u1 * ca - u2 * sa;
					u2 = u1 * sa + u2 * ca;
					u1 = z;
				}
				sa = sqrt((1.0 - ca) / 2.0);
				if(!inverse) sa = -sa;
				ca = sqrt((1.0 + ca) / 2.0);
			}
		}

		int d;
		if(inverse) d = n; else d = m;
		for(int x = 0; x < n; x++) for(int y = 0; y < m; y++) for(int c = 0; c < 3; c++) //for every value of the buffers
		{
			GRe[3 * m * x + 3 * y + c] /= d;
			GIm[3 * m * x + 3 * y + c] /= d;
		}
}

// --------------------------------------------------------------------------------------------------------------------------------
void TerrainGenerator::makeFaultFormation(int iterations, int filterIterations, float firValue)
{
	for ( int i = 0; i < iterations; i ++ )
	{
		// Höhendifferenz für diese Interation berechnen
		float	heightDifference = 100.0f * ( 1.0f - (float)i / (float)iterations );

		// zwei zufällige, ungleiche punkte auf der heightmap suchen
		int	x1, x2, y1, y2;

		x1 = rand()&m_iSizeMask;
		y1 = rand()&m_iSizeMask;

		do
		{
			x2 = rand()&m_iSizeMask;
			y2 = rand()&m_iSizeMask;
		} while ( x2==x1 && y2==y1 );

		float	dx = (float)(x2 - x1);
		float	dy = (float)(y2 - y1);

		int	upDown = (dx>0&&dy<0)||(dx>0&&dy>0);

		if ( dx )
			dy /= dx; else dy = 0.0f;

		float	x, y;

		x = 0;
		y = y1 - x1 * dy;

		// und alle Punkte einer Seite (upDown entscheidet !) erhöhen
		for ( x2 = 0; x2 < m_iSize; x2 ++, y += dy )
		{
			for ( y2 = 0; y2 < m_iSize; y2 ++ )
				if ( (  upDown && y2 < y ) /*up*/ ||
					( !upDown && y2 > y ) /*down*/ )
					m_fTerrain[ x2 + y2 * m_iSize ] += heightDifference;
		}

		// und wenn gewünscht filtern
		if(firValue>0.0f)
			if ( (i % filterIterations) == 0 && filterIterations )
				filterTerrain(firValue);
	}

	normalize();
} // makeFaultFormation

// --------------------------------------------------------------------------------------------------------------------------------
void TerrainGenerator::makeMidpointDisplacement(float roughness, int seed, float firValue)
{
	int i, j;

	int	rectSize = m_iSize;
	float	dh = (float)m_iSize * 0.25f;

	srand( seed );

	float it = 1.0f;

	while( rectSize > 0 )
	{
		dh  = pow(roughness, it);			//r^i

		// Square
		for ( i = 0; i < m_iSize; i += rectSize )
			for ( j = 0; j < m_iSize; j += rectSize )
			{
				int x = ( i + rectSize ) & m_iSizeMask;
				int y = ( j + rectSize ) & m_iSizeMask;

				int mx = i + rectSize / 2;
				int my = j + rectSize / 2;

				int sx = ( i - rectSize / 2 + m_iSize ) & m_iSizeMask;
				int sy = ( j - rectSize / 2 + m_iSize ) & m_iSizeMask;

				// mitte
				m_fTerrain[ computeIndex( mx, my ) ] =
					0.25f * ( m_fTerrain[ computeIndex( i, j ) ] + m_fTerrain[ computeIndex( x, j ) ] +
					m_fTerrain[ computeIndex( i, y ) ] + m_fTerrain[ computeIndex( x, y ) ] ) +
					random( -dh, +dh );

				// oben
				m_fTerrain[ computeIndex( mx, j ) ] =
					0.5f * ( m_fTerrain[ computeIndex( i, j ) ] + m_fTerrain[ computeIndex( x, j ) ]) +
					random( -dh, +dh );

				// links
				m_fTerrain[ computeIndex( i, my ) ] =
					0.5f * ( m_fTerrain[ computeIndex( i, j ) ] + m_fTerrain[ computeIndex( i, y ) ]) +
					random( -dh, +dh );
			} // for

			rectSize >>= 1;
			it++;
	} // while

	if(firValue>0.0f)
		filterTerrain(firValue);
	normalize();
} // makeMidpointDisplacement

void TerrainGenerator::makeParticleDeposition(int nMountain, int moveDrop, int particle, float caldera, float firValue)
{
	const int DX[] = { 0, 1, 0, m_iSize-1, 1, 1, m_iSize-1, m_iSize-1 };
	const int DY[] = { 1, 0, m_iSize-1, 0, m_iSize-1, 1, m_iSize-1, 1 };

	for ( int m = 0; m < nMountain; m++ )
	{
		unsigned int x = rand() & m_iSizeMask;
		unsigned int y = rand() & m_iSizeMask;

		int topX = x, topY = y;

		int nParticles = particle;

		for ( int i = 0; i < nParticles; i++ )
		{
			// eventuell den Particle Startpunkt verschieben
			if ( moveDrop && ( i % moveDrop ) == 0 )
			{
				int dir = rand() & 7;
				x = ( x + DX[ dir ] ) & m_iSizeMask;
				y = ( y + DY[ dir ] ) & m_iSizeMask;
			}

			m_fTerrain[ computeIndex( x, y ) ] ++;

			int	px = x;
			int	py = y;
			int ok = 0;

			while ( !ok++ )
			{
				int dir = rand();

				for ( int j = 0; j < 8; j++ )
				{
					int	ofs = ( j + m ) & 7;
					int	tx = ( px + DX[ ofs ] ) & m_iSizeMask;
					int	ty = ( py + DY[ ofs ] ) & m_iSizeMask;

					// wenn der Nachbarpunkt niedriger liegt, Partikel bewegen
					if ( m_fTerrain[ computeIndex( px, py ) ] > m_fTerrain[ computeIndex( tx, ty ) ] + 1.0f )
					{
						m_fTerrain[ computeIndex( tx, ty ) ] ++;
						m_fTerrain[ computeIndex( px, py ) ] --;
						// neue Partikelposition
						px = tx;
						py = ty;
						ok = 0;
						break;
					}
				}
			}

			// den Gipfelpunkt speichern
			if ( m_fTerrain[ computeIndex( px, py ) ] > m_fTerrain[ computeIndex( topX, topY ) ] )
			{
				topX = px;
				topY = py;
			}
		}

		// nun alle Punkte um den Gipfel, die höher als calderaLine liegen
		// umklappen !
		float	calderaLine = m_fTerrain[ computeIndex( topX, topY ) ] - caldera;

		if ( calderaLine > 0 )
			createCalderas(topX, topY, calderaLine);

	}

	if(firValue>0.0f)
		filterTerrain(firValue);
	normalize();
} // makeParticleDeposition

// --------------------------------------------------------------------------------------------------------------------------------
void TerrainGenerator::makePerlinNoise(float persistence, float frequency, float amplitude, int octaves, int seed, float firValue)
{
	PerlinNoise p(static_cast<double>(persistence), static_cast<double>(frequency), static_cast<double>(amplitude), octaves, seed);
	
	for(int i=0; i<m_iSize; ++i)
		for(int j=0; j<m_iSize; ++j)
			m_fTerrain[computeIndex(i,j)] = static_cast<float>(p.GetHeight(static_cast<double>(i),static_cast<double>(j)));

	if(firValue>0.0f)
		filterTerrain(firValue);

	normalize();
} // makePerlinNoise

void TerrainGenerator::makeVoronoiDiagram(int points, int seed, float firValue)
{
	srand(seed);

	std::vector<vpoint> voronoi_points;

	// set random points
	for(int i=0; i<points; ++i)
	{
		vpoint p;
		p.x = irand(0,m_iSize);
		p.y = irand(0,m_iSize);
		voronoi_points.push_back(p);
	}

	float c1 = 1.0f;
	float c2 = -1.0f;

	for(int i=0; i<m_iSize; ++i)
		for(int j=0; j<m_iSize; ++j)
		{
			if(m_fTerrain[computeIndex(i,j)] == 0.0f)
			{
				// search next point and set color
				float fmax = 999999.0f;
				float closest_Point = fmax;
				float second_closest_Point = fmax;

				for(std::vector<vpoint>::iterator it = voronoi_points.begin(); it!=voronoi_points.end(); ++it)
				{
					float distance = static_cast<float>(abs(i - (*it).x) + abs(j - (*it).y));
					if(distance<closest_Point)
					{
						second_closest_Point = closest_Point;
						closest_Point = distance;		
					}
					else
					{
						if(distance<second_closest_Point)
							second_closest_Point = distance;
					}
				}
				
				float d1 = MAX_HEIGHT - closest_Point;
				float d2 = MAX_HEIGHT - second_closest_Point;

				m_fTerrain[computeIndex(i,j)] = c1*d1 + c2*d2;
			}
		}

	if(firValue>0.0f)
		filterTerrain(firValue);
	normalize();
}

// --------------------------------------------------------------------------------------------------------------------------------
// ein einfacher FIR Lowpass-Filter
void TerrainGenerator::filterFIR(float *value, int ofs, float filter)
{
	float v  =  value[ 0 ];
	float *b = &value[ ofs ];
	float ifilter = 1.0f - filter;

	for ( int i = 0; i < m_iSizeMask; i++, b += ofs )
	{
		*b = filter * v + ifilter * *b;
		v = *b;
	}
}

// --------------------------------------------------------------------------------------------------------------------------------
// die Landschaft in alle 4 Richtungen filtern
void TerrainGenerator::filterTerrain(float filter)
{
	// die Zeilen von links nach rechts und umgekehrt filtern
	for ( int i = 0; i < m_iSize; i++ )
	{
		filterFIR( &m_fTerrain[ m_iSize * i ], 1, filter);
		filterFIR( &m_fTerrain[ m_iSize * i + m_iSize - 1 ], -1, filter );
	}

	// die Spalten von oben nach unten und umgekehrt filtern
	for ( int i = 0; i < m_iSize; i++ )
	{
		filterFIR( &m_fTerrain[i], m_iSize, filter );
		filterFIR( &m_fTerrain[ m_iSize * ( m_iSize - 1 ) + i ], -m_iSize, filter );
	}
}

// --------------------------------------------------------------------------------------------------------------------------------
void TerrainGenerator::createCalderas(int x, int y, float height)
{
	if ( x < 0 || x >= m_iSize || y < 0 || y >= m_iSize )
		return;

	if ( m_fTerrain[ computeIndex( x, y ) ] > height )
	{
		m_fTerrain[ computeIndex( x, y ) ] -= ( m_fTerrain[ computeIndex( x, y ) ] - height ) * 2.0f;

		createCalderas(x + 1, y, height);
		createCalderas(x - 1, y, height);
		createCalderas(x, y + 1, height);
		createCalderas(x, y - 1, height);
	}
}

// --------------------------------------------------------------------------------------------------------------------------------
void	TerrainGenerator::normalize()
{
	float	minHeight;
	float	maxHeight;
	float	*t = m_fTerrain;
	int i;

	minHeight = maxHeight = *t;

	// Minima und Maxima suchen
	for ( i = 1; i < m_iSize2; i++, t++ )
	{
		if ( *t < minHeight )
			minHeight = *t; else
			if ( *t > maxHeight )
				maxHeight = *t;
	}

	// und Landschaft entsprechend skalieren
	float	scale = 1.0f / ( maxHeight - minHeight );

	t = m_fTerrain;
	for ( i = 0; i < m_iSize2; i++, t++ )
		*t = ( *t - minHeight ) * scale;
}

// --------------------------------------------------------------------------------------------------------------------------------
void TerrainGenerator::makeThermalErosion(float talus, int iterations)
{
	for(int k = 0; k < iterations; ++k)
	{
		// step through the heightmap
		for(int i = 0; i < m_iSize; ++i)
			for(int j = 0; j < m_iSize; ++j)
			{
				float height = m_fTerrain[computeIndex(i,j)];
				// there is material to move
				if(height > talus)
				{
					// find out how much material should move
					float material = height - talus;
					material = material * random(0.01f,0.09f);			// move between 1% and 9%

					// bring that material to the neighbour cells
					// use Von Neumann neighbourhood
					
					// count lower neighbours:
					int lower_neighbours = 0;
					bool up = false,down = false,left = false,right = false;
					// up
					if(i > 0)
						if(m_fTerrain[computeIndex(i-1,j)] < height)
						{
							up = true;
							lower_neighbours++;
						}
					// down
					if(i < m_iSize-1)
						if(m_fTerrain[computeIndex(i+1,j)] < height)
						{
							down = true;
							lower_neighbours++;
						}
					// left
					if(j > 0)
						if(m_fTerrain[computeIndex(i,j-1)] < height)
						{
							left = true;
							lower_neighbours++;
						}
					// right
					if(j < m_iSize-1)
						if(m_fTerrain[computeIndex(i,j+1)] < height)
						{
							right = true;
							lower_neighbours++;
						}

					// only if there are lower neighbours
					if(lower_neighbours > 0)
					{
						float material_part = static_cast<float>(material / lower_neighbours);

						// up
						if(up)
							m_fTerrain[computeIndex(i-1,j)] += material_part;
						// down
						if(down)
							m_fTerrain[computeIndex(i+1,j)] += material_part;
						// left
						if(left)
							m_fTerrain[computeIndex(i,j-1)] += material_part;
						// right
						if(right)
							m_fTerrain[computeIndex(i,j+1)] += material_part;

						m_fTerrain[computeIndex(i,j)] -= material;
					}
				}
			}
	}
	
} // makeThermalErosion

// --------------------------------------------------------------------------------------------------------------------------------
void TerrainGenerator::makeHydraulicErosion(float water, float sediment, float evaporation, float capacity, int iterations)
{
	// create maps for water and sediments
	float *watermap = new float[ m_iSize2 ];
	memset( watermap, 0, m_iSize2 * sizeof( float ) );

	float *sedimentmap = new float[ m_iSize2 ];
	memset( sedimentmap, 0, m_iSize2 * sizeof( float ) );

	for(int k = 0; k < iterations; ++k)
	{
		// 1.
		// generate new water
		for(int i = 0; i < m_iSize2; ++i)
			watermap[i] += water;

		// 2.
		// convert height material to sediment proportional to water
		for(int i = 0; i < m_iSize2; ++i)
		{
			float amount = watermap[i] * sediment;
			// if more sediment should be generated as height, only use height
			if(amount > m_fTerrain[i])
				amount = m_fTerrain[i];
			m_fTerrain[i] -= amount;
			sedimentmap[i] += amount;
		}

		// 3.
		// water and sediment distributing to neighbour cells
		for(int i = 0; i < m_iSize; ++i)
			for(int j = 0; j < m_iSize; ++j)
			{
				float height = m_fTerrain[computeIndex(i,j)];
				float water = watermap[computeIndex(i,j)];

				float a = height + water;
				float a_average = a;

				// count neighbours
				// use Von Neumann neighbourhood
				int lower_neighbours = 0;
				bool up = false,down = false,left = false,right = false;
				float d1 = 0.0f, d2 = 0.0f, d3 = 0.0f, d4 = 0.0f;

				// up
				if(i > 0)
				{
					d1 = a - (m_fTerrain[computeIndex(i-1,j)] + watermap[computeIndex(i-1,j)]);
					if(d1 > 0.0f)
					{
						a_average += m_fTerrain[computeIndex(i-1,j)] + watermap[computeIndex(i-1,j)];
						up = true;
						lower_neighbours++;
					}
				}
				// down
				if(i < m_iSize-1)
				{
					d2 = a - (m_fTerrain[computeIndex(i+1,j)] + watermap[computeIndex(i+1,j)]);
					if(d2 > 0.0f)
					{
						a_average += m_fTerrain[computeIndex(i+1,j)] + watermap[computeIndex(i+1,j)];
						down = true;
						lower_neighbours++;
					}
				}
				// left
				if(j > 0)
				{
					d3 = a - (m_fTerrain[computeIndex(i,j-1)] + watermap[computeIndex(i,j-1)]);
					if(d3 > 0.0f)
					{
						a_average += m_fTerrain[computeIndex(i,j-1)] + watermap[computeIndex(i,j-1)];
						left = true;
						lower_neighbours++;
					}
				}
				// right
				if(j < m_iSize-1)
				{
					d4 = a - (m_fTerrain[computeIndex(i,j+1)] + watermap[computeIndex(i,j+1)]);	
					if(d4 > 0.0f)
					{
						a_average += m_fTerrain[computeIndex(i,j+1)] + watermap[computeIndex(i,j+1)];
						right = true;
						lower_neighbours++;
					}
				}

				// move water only if there are lower_neighbours
				if(lower_neighbours > 0)
				{
					a_average = a_average / (lower_neighbours + 1);

					// part of the water moves
					float d_total = 0.0f;
					float water_sum = 0.0f;
					float water_part = 0.0f;
					float sediment_to_move = 0.0f;
					float sediment_part = 0.0f;

					float a_delta = a - a_average;

					if(d1 > 0.0f) d_total += d1;
					if(d2 > 0.0f) d_total += d2;
					if(d3 > 0.0f) d_total += d3;
					if(d4 > 0.0f) d_total += d4;

					// up
					if(up)
					{
						water_part = std::min(water, a_delta) * (d1 / d_total);
						water_sum += water_part;
						sediment_part = sedimentmap[computeIndex(i,j)] * (water_part / water);
						watermap[computeIndex(i-1,j)] += water_part;
						sedimentmap[computeIndex(i-1,j)] += sediment_part;
						sediment_to_move += sediment_part;
					}
					// down
					if(down)
					{
						water_part = std::min(water, a_delta) * (d2 / d_total);
						water_sum += water_part;
						sediment_part = sedimentmap[computeIndex(i,j)] * (water_part / water);
						watermap[computeIndex(i+1,j)] += water_part;
						sedimentmap[computeIndex(i+1,j)] += sediment_part;
						sediment_to_move += sediment_part;
					}
					// left
					if(left)
					{
						water_part = std::min(water, a_delta) * (d3 / d_total);
						water_sum += water_part;
						sediment_part = sedimentmap[computeIndex(i,j)] * (water_part / water);
						watermap[computeIndex(i,j-1)] += water_part;
						sedimentmap[computeIndex(i,j-1)] += sediment_part;
						sediment_to_move += sediment_part;
					}
					// right
					if(right)
					{
						water_part = std::min(water, a_delta) * (d4 / d_total);
						water_sum += water_part;
						sediment_part = sedimentmap[computeIndex(i,j)] * (water_part / water);
						watermap[computeIndex(i,j+1)] += water_part;
						sedimentmap[computeIndex(i,j+1)] += sediment_part;
						sediment_to_move += sediment_part;
					}

					if(water_sum > watermap[computeIndex(i,j)])
						water_sum = watermap[computeIndex(i,j)];
					watermap[computeIndex(i,j)] -= water_sum;
					sedimentmap[computeIndex(i,j)] -= sediment_to_move;
					if(sedimentmap[computeIndex(i,j)] < 0.0f)
						sedimentmap[computeIndex(i,j)] = 0.0f;
				}
			}

			// 4.
			// amount of water evaporates and sediment drops out
			for(int i = 0; i < m_iSize2; ++i)
			{
				if(watermap[i] < 0.0f)
					watermap[i] = 0.0f;

				// evaporation of water
				watermap[i] = watermap[i] * (1.0f - evaporation);

				// maximum amount sediment can be carried by water
				float sediment_max = capacity * watermap[i];

				// if it is too much sediment --> move back to heightmap
				float sediment_overhead = std::max(0.0f, sedimentmap[i] - sediment_max);
				if(sediment_overhead > sedimentmap[i])
					sediment_overhead = sedimentmap[i];
				sedimentmap[i] -= sediment_overhead;
				m_fTerrain[i] += sediment_overhead;
			}
		}

		// now write sediments back to heightmap
		for(int i = 0; i < m_iSize2; ++i)
			m_fTerrain[i] += sedimentmap[i];

} // makeHydraulikErosion

void TerrainGenerator::createGameMap(float threshold, float *data)
{
	// use Von Neumann neighbours
	float slopevalue;
	float height;
	float help;

	for(int i = 0; i < m_iSize; ++i)
		for(int j = 0; j < m_iSize; ++j)
		{
			height = m_fTerrain[computeIndex(i,j)];
			slopevalue = 0.0f;

			// up
			if(i > 0)
			{
				help = std::abs(height - m_fTerrain[computeIndex(i-1,j)]);
				if(slopevalue < help)
					slopevalue = help;
			}
			// down
			if(i < m_iSize-1)
			{
				help = std::abs(height - m_fTerrain[computeIndex(i+1,j)]);
				if(slopevalue < help)
					slopevalue = help;
			}
			// left
			if(j > 0)
			{
				help = std::abs(height - m_fTerrain[computeIndex(i,j-1)]);
				if(slopevalue < help)
					slopevalue = help;
			}
			// right
			if(j < m_iSize-1)
			{
				help = std::abs(height - m_fTerrain[computeIndex(i,j+1)]);
				if(slopevalue < help)
					slopevalue = help;
			}

			// is cell accessible?
			if(slopevalue < threshold)
				data[computeIndex(i,j)] = 1.0f;
		}
}

void TerrainGenerator::makeAccessabilityMap(std::string filename, float threshold)
{
	float *accessibilitymap = new float[ m_iSize2 ];
	memset( accessibilitymap, 0, m_iSize2 * sizeof( float ) );

	clock_t start, end;
	start = clock();
	createGameMap(threshold, accessibilitymap);
	end = clock();
	std::cout << "Runtime (AccessabilityMap): " << start << " to " << end << " --> " << end-start << " / " << CLOCKS_PER_SEC << " seconds" << std::endl;

	std::string file = filename.append("_accessabilityMap");
	writeToFile(file, accessibilitymap);
} // makeAccessabilityMap

void TerrainGenerator::makeFlatnessMap(std::string filename, float threshold)
{
	float *flatnessmap = new float[ m_iSize2 ];
	memset( flatnessmap, 0, m_iSize2 * sizeof( float ) );

	clock_t start, end;
	start = clock();
	createGameMap(threshold, flatnessmap);
	end = clock();
	std::cout << "Runtime (FlatnessMap): " << start << " to " << end << " --> " << end-start << " / " << CLOCKS_PER_SEC << " seconds" << std::endl;

	std::string file = filename.append("_flatnessMap");
	writeToFile(file, flatnessmap);
} // makeFlatnessMap

void TerrainGenerator::makeProceduralTexture(std::string filename)
{
	std::string file = filename.append(".bmp");

	// load Textures
	BMP img_sand;
	img_sand.ReadFromFile("Data/sand.bmp");
	BMP img_gras;
	img_gras.ReadFromFile("Data/gras.bmp");
	BMP img_rock;
	img_rock.ReadFromFile("Data/rock.bmp");
	BMP img_snow;
	img_snow.ReadFromFile("Data/snow.bmp");

	// generate new texture
	BMP oFile;
	oFile.SetBitDepth(img_sand.TellBitDepth());
	oFile.SetSize(m_iSize,m_iSize);
	//CreateGrayscaleColorTable(oFile);

	clock_t start, end;
	start = clock();

	int x,y;
	for(int i = 0; i < m_iSize; i++)
	{
		for(int j = 0; j < m_iSize; j++)
		{
			float height = m_fTerrain[ i + j * m_iSize ] * MAX_HEIGHT;
			x=j;
			y=i;
			if(height<=53.0f)
			{
				wrapAround(x,y,img_sand.TellWidth(),img_sand.TellHeight());

				oFile(j,i)->Red = img_sand(x,y)->Red;
				oFile(j,i)->Green = img_sand(x,y)->Green;
				oFile(j,i)->Blue = img_sand(x,y)->Blue;
			}
			// blend Textures
			if(height>53.0f && height<=73.0f)
			{
				float factorA = (height - 53.0f)/20.0f;
				float factorB = 1.0f - factorA;

				// Texture 1
				wrapAround(x,y,img_sand.TellWidth(),img_sand.TellHeight());

				oFile(j,i)->Red = img_sand(x,y)->Red * factorB;
				oFile(j,i)->Green = img_sand(x,y)->Green * factorB;
				oFile(j,i)->Blue = img_sand(x,y)->Blue * factorB;

				x=j;
				y=i;

				// Texture 2
				wrapAround(x,y,img_gras.TellWidth(),img_gras.TellHeight());

				oFile(j,i)->Red += img_gras(x,y)->Red * factorA;
				oFile(j,i)->Green += img_gras(x,y)->Green * factorA;
				oFile(j,i)->Blue += img_gras(x,y)->Blue * factorA;
			}
			if(height>73.0f && height<=117.0f)
			{
				wrapAround(x,y,img_gras.TellWidth(),img_gras.TellHeight());

				oFile(j,i)->Red = img_gras(x,y)->Red;
				oFile(j,i)->Green = img_gras(x,y)->Green;
				oFile(j,i)->Blue = img_gras(x,y)->Blue;
			}
			// blend Textures
			if(height>117.0f && height<=137.0f)
			{
				float factorA = (height - 117.0f)/20.0f;
				float factorB = 1.0f - factorA;

				//Texture 1
				wrapAround(x,y,img_gras.TellWidth(),img_gras.TellHeight());

				oFile(j,i)->Red = img_gras(x,y)->Red * factorB;
				oFile(j,i)->Green = img_gras(x,y)->Green * factorB;
				oFile(j,i)->Blue = img_gras(x,y)->Blue * factorB;

				x=j;
				y=i;

				//Texture 2
				wrapAround(x,y,img_rock.TellWidth(),img_rock.TellHeight());

				oFile(j,i)->Red += img_rock(x,y)->Red * factorA;
				oFile(j,i)->Green += img_rock(x,y)->Green * factorA;
				oFile(j,i)->Blue += img_rock(x,y)->Blue * factorA;
			}
			if(height>127.0f && height<=181.0f)
			{
				wrapAround(x,y,img_rock.TellWidth(),img_rock.TellHeight());

				oFile(j,i)->Red = img_rock(x,y)->Red;
				oFile(j,i)->Green = img_rock(x,y)->Green;
				oFile(j,i)->Blue = img_rock(x,y)->Blue;
			}
			// blend Textures
			if(height>181.0f && height<=201.0f)
			{
				float factorA = (height - 186.0f)/20.0f;
				float factorB = 1.0f - factorA;

				//Texture 1
				wrapAround(x,y,img_rock.TellWidth(),img_rock.TellHeight());

				oFile(j,i)->Red = img_rock(x,y)->Red * factorB;
				oFile(j,i)->Green = img_rock(x,y)->Green * factorB;
				oFile(j,i)->Blue = img_rock(x,y)->Blue * factorB;

				x=j;
				y=i;

				//Texture 2
				wrapAround(x,y,img_snow.TellWidth(),img_snow.TellHeight());

				oFile(j,i)->Red += img_snow(x,y)->Red * factorA;
				oFile(j,i)->Green += img_snow(x,y)->Green * factorA;
				oFile(j,i)->Blue += img_snow(x,y)->Blue * factorA;
			}
			if(height>201.0f)
			{
				wrapAround(x,y,img_snow.TellWidth(),img_snow.TellHeight());

				oFile(j,i)->Red = img_snow(x,y)->Red;
				oFile(j,i)->Green = img_snow(x,y)->Green;
				oFile(j,i)->Blue = img_snow(x,y)->Blue;
			}

			// range validation
			if(oFile(j,i)->Red<0) oFile(j,i)->Red = 0;
			if(oFile(j,i)->Red>255) oFile(j,i)->Red = 255;
			if(oFile(j,i)->Green<0) oFile(j,i)->Green = 0;
			if(oFile(j,i)->Green>255) oFile(j,i)->Green = 255;
			if(oFile(j,i)->Blue<0) oFile(j,i)->Blue = 0;
			if(oFile(j,i)->Blue>255) oFile(j,i)->Blue = 255;
		}
	}

	end = clock();
	std::cout << "Runtime (ProceduralTexture): " << start << " to " << end << " --> " << end-start << " / " << CLOCKS_PER_SEC << " seconds" << std::endl;

	oFile.WriteToFile(file.c_str());

}// makeProceduralTexture

void TerrainGenerator::wrapAround(int &i, int &j, int max_i, int max_j)
{
	if(i<0)
	{
		do 
		{
			i+=max_i;
		} while (i<0);
	}
	if(j<0)
	{
		do 
		{
			j+=max_j;
		} while (j<0);
	}
	if(i>=max_i)
	{
		do 
		{
			i-=max_i;
		} while (i>=max_i);
	}
	if(j>=max_j)
	{
		do 
		{
			j-=max_j;
		} while (j>=max_j);
	}
}

// --------------------------------------------------------------------------------------------------------------------------------
float TerrainGenerator::random(float a, float b)
{
	return (float)rand() / 32768.0f * (b - a) + a;
} // random

// --------------------------------------------------------------------------------------------------------------------------------
int TerrainGenerator::irand( int a, int b) {
	double r = b - a + 1;
	return a + static_cast<int>(r * rand()/(RAND_MAX+1.0));
} // irand

// --------------------------------------------------------------------------------------------------------------------------------
int TerrainGenerator::computeIndex(int a, int b)
{
	if(a<0 || b<0 || a>=m_iSize || a>= m_iSize)
	{
		std::cout << "range error in computeIndex\n";
		return 0;
	}
	else
		return a + b * m_iSize;
	//return ((((int)a)&m_iSizeMask) + (((int)b)&m_iSizeMask) * m_iSize);
} // computeIndex