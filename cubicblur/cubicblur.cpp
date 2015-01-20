// cubicblur.cpp : Defines the entry point for the console application.
//

#include "cubicblur.h"

using namespace TCLAP;
using namespace std;


float** splitCrossmap( HDRLoaderResult* hdrRes );
float* assembleCrossmap( float** faces, int faceSize );
unsigned char* f32toRgbe( float* faces, int w, int h, double base );
float* subscale( float* pixels, int w, int h, int level );
float** subscaleFaces(float** faces, int w, int h, int level );

int numMipForSize( int inSize ) {
	
	int mip = 1;
	int size = inSize;
	while( size > 2 ) {
		size = size >> 1;
		mip++;
	}

	return mip;
		
}

int getMipForSize( int reqSize, int inSize ) {
	int nummips = numMipForSize( inSize );

	for( int m = 0; m<nummips;m++ ) {
		if(( inSize >> m ) <= reqSize ) return m;
	}
	return -1;
}

bool IsPowerOfTwo(unsigned int x)
{
    return (x != 0) && ((x & (x - 1)) == 0);
}

double getNewBase( char min, char max ) {
	double newbaseMax = pow( pow( 2.0, (double)max ), 1.0/128.0 );
	double newbaseMin = pow( pow( 2.0, (double)min ), -1.0/128.0 );
    
    if( newbaseMax > newbaseMin)
        return newbaseMax;
	return newbaseMin;
}


int main(int argc, char* argv[])
{

	int faceSize;

	CmdLine cmd("diffusefilter - generate blured cube map", ' ', "0.1");

	ValueArg<string> a_input("i","input","input raw cube texture",true,"","string", cmd );
	ValueArg<string> a_output("o","output","output raw file",true,"","string", cmd );
	
    ValueArg<float> a_power("p","power","fresnel power",false,1.0,"float", cmd );
    ValueArg<float> a_curve("c","curve","dot curve",false,1.0,"float", cmd );
    ValueArg<float> a_mix("m","mix","dot curve",false,0.0,"float", cmd );
	ValueArg<int> a_size("s","size","output size",false,64,"int", cmd );
	
	SwitchArg a_rebase("r","rebase", "maximize exponent range", true);
	cmd.add( a_rebase );

	cmd.parse( argc, argv );
	
	const char* input = a_input.getValue().c_str();
	const char* output = a_output.getValue().c_str();
	
    const float power = a_power.getValue();
    const float curve = a_curve.getValue();
    const float mix = a_mix.getValue();
	int size = a_size.getValue();
	const bool rebase = a_rebase.getValue();

	

	if( ! IsPowerOfTwo( size ) ) {
		printf( "error - output size %i must be POT", size );
		return 2;
	}
	
	printf( "blurring input cube texture \n    input %s \n    output %s \n" ,  input, output );


	//==================================================
	//								       Load hdr file
	//==================================================

	HDRLoaderResult* hdrData = new HDRLoaderResult();

	HDRLoader* hdrLoader = new HDRLoader();
	if( ! hdrLoader->load( input, *hdrData ) ) {
		printf( "error loading %s \n", input );
		return 1;
	}

	faceSize = hdrData->height/4;

	printf( "input loaded \n   size : %i*%i \n   range %i>%i \n", hdrData->width, hdrData->height , hdrData->eMin, hdrData->eMax );
	
	double base;
	if( rebase ) 
		base = getNewBase(  hdrData->eMin, hdrData->eMax );
	else
		base = 2.0f;

	//==================================================
	//								       extract faces
	//==================================================

	int mlevel = 1;

	printf( "splitCrossmap \n" );
	float** faces = splitCrossmap( hdrData );


	
	//==================================================
	//								       process blur
	//==================================================

	int rescaleMul = -1; // subscale source to have source double size of output
	
	int mip = getMipForSize( size, faceSize ) - rescaleMul;

	if( mip > 0 ) {
		float** sfaces = subscaleFaces( faces, faceSize, faceSize, mip );
		faceSize = faceSize >> mip;

		for (int j = 0; j < 6; ++j ) 
			free( faces[j] );
		free( faces );

		faces = sfaces;
	}

	float** blurred = blurFaces( faces, faceSize, size, power, curve, mix );
	
	
	//size =  hdrData->height/4;

	//==================================================
	//								       reassembling
	//==================================================

	printf( "assembleCrossmap \n" );
	float* assembled = assembleCrossmap( blurred, size );
	//float* assembled = assembleCrossmap( faces, size );
	
	printf( "f32toRgbe (base : %f)\n", base );
	const unsigned char* rgbe = f32toRgbe( assembled, size*3, size*4, base );

	printf( "encode png \n" );
	/*Encode the image*/
	unsigned error = lodepng_encode32_file( output, rgbe, size*3, size*4 );

	/*if there's an error, display it*/
	if(error) printf("error %u: %s \n", error, lodepng_error_text(error));


	delete hdrLoader;
	delete hdrData;

	free( faces );
	free( assembled );

	return 0;
}

float** subscaleFaces(float** faces, int w, int h, int level ) {
	float ** sfaces = (float **) malloc( 6 * sizeof( float* ) );

	for (int j = 0; j < 6; ++j) 
		sfaces[j] = subscale( faces[j], w, h, level );

	return sfaces;
}

float* subscale( float* pixels, int w, int h, int level ) {

	int rw = w >> level;
	int rh = h >> level;

	int sx, sy;

	int w3 = w*3;
	int rw3 = rw*3;

	int p;
	
	int scale = 1<<(level);
	int scale2 = scale*scale;

	int plen = w*h;

	float* res = (float*) malloc( plen * sizeof(float)*3  );
	
	double pbuffR;
	double pbuffG;
	double pbuffB;


	for( int y = 0; y< rh; y++ ) {
		for( int x = 0; x< rw; x++ ) {
			pbuffR = 0.0f;
			pbuffG = 0.0f;
			pbuffB = 0.0f;
			sx = x*scale;
			sy = y*scale;

			for( int i = sy; i< sy+scale; i++ ) {
				for( int j = sx; j< sx+scale; j++ ) {
					p = (i*w3 + j*3);
					pbuffR += pixels[ p ];
					pbuffG += pixels[ p + 1 ];
					pbuffB += pixels[ p + 2 ];
				}
			}
			p =  (y*rw3 + x*3);
			res[ p ] = pbuffR/scale2;
			res[ p + 1 ] = pbuffG/scale2;
			res[ p + 2 ] = pbuffB/scale2;
		}

	}

	return res;

}

unsigned char* f32toRgbe( float* pixels, int w, int h, double base ) {

	//base = 2.0;

	int j;

	int resSize = w*h;

	unsigned char* rgbe = ( unsigned char* ) malloc( resSize*4*sizeof( unsigned char* ) );

	float r, g, b; 
	double e, re;
	int f;

	double logbase = log( base ); 
	

	int c = 0;
	int fc = 0;
	for (j = 0; j < resSize; j++) {
		
		fc = j*3;
		c = j*4;

		r = pixels[fc];
		g = pixels[fc+1];
		b = pixels[fc+2];

		re = max( r, max( g, b ) );

		f = int( ceil( log( re ) / logbase ) );

		if( f < -128.0f ) f = -128.0f;
		if( f > 127.0f ) f = 127.0f;

		e = pow( base, f );
		
		r = r*255.0f / e;
		g = g*255.0f / e;
		b = b*255.0f / e;

		f += 128.0f;

		
		rgbe[c] = char( r );
		rgbe[c+1] = char( g );
		rgbe[c+2] = char( b );
		rgbe[c+3] = char( f );
	}


	return rgbe;
}

float* assembleCrossmap( float** faces, int faceSize) {

	int i;
	int j;
	int p;
	float* face;
	float* startPos;

	

	int fsize = sizeof( float );
	int facesize3 = faceSize*3;

	const int linelen = fsize*3*faceSize;
	
	int twidth = faceSize * 3;
	int theight = faceSize * 4;

	int resSize = faceSize*faceSize*12* 3;

	float* source = (float*) malloc( resSize * sizeof(float)  );

	for (j = 0; j < resSize; ++j) {
		source[j] = 0.0;
	}

	// left face

	face = faces[0];

	for (j = 0; j < faceSize; ++j) {
		startPos = &source[ twidth * (faceSize+j)*3 ];
		memcpy( startPos, face, linelen );
		face += facesize3;
	}

	// right face

	face = faces[1];

	for (j = 0; j < faceSize; ++j) {
		startPos = &source[ (twidth * (faceSize+j) + (2 * faceSize))*3 ];
		memcpy( startPos, face, linelen );
		face += facesize3;
	}

	// top face

	face = faces[2];

	for (j = 0; j < faceSize; ++j) {
		startPos = &source[ (twidth * j + (faceSize))*3 ];
		memcpy( startPos, face, linelen );
		face += facesize3;
	}

	// bottom face

	face = faces[3];

	for (j = 0; j < faceSize; ++j) {
		startPos = &source[ (twidth * (2*faceSize+j) + (faceSize))*3 ];
		memcpy( startPos, face, linelen );
		face += facesize3;
	}

	// front face

	face = faces[4];

	for (j = 0; j < faceSize; ++j) {
		startPos = &source[ (twidth * (faceSize+j) + (faceSize))*3 ];
		memcpy( startPos, face, linelen );
		face += facesize3;
	}

	// back face (vertical revert)

	face = faces[5];
	int c = 0;
	for (j = 0; j < faceSize; ++j) {
		for ( i = 0; i < faceSize; ++i) {
			p = (( twidth * (theight-j-1) )+(2 * faceSize - (i + 1) ))*3;
			source[ p ] = face[c++];
			source[ p+1 ] = face[c++];
			source[ p+2 ] = face[c++];
		}
	}

	return source;

}

float** splitCrossmap( HDRLoaderResult* hdrRes ){

	int i;
	int j;
	int p;
	int faceSize;
	int faceSize3;
	int lineByteslen;
	float* startPos;
	int fsize = sizeof( float );
	


	float* face;
	float ** faces;
	float* source;

	source = hdrRes->cols;

	faceSize = hdrRes->height/4;
	faceSize3 = faceSize*3;
	const int cpylen = faceSize*fsize*3;

	faces = (float **) malloc( 6 * sizeof( float* ) );

	lineByteslen = hdrRes->width * 12;
	
	for (j = 0; j < 6; ++j) 
		faces[j] = (float *) malloc( faceSize*faceSize*3 * sizeof( float ) );
	
	
	
	// left face

	face = faces[0];

	for (j = 0; j < faceSize; ++j) {
		startPos = &source[ (hdrRes->width * (faceSize+j))*3 ];
		memcpy( face, startPos, cpylen );
		face += faceSize3;
	}

	// right face

	face = faces[1];

	for (j = 0; j < faceSize; ++j) {
		startPos = &source[ (hdrRes->width * (faceSize+j) + (2 * faceSize))*3 ];
		memcpy( face, startPos, cpylen );
		face += faceSize3;
	}

	// top face

	face = faces[2];

	for (j = 0; j < faceSize; ++j) {
		startPos = &source[ (hdrRes->width * j + (faceSize))*3 ];
		memcpy( face, startPos, cpylen );
		face += faceSize3;
	}

	// bottom face

	face = faces[3];

	for (j = 0; j < faceSize; ++j) {
		startPos = &source[ (hdrRes->width * (2*faceSize+j) + (faceSize))*3 ];
		memcpy( face, startPos, cpylen );
		face += faceSize3;
	}

	// front face

	face = faces[4];

	for (j = 0; j < faceSize; ++j) {
		startPos = &source[ (hdrRes->width * (faceSize+j) + (faceSize))*3 ];
		memcpy( face, startPos, cpylen );
		face += faceSize3;
	}

	// back face (vertical revert)

	face = faces[5];
	int c = 0;
	for (j = 0; j < faceSize; ++j) {
		for ( i = 0; i < faceSize; ++i) {
			p = (( hdrRes->width * (hdrRes->height-j-1) )+(2 * faceSize - (i + 1) ))*3;
			face[c++] = source[ p ];
			face[c++] = source[ p+1 ];
			face[c++] = source[ p+2 ];
		}
	}


	return faces;

}


