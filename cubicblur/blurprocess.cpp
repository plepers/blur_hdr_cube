
#include "cubicblur.h"

typedef struct core_blur_params {
	float x;
	float y;
	float z;
	CubePixel ** faces;
	float power;
	float curve;
	int facesize;
	float* res;
} BLURPARAM;

float** computeBlur( CubePixel ** faces, int inputSize, int outputSize, float power, float curve );


void getBlur( void* in ) {
	
	BLURPARAM* params = (BLURPARAM*) in;

	float x = params->x;
	float y = params->y;
	float z = params->z;
	float power = params->power;
	float curve = params->curve;
	float mcurve = 1.0f-curve;
	float invCurve = 1.0f/curve;

	CubePixel ** faces = params->faces;


	// normalize
	float len = sqrt( x*x+y*y+z*z );
	x/=len;
	y/=len;
	z/=len;

	int f;
	int flen = params->facesize;
	flen = flen*flen;

	size_t psize;
	float dot, powdot; 
	double trR, trG, trB;
	double bias;
	float px, py, pz;
	
	CubePixel* face;

	psize = sizeof( CubePixel );
	
	bias = 0.0;
	trR = 0.0;
	trG = 0.0;
	trB = 0.0;
	

	for (int i = 0; i < 6; ++i)  {

		// compute face
		f = 0;
		face = faces[i];

		

		while( f < flen ) {
			px = face[f].x;
			py = face[f].y;
			pz = face[f].z;
			dot = ((x*px+y*py+z*pz)-mcurve)/invCurve;
			if( dot > 0.0 ) {
				powdot = pow( dot, power );
				trR += powdot * face[f].r;
				trG += powdot * face[f].g;
				trB += powdot * face[f].b;
				bias += powdot;
			}
			f++;
		}
		
		

	}

	
	trR /= bias;
	trG /= bias;
	trB /= bias;
	
	params->res[0] = float( trR );
	params->res[1] = float( trG );
	params->res[2] = float( trB );

}

float** blurFaces( float ** faces, int inputSize, int outputSize, float power, float curve ) {

	CubePixel** cfaces = createCpMap( faces, inputSize );


	float** blur = computeBlur( cfaces, inputSize, outputSize, power, curve );


	// free

	for (int i = 0; i < 6; ++i) 
		free( cfaces[i] );

	free( cfaces );

	return blur;

}


float** computeBlur( CubePixel ** faces, int inputSize, int outputSize, float power, float curve ) {

	int j;
	int y,x,c;
	float px, py, pz;

	float* bface;


	float** bfaces = (float **) malloc( 6 * sizeof( float* ) );

	float invSize = 1.0f / outputSize;
	float invSize2 = invSize * 2.0f;

	for (j = 0; j < 6; ++j) 
		bfaces[j] = (float *) malloc( outputSize*outputSize*3 * sizeof( float ) );

	HANDLE* handles = (HANDLE*) malloc( outputSize * sizeof( HANDLE ) );
	BLURPARAM* params = (BLURPARAM*) malloc( outputSize * sizeof( BLURPARAM ) );

	// --------------------------  left


	bface = bfaces[0];
	c = 0;

	for (y = 0; y < outputSize; ++y) {
		for (x = 0;x < outputSize; ++x) {

			px = -1.0f;
			py = -(y+.5f) * invSize2 + 1.0f;
			pz = (x+.5f) * invSize2 - 1.0f;
			
			params[x].x = px;
			params[x].y = py;
			params[x].z = pz;
			params[x].faces= faces;
			params[x].power= power;
			params[x].curve = curve;
			params[x].facesize= inputSize;
			params[x].res= &bface[c*3];

			handles[x] = (HANDLE) _beginthread( getBlur, 0, &params[x]  );
			c++;

		}
		WaitForMultipleObjects( outputSize, handles, true, INFINITE );
	}
	printf( "." );

	// --------------------------  right

	bface = bfaces[1];
	c = 0;

	for (y = 0; y < outputSize; ++y) {
		for (x = 0;x < outputSize; ++x) {

			px = 1.0f;
			py = -(y+.5f) * invSize2 + 1.0f;
			pz = -(x+.5f) * invSize2 + 1.0f;

			params[x].x = px;
			params[x].y = py;
			params[x].z = pz;
			params[x].faces= faces;
			params[x].power= power;
			params[x].curve = curve;
			params[x].facesize= inputSize;
			params[x].res= &bface[c*3];

			handles[x] = (HANDLE) _beginthread( getBlur, 0, &params[x]  );
			c++;
		}
		WaitForMultipleObjects( outputSize, handles, true, INFINITE );
	}
	printf( "." );

	
	// --------------------------  top

	bface = bfaces[2];
	c = 0;

	for (y = 0; y < outputSize; ++y) {
		for (x = 0;x < outputSize; ++x) {

			px = (x+.5f) * invSize2 - 1.0f;
			py = 1.0f;
			pz = (y+.5f) * invSize2 - 1.0f;

			params[x].x = px;
			params[x].y = py;
			params[x].z = pz;
			params[x].faces= faces;
			params[x].power= power;
			params[x].curve = curve;
			params[x].facesize= inputSize;
			params[x].res= &bface[c*3];

			handles[x] = (HANDLE) _beginthread( getBlur, 0, &params[x]  );
			c++;
		}
		WaitForMultipleObjects( outputSize, handles, true, INFINITE );
	}
	printf( "." );

	// --------------------------  bottom

	bface = bfaces[3];
	c = 0;

	for (y = 0; y < outputSize; ++y) {
		for (x = 0;x < outputSize; ++x) {

			px = (x+.5f) * invSize2 - 1.0f;
			py = -1.0f;
			pz = -(y+.5f) * invSize2 + 1.0f;

			params[x].x = px;
			params[x].y = py;
			params[x].z = pz;
			params[x].faces= faces;
			params[x].power= power;
			params[x].curve = curve;
			params[x].facesize= inputSize;
			params[x].res= &bface[c*3];

			handles[x] = (HANDLE) _beginthread( getBlur, 0, &params[x]  );
			c++;
		}
		WaitForMultipleObjects( outputSize, handles, true, INFINITE );
	}
	printf( "." );

	
	// --------------------------  back

	bface = bfaces[5];
	c = 0;

	for (y = 0; y < outputSize; ++y) {
		for (x = 0;x < outputSize; ++x) {

			px = -(x+.5f) * invSize2 + 1.0f;
			py = -(y+.5f) * invSize2 + 1.0f;
			pz = -1.0f;

			params[x].x = px;
			params[x].y = py;
			params[x].z = pz;
			params[x].faces= faces;
			params[x].power= power;
			params[x].curve = curve;
			params[x].facesize= inputSize;
			params[x].res= &bface[c*3];

			handles[x] = (HANDLE) _beginthread( getBlur, 0, &params[x]  );
			c++;
		}
		WaitForMultipleObjects( outputSize, handles, true, INFINITE );
	}
	printf( "." );



	// --------------------------  front

	bface = bfaces[4];
	c = 0;

	for (y = 0; y < outputSize; ++y) {
		for (x = 0;x < outputSize; ++x) {

			px = (x+.5f) * invSize2 - 1.0f;
			py = -(y+.5f) * invSize2 + 1.0f;
			pz = 1.0f;

			params[x].x = px;
			params[x].y = py;
			params[x].z = pz;
			params[x].faces= faces;
			params[x].power= power;
			params[x].curve = curve;
			params[x].facesize= inputSize;
			params[x].res= &bface[c*3];

			handles[x] = (HANDLE) _beginthread( getBlur, 0, &params[x]  );
			c++;
		}
		WaitForMultipleObjects( outputSize, handles, true, INFINITE );

		//printf( "frontface line %i \n", y );
	}
	printf( ".\n" );


	Sleep( 50 );
	free( params );
	free( handles );

	return bfaces;
}

CubePixel** createCpMap( float ** faces, int inputSize ) {

	int y;
	int x;
	int c;
	
	CubePixel* cp;
	CubePixel* cface;
	float* face;

	float px, py, pz, l;

	CubePixel** cfaces = (CubePixel**) malloc( 6*sizeof( CubePixel* ) );

	float invSize = 1.0f / inputSize;
	float invSize2 = invSize * 2;

	for (y = 0; y < 6; ++y) 
		cfaces[y] = (CubePixel*) malloc( ( inputSize*inputSize ) * sizeof( CubePixel ) );
	

	// --------------------------  left

	cface = cfaces[0];
	face = faces[0];
	c = 0;

	for (y = 0; y < inputSize; ++y) {
		for (x = 0;x < inputSize; ++x) {
			cp = &cface[c];

			px = -1.0f;
			py = -(y+.5f) * invSize2 + 1.0f;
			pz = (x+.5f) * invSize2 - 1.0f;

			l = sqrt( px*px+py*py+pz*pz );

			cp->r = face[c*3];
			cp->g = face[c*3+1];
			cp->b = face[c*3+2];

			cp->x = px/l;
			cp->y = py/l;
			cp->z = pz/l;

			c++;

		}
	}
	

	// --------------------------  right

	cface = cfaces[1];
	face = faces[1];
	c = 0;

	for (y = 0; y < inputSize; ++y) {
		for (x = 0;x < inputSize; ++x) {
			cp = &cface[c];

			px = 1.0f;
			py = -(y+.5f) * invSize2 + 1.0f;
			pz = -(x+.5f) * invSize2 + 1.0f;

			l = sqrt( px*px+py*py+pz*pz );

			cp->r = face[c*3];
			cp->g = face[c*3+1];
			cp->b = face[c*3+2];

			cp->x = px/l;
			cp->y = py/l;
			cp->z = pz/l;

			c++;
		}
	}

	
	// --------------------------  top

	cface = cfaces[2];
	face = faces[2];
	c = 0;

	for (y = 0; y < inputSize; ++y) {
		for (x = 0;x < inputSize; ++x) {
			cp = &cface[c];

			px = (x+.5f) * invSize2 - 1.0f;
			py = 1.0f;
			pz = (y+.5f) * invSize2 - 1.0f;

			l = sqrt( px*px+py*py+pz*pz );

			cp->r = face[c*3];
			cp->g = face[c*3+1];
			cp->b = face[c*3+2];

			cp->x = px/l;
			cp->y = py/l;
			cp->z = pz/l;

			c++;
		}
	}

	// --------------------------  bottom

	cface = cfaces[3];
	face = faces[3];
	c = 0;

	for (y = 0; y < inputSize; ++y) {
		for (x = 0;x < inputSize; ++x) {
			cp = &cface[c];

			px = (x+.5f) * invSize2 - 1.0f;
			py = -1.0f;
			pz = -(y+.5f) * invSize2 + 1.0f;

			l = sqrt( px*px+py*py+pz*pz );

			cp->r = face[c*3];
			cp->g = face[c*3+1];
			cp->b = face[c*3+2];

			cp->x = px/l;
			cp->y = py/l;
			cp->z = pz/l;

			c++;
		}
	}


	// --------------------------  front

	cface = cfaces[4];
	face = faces[4];
	c = 0;

	for (y = 0; y < inputSize; ++y) {
		for (x = 0;x < inputSize; ++x) {
			cp = &cface[c];

			px = (x+.5f) * invSize2 - 1.0f;
			py = -(y+.5f) * invSize2 + 1.0f;
			pz = 1.0f;

			l = sqrt( px*px+py*py+pz*pz );

			cp->r = face[c*3];
			cp->g = face[c*3+1];
			cp->b = face[c*3+2];

			cp->x = px/l;
			cp->y = py/l;
			cp->z = pz/l;

			c++;
		}
	}


	// --------------------------  back

	cface = cfaces[5];
	face = faces[5];
	c = 0;

	for (y = 0; y < inputSize; ++y) {
		for (x = 0;x < inputSize; ++x) {
			cp = &cface[c];

			px = -(x+.5f) * invSize2 + 1.0f;
			py = -(y+.5f) * invSize2 + 1.0f;
			pz = -1.0f;

			l = sqrt( px*px+py*py+pz*pz );

			cp->r = face[c*3];
			cp->g = face[c*3+1];
			cp->b = face[c*3+2];

			cp->x = px/l;
			cp->y = py/l;
			cp->z = pz/l;

			c++;
		}
	}

	return cfaces;

}