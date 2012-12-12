/*
Project code from a computer science graphics course.  Algorithms process images in the 
Targa file format to perform quantization and painterly rendering.  Code excerpt from
TargaImage.cpp file.

Samuel Heinith
October 2012
*/




///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
	unsigned char *rgb = To_RGB_All();
	typedef tuple<int, int, int> triplet;
	int factor = 8; //256/32 levels = 8
	int i, j;
	for (i=0; i<width*height*4; i+=4) {
		for (int c=0; c<3; ++c) {
			data[i+c] = (rgb[i+c]/factor)*factor;
		}
	}
	
	//Find 256 most popular colors
	map < triplet, int > c_counts;
	triplet key;
	//Populate map with color counts use RGB to encode color channels to integer key
	for (i=0; i<width*height*4; i+=4) {
		key = make_tuple(data[i], data[i+1], data[i+2]);
		if(c_counts.find(key) == c_counts.end()) { //Not already in map
			c_counts[triplet(key)] = 1; //Initialize new color
		} else {
			c_counts[triplet(key)]++;
		}
	}

	if(c_counts.size() > 256) { //If there are fewer than 256 colors then finished

		vector < pair<triplet, int> > colors(c_counts.begin(), c_counts.end()); //Copy map key/value pairs to a vector
		sort(colors.begin(), colors.end(), cmp); //Sort the vector by the value (color count) of the pair
		vector <triplet> reduced_colors(256);

		int red, green, blue;
		for(i=0; i<256; ++i) {
			tie (red, green, blue) = colors[i].first;
			reduced_colors[i] = make_tuple (red, green, blue);
			//cout << i << " " << red << " " << green << " " << blue << endl;
		}
		
		int r0, g0, b0, r1, g1, b1;
		double min_d, d;
		triplet closest;
		map <triplet, triplet> q_map;

		for (i = 0; i < width*height*4; i+=4) { //for each pixel
			key = make_tuple(data[i], data[i+1], data[i+2]);
			tie (r0, g0, b0) = key;
			if( q_map.find(key) == q_map.end() ) { //Search for correct target color

				//Not in map, find closest of the 256
				min_d = INT_MAX; //Reset minimum distance for each source color
				for (j = 0; j < 256; ++j) { //for each reduced color
					tie (r1, g1, b1) = reduced_colors[j]; //unpack the tuple
					d = sqrt( (double) (r1-r0)*(r1-r0) + (g1-g0)*(g1-g0) + (b1-b0)*(b1-b0) ); //Calculate euclidean distance
					if(d < min_d) {
						min_d = d;
						closest = reduced_colors[j];
					}
				}				
				q_map[triplet(key)] = closest; //Store closest color to the quantization map
			}
			tie (data[i], data[i+1], data[i+2]) = q_map[triplet(key)]; 	//Overwrite the pixel color with the quantized color
		}//for each pixel
		

	}//If fewer than 256


	delete rgb;
    return true;
}// Quant_Populosity





///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
	assert(N > 1);
	//assert(N < 15); //Overflows value
	unsigned char *rgb = To_RGB_All();

	int n = (int) N;
	int i, j, b, c;
	//Overwrite original image so that unfiltered edge pixels have their alpha multiplied out
	memcpy(data, rgb, sizeof(unsigned char) * width * height * 4);

	vector < int > coeffs(n); //Binomial coefficients
	for(i = 0; i<=n/2; ++i) {
		b = Binomial((n-1),i);
		coeffs[i] = b; //Fill first half
		coeffs[n - i - 1] = b; //Fill second half
	}

	vector< vector<int> > filter(n, vector<int>(n)); //nxn gaussian filter
	int x, y;
	for(y=0; y<n; ++y) {
		for(x=0; x<n; ++x) {
			filter[y][x] = coeffs[x]*coeffs[y];
		}
	}
	
	double scalar = 1/pow(2.0, 2*(n-1)); //Used to keep sum of filter elements equal to 1

	//Apply the filter to the image, read from rgb[], write to data[]
	int offset = n/2; //Closest to the edge the filter can be centered
	unsigned __int64 value; //Large type to prevent common overflows

	for (c=0; c<3; ++c) { //For each color channel
		for (y = offset; y < (height - offset); ++y) {
			for (x = offset; x < (width - offset); ++x) { //For each point in the image

				value = 0; //Reset for each filtering
				for (j=0; j<n; ++j) {  //For each value in the filter
					for (i=0; i<n; ++i) {
						value += rgb[((y - offset+j)*width +(x - offset+i))*4 + c] * filter[j][i];
					}//i
				}//j

				//Assign value to result image
				data[(y*width + x)*4 + c] = (int) floor(scalar*value + 0.5);
			}//x
		}//y
	}//c

	delete rgb;
	return true;
}// Filter_Gaussian_N





///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
//      Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
	TargaImage ref_image = TargaImage(*this);
	TargaImage canvas = TargaImage(width, height); //Constant color that will get painted on
	
	int x, y, i, j, g, m;
	int r0, g0, b0, r1, g1, b1;
	typedef tuple<int, int, double> error_point;

	double d, area_error;
	int T = 25; //Error threshold
	static const unsigned int arr[] = {7, 3, 1}; //Brush sizes
	vector<int> rs (arr, arr + sizeof(arr) / sizeof(arr[0]) );

	int rad;
	//*****
	//Paint
	//*****
	for(size_t b=0; b < rs.size(); b++) { //For each brush size
		ref_image.Filter_Gaussian_N(2*rs[b]+1); //referenceImage
		g = 2*rs[b]+1; //Could also be 2*rs[b]+1
		rad = 3*rs[b];
		//***Paint a layer***
		vector <Stroke> strokes;
		//Find all of the grid block centers
		vector < tuple<int, int> > centers;
		centers.reserve( (width/g + 1)* (height/g + 1) );
		for( x = g/2; x < g*(width/g); x+=g ) {
			for( y = g/2; y < g*(height/g); y+=g ) {
				centers.push_back( make_tuple(x, y) );
			}
		}
		//Construct overlapping blocks around right and bottom edge, Possibly Unnecessary
		if (width % g) { //Hold x constant equal to width-g/2-1, loop down the column increasing y			
			x = width - g/2 - 1;
			for( y = g/2; y < g*(height/g); y+=g ) {
				centers.push_back( make_tuple(x, y) );
			}
		}
		if (height % g) {
			y = height - g/2 - 1; //Hold y constant equal to height-g/2-1, loop across the row increasing x
			for( x = g/2; x < g*(width/g); x+=g ) {
				centers.push_back( make_tuple(x, y) );
			}
		}
		//if both x_extra and y_extra then need final bottom right corner block, probably unnecessary		
		

		for( size_t c=0; c < centers.size(); c++ ) { //For each grid center
			tie (x, y) = centers[c];
				//Find the error in the area/region/block
				vector < error_point > errors;
				errors.reserve(g*g);
				area_error = 0;
				for (j = -g/2; j<=g/2; ++j) {
					for (i = -g/2; i<=g/2; ++i) {
						m = (width*(y+j) + (x+i))*4; //Location of current pixel
						tie (r0, g0, b0) = make_tuple(ref_image.data[m], ref_image.data[m+1], ref_image.data[m+2]);
						tie (r1, g1, b1) = make_tuple(canvas.data[m], canvas.data[m+1], canvas.data[m+2]);
						//Calculate euclidean distance
						d = sqrt( (double) (r1-r0)*(r1-r0) + (g1-g0)*(g1-g0) + (b1-b0)*(b1-b0) ); 
						//Save the coordinate and its error
						errors.push_back(error_point(x+i, y+j, d));
						area_error += d;
					}//i
				}//j
				if( area_error > T ) {
					//Find largest error point
					sort(errors.begin(), errors.end(), cmp_errors);					
					tie (i, j, d) = errors.front();
					//Build a Stroke at loacation x,y with radius and value of r,g,b,alpha
					m = (width*j + i)*4; //Location of current pixel
					Stroke s = Stroke(rad, i, j, ref_image.data[m], ref_image.data[m+1], ref_image.data[m+2], ref_image.data[m+3]);
					strokes.push_back(s);
				}

		}//c

		//*******
		//Paint all strokes in random order onto canvas
		//*******
		random_shuffle ( strokes.begin(), strokes.end() );
		vector <Stroke>::iterator it;
		for (it = strokes.begin(); it != strokes.end(); it++) {
			canvas.Paint_Stroke(*it);
		}

		//Reset the reference image to the unchanged image data for the next smaller pass
		memcpy(ref_image.data, data, sizeof(unsigned char) * width * height * 4); 
	}//b	

	//Write the canvas data to the displayed image
	memcpy(data, canvas.data, sizeof(unsigned char) * width * height * 4);
	//Overwrite alpha values
	for(i=3; i<width*height*4; i+=4) {
		data[i] = 255;
	}

    return true;
}

