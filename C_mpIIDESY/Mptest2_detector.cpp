/*

This source file contains definitions of various functions in the Detector class. 

*/  

#include "Mptest2_detector.h"

using namespace std;

// Set up empty pointer for instance of class.
Detector* Detector::s_instance = NULL; 

/** 
	Empty constructor for detector class.
 */
Detector::Detector() {
    resolution = 0.002;

}

/** 
	Empty destructor for detector class.
*/
Detector::~Detector() {
}


/**
   Get pointer to only instance of detector class, creating this instance if it doesn't already exist.

   @return Pointer to detector instance.
 */
Detector* Detector::instance() {

	// Create pointer to class instance if one doesn't exist already, then return that pointer.
	if (s_instance == NULL) s_instance = new Detector();
	return s_instance;
}

/**
   Simulate the passage of a randomly generated linear track through the detector, calculating properties of hits by the track on detector planes.
  
   @return LineData struct containing data about detector hits.
*/
//Source code for genlin courtesy of John. 
// Function to simulate a linear track through the detector, returning data about detector hits.

LineData Detector::genlin2(float scatterError, ofstream& debug_calc, ofstream& debug_off, ofstream& debug_mc, bool debugBool) {
  
	// Set up new container for track data, with hit count set to zero`
	LineData line;
    line.hit_count = 0;

    float rand_num, rand_num1, rand_num2, rand_num3, rand_num4;
    float rand_gaus;

    // Track parameters for rand-generated line MC 
    rand_num1 = (( RandomBuffer::instance()->get_uniform_number()+ RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max())); //- 0.5;
    rand_num1 -= 0.5;
    float x_0 = layerSize * (rand_num1); //uniform vertex
    rand_num2 = ((RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max())); //- 0.5;
    rand_num2 -= 0.5;
    float y_0 = layerSize * (rand_num2); //uniform vertex
    rand_num3 = ((RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max())); //- 0.5;
    rand_num3 -= 0.5;
    float x_1 = layerSize * (rand_num3); //uniform exit point: so fitting a line to these two points
    rand_num4 = ((RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max())); //- 0.5;
    float extraRand = rand_num4;
    float y_1 = layerSize * (extraRand - 0.5); //uniform exit point:
    float x_slope=(x_1-x_0)/distance[layerN-1];
    float y_slope=(y_1-y_0)/distance[layerN-1];
    

    float x = x_0;
    float dx = x_slope;
    float y = y_0;
    float dy = y_slope;
    float s_old = 0.0;  // previous position in "z"

    if (debugBool){
        debug_calc << x_0<< "  " << y_0 << " " << x_1<< " " << y_1 <<" " << x_slope<< " "<< y_slope << "  " << (rand_num1) << " " << (rand_num2) << " " << (rand_num3) << " " << (rand_num4) << " " << layerSize << " " << distance[layerN-1] << endl;
    }

    
    for(int i=0; i<layerN; i++){
        //distance between cons. layers 
        float ds = distance[i] - s_old;
        s_old = distance[i];

        //positions on the detector plane  
        float xs=x_0 + distance[i] * x_slope;
        float ys=y_0 + distance[i] * y_slope;
        line.x_det.push_back(xs);
        line.y_det.push_back(ys);
        
        //true track position [for MS]
        x=x+dx*ds;
        y=y+dy*ds;
        line.x_true.push_back(x);
        line.y_true.push_back(y);

        //multiple scattering
        rand_gaus= RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
        dx = dx+ rand_gaus * scatterError;
       
        rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
        dy = dy+ rand_gaus * scatterError;
        
        // which pixel was hit [0,5 Y; 0,10 X] MC rejection if beyond plane geometry
        int imx=int((x+layerSize*0.5)/layerSize*float(pixelXN));
        if (imx < 0 || imx >= pixelXN){ //[between (0,10]       
            continue;
        } 
        int imy=int((y+layerSize*0.5)/layerSize*float(pixelYN));
         if (imy < 0 || imy >= pixelYN){  //[between (0,5] 
            continue;
        }
       
        int ihit= ((i)*pixelYN+imy)*pixelXN+imx; // i from 0 to 13 (incl.)
        line.i_hits.push_back(ihit); // vector of planes that were actually hit

        float xl=x-sdevX[layer[i]-1][imy][imx]; //misalign. 
        float yl=y-sdevY[layer[i]-1][imy][imx];
        line.x_mis.push_back(xl);
        line.y_mis.push_back(yl);

        debug_mc << xl << " " << yl << " " << xs << " "  << ys << " " << dx << " " << dy << " " << x << " " << y << endl;  

       
        // we seem to now redefine the coordinates so that x is now the distance and y is a measure of the residual
        line.x_hits.push_back(distance[i]);
        // the residual looks to be deltaX + deltaY rather than the magnitude of the distance... worth noting?
        // projection Y is always 0 for non-stero modules?? what is the motivation? 
        rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
        float yhit = (xl-xs)*projectionX[i]+(yl-ys)*projectionY[i]+ rand_gaus *resolution; 
        line.y_hits.push_back(yhit);
        line.hit_sigmas.push_back(resolution);
        line.hit_count++;

        if (debugBool){ 
            debug_off << projectionX[i] << " " << projectionY[i] << " " << rand_gaus << " " << resolution << " " << sdevX[layer[i]-1][imy][imx] <<  endl;
        } 
   
    }// end of looping over detector layers
    
     
   
    return line; // Return data from simulated track
    
} // end of genlin2

//Geometry of detecor arrangement 
void Detector::setGeometry(ofstream& debug_geom, bool debugBool){
	float s=startingDistancePlane1;
    int i_counter = 0;
    float sign = 1.0;

    // Geometry of detecor arrangement 
    for (int layer_i=1; layer_i<=10; layer_i++){
        
        layer.push_back(layer_i);  // layer [starting from 1st layer]
        distance.push_back(s);  //distance between planes  [14]
        resolutionLayer.push_back(resolution); //resolution
        projectionX.push_back(float(1.0));  // x
        projectionY.push_back(float(0.0));  // y
        if (debugBool){
              // debug_geom << "layer_i= " << layer_i << " layer[i_counter]= " << layer[i_counter]  << endl;
              // debug_geom << "i_counter= " << i_counter << " distance[i_counter] " << distance[i_counter] << endl;
              // debug_geom << "projectionX= " << projectionX[i_counter] << " projectionY= " << projectionY[i_counter] << endl;
              // debug_geom << endl; 
        } 
        i_counter++; 
        //taking care of stereo planes [have no pixels] 1, 4, 7, 10
        if (((layer_i) % 3) == 1){
            layer.push_back(layer_i);  // layer
            distance.push_back(s+offset);  //distance between planes  [14]
            resolutionLayer.push_back(resolution); //resolution
            projectionX.push_back(std::sqrt(1.0-stereoTheta*stereoTheta));  // x
            projectionY.push_back(stereoTheta*sign);  // y
            if (debugBool){
              // debug_geom << "S_layer_i= " << layer_i << " layer[i_counter]= " << layer[i_counter]  << endl;
              // debug_geom << "i_counter= " << i_counter << " distance[i_counter] " << distance[i_counter] << endl;
              // debug_geom << "projectionX= " << projectionX[i_counter] << " projectionY= " << projectionY[i_counter] << endl;
              // debug_geom << endl; 
                debug_geom << std::sqrt(1.0-stereoTheta*stereoTheta) << endl;
            }  
            sign=-sign;
            i_counter++;    
        }

        s=s+planeDistance;  // incrimenting distance between detecors

    }  // end of looping over layers

if (debugBool){
                for (int i=0; i<14; i++){

                debug_geom << projectionX[i] << endl;
            }
    }  

} // end of geom



// MC misalignment of detecors 
void Detector::misalign(ofstream& debug_mis, bool debugBool){

    float rand_gaus;

    //bool debugBool = true; // print out to debug files 
    int counterMis = 0;
	//Now misaligning detecors
    float dispX = 0.01; // plane displacement in X .05 mm * N(0,1)
    float dispY = 0.01; // plane displacement in Y .05 mm * N(0,1)

    
    for (int i=0; i<=detectorN-1; i++){   
        for(int k=0; k<=pixelYN-1; k++){
            for(int l=0; l<=pixelXN-1; l++){
                
                rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
                sdevX[i][k][l] = dispX * rand_gaus;
                if (debugBool){debug_mis << sdevX[i][k][l] << "  " << dispX <<  " " << rand_gaus << endl;}  
                rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
                sdevY[i][k][l] = dispY * rand_gaus;
                counterMis++;

                
            
            } // // end of number of pixel in x
        } // end of number of pixels in y 
    } // end of layers
if (debugBool){
debug_mis << "counterMis= " << counterMis; 
}
}//end of misalign


/**
   Write a constraint file to the supplied file-stream.

    @param constraint_file Reference to ofstream to write constraint file to. 
 */
void Detector::write_constraint_file(ofstream& constraint_file, ofstream& debug_con, bool debugBool) {

    // constraints: fix center pixels in first/last layer
	// Check constraints file is open, then write. 
	if (constraint_file.is_open()) {
		
		cout << "Writing constraint file..." << endl;
		

		//Evaluation of constraints
    	int ncx = (pixelXN+1)/2; 
    	int lunt = 9;
    	float one = 1.0;
         
        //XXX looks like plane 1 and 10 are being fixed for the central x-row of 5 pixels
        //This then forms contraints on alignment (?)
		for (int i = 1; i <= detectorN; i=i+(detectorN-1)){ 
            constraint_file << "Constraint 0.0" << endl;
            for (int k=0; k<=pixelYN-1; k++){
                int labelt=(i*pixelYN+k)*pixelXN+ncx-1;
                constraint_file << labelt << " " << fixed << setprecision(5) << one<< endl;
              //sdevX[detectorN][pixelYN][pixelXN]=0.0;      // fix center pixels at 0.  ///TODO fix tis!!!!
                sdevX[i-1][k][ncx-1]=0.0;
                if (debugBool){
                    debug_con << sdevX[i-1][ncx-1][k] << " " << i << " "  <<  k << " " << ncx <<endl;
                }
            } // end of y loop
            constraint_file << "Constraint 0.0" << endl;
            for(int k=0; k<=pixelYN-1; k++){
                int labelt=(i*pixelYN+k)*pixelXN+ncx+1000-1;
                constraint_file << labelt << " " << fixed << setprecision(5) << one<< endl;
                //sdevY[detectorN][pixelYN][pixelXN]=0.0;                         ///TODO fix tis!!!!
                sdevY[i-1][k][ncx-1]=0.0; 

            } // end of x loop
        } // end of detecors loop 

	} // constraing file open 
} // end of wrting cons file 

/**
   Set filename to read uniform random numbers from.

   @param uniform_filename String for filename of file of uniform random numbers.
 */
void Detector::set_uniform_file(string uniform_filename) {
	RandomBuffer::instance()->open_uniform_file(uniform_filename);
}


/**
   Set filename to read gaussian random numbers from.

   @param uniform_filename String for filename of file of gaussian random numbers.
 */
void Detector::set_gaussian_file(string gaussian_filename) {
	RandomBuffer::instance()->open_gaussian_file(gaussian_filename);
}







