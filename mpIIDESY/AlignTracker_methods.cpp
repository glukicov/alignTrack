/*
This source file contains definitions of various functions in the method class. 
*/  

#include "AlignTracker_methods.h"

using namespace std;

// Set up empty pointer for instance of the class.
Tracker* Tracker::s_instance = NULL; 

/** 
Empty constructor for detector class.
 */
Tracker::Tracker() {
  //XXX non-static variables definition here 
    resolution=0.002;  // 20um = 0.002 cm setting this in the constructor //TODO decide if this is a constant 
    dispX = 0.01;
}

/** 
    Empty destructor for detector class.
*/
Tracker::~Tracker() {
}


/**
   Get pointer to only instance of detector class, creating this instance if it doesn't already exist.

   @return Pointer to detector instance.
 */
Tracker* Tracker::instance() {

    // Create pointer to class instance if one doesn't exist already, then return that pointer.
    if (s_instance == NULL) s_instance = new Tracker();
    return s_instance;
}

/**
   Simulate the passage of a randomly generated linear track through the detector, calculating properties of hits by the track on detector planes.
  
   @return LineData struct containing data about detector hits.
*/
LineData Tracker::MC(float scatterError, ofstream& debug_calc, ofstream& debug_off, ofstream& debug_mc, bool debugBool) {
  
    // Set up new container for track data, with hit count set to zero
    LineData line;
    line.hit_count = 0; //count only counts through detector

    float rand_num;
    float rand_gaus;

   // Track parameters for rand-generated line MC 
    rand_num = (( RandomBuffer::instance()->get_uniform_number()+ RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max()))- 0.5;
    float x_0 = layerSize * (rand_num); //uniform vertex
    rand_num = ((RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max()))- 0.5;
    float x_1 = layerSize * (rand_num); //uniform exit point: so fitting a line to these two points
    float x_slope=(x_1-x_0)/distance[layerN-1];
    
    float x = x_0;
    float dx = x_slope;
    float s_old = 0.0;  // previous position in "z"
   
    for(int i=0; i<layerN; i++){
        //distance between consecutive layers 
        float ds = distance[i] - s_old;
        s_old = distance[i];

        //positions on the detector plane  
        float xs=x_0 + distance[i] * x_slope;
        line.x_det.push_back(xs);
        
        //true track position [for MS]
        x=x+dx*ds;
        line.x_true.push_back(x);

        //multiple scattering
        rand_gaus= RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
        dx = dx+ rand_gaus * scatterError;
                
        // which straw was hit [0,5 Y; 0,10 X] MC rejection if beyond plane geometry
         //Rejection of hits
        //TODO rewrite this ugly part 
        int imx=int((x+layerSize*0.5)/layerSize*float(pixelXN));

        if (imx < 0 || imx >= pixelXN){ //[between (0,10]       
            continue;
        } 
        
        //This will be straw number 
        int ihit= (i)*pixelXN+imx; // i from 0 to 15 (incl.)
        line.i_hits.push_back(ihit); // vector of planes that were actually hit

        float xl=x-sdevX[layer[i]-1][imx]; //misalign. 
        line.x_mis.push_back(xl);
       
        // we seem to now redefine the coordinates so that x is now the distance and y is a measure of the residual
        line.x_hits.push_back(distance[i]);
        // the residual looks to be deltaX + deltaY rather than the magnitude of the distance... worth noting?
        // projection Y is always 0 for non-stereo modules?? what is the motivation? 
        rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
        float hit = (xl-xs)*projectionX[i] + rand_gaus *resolution; 
        line.y_hits.push_back(hit);
        line.hit_sigmas.push_back(resolution);
        line.hit_count++;
   
    }// end of looping over detector layers
    return line; // Return data from simulated track
    
} // end of MC

//Geometry of detector arrangement 
void Tracker::setGeometry(ofstream& debug_geom, bool debugBool){
    float s=startingDistancePlane1;
    int i_counter = 0;
    float sign = 1.0;

    // Geometry of detecor arrangement 
    for (int layer_i=1; layer_i<=10; layer_i++){
        layer.push_back(layer_i);  // layer [starting from 1st layer]
        distance.push_back(s);  //distance between planes  [14]
        resolutionLayer.push_back(resolution); //resolution
        projectionX.push_back(float(1.0));  // x
        i_counter++; 
        //taking care of stereo planes [have no pixels] 1, 4, 7, 10
        if (((layer_i) % 3) == 1){
            layer.push_back(layer_i);  // layer
            distance.push_back(s+offset);  //distance between planes  [14]
            resolutionLayer.push_back(resolution); //resolution
            projectionX.push_back(sqrt(1.0-stereoTheta*stereoTheta));  // x 
            sign=-sign;
            i_counter++;    
        }
        s=s+planeDistance;  // incrementing distance between detectors
    }  // end of looping over layers
} // end of geom

// MC misalignment of detectors 
void Tracker::misalign(ofstream& debug_mis, bool debugBool){
    float rand_gaus;
    int counterMis = 0;
    //Now misaligning detectors
    for (int i=0; i<=detectorN-1; i++){   
            for(int l=0; l<=pixelXN-1; l++){
                rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
                sdevX[i][l] = dispX * rand_gaus;
                counterMis++;
            } // // end of number of pixel in x
    } // end of layers
}//end of misalign


/**
   Write a constraint file to the supplied file-stream.

    @param constraint_file Reference to ofstream to write constraint file to. 
 */
void Tracker::write_constraint_file(ofstream& constraint_file, ofstream& debug_con, bool debugBool) {
    // constraints: fix center pixels in first/last layer
    // Check constraints file is open, then write. 
    if (constraint_file.is_open()) {
       
        //Evaluation of constraints
        int ncx = (pixelXN+1)/2; 
        int lunt = 9;
        float one = 1.0;
        
        for (int i = 1; i <= detectorN; i=i+(detectorN-1)){ 
            constraint_file << "Constraint 0.0" << endl;
                int labelt=(i)*pixelXN+ncx-1;
                constraint_file << labelt << " " << fixed << setprecision(5) << one<< endl;
                sdevX[i-1][ncx-1]=0.0;     // fix center pixels at 0.
            constraint_file << "Constraint 0.0" << endl;
        } // end of detectors loop 
    } // constrain file open 
} // end of writing cons file 

/**
   Set filename to read uniform random numbers from.
   @param uniform_filename String for filename of file of uniform random numbers.
 */
void Tracker::set_uniform_file(string uniform_filename) {
    RandomBuffer::instance()->open_uniform_file(uniform_filename);
}

/**
   Set filename to read Gaussian random numbers from.
   @param uniform_filename String for filename of file of Gaussian random numbers.
 */
void Tracker::set_gaussian_file(string gaussian_filename) {
    RandomBuffer::instance()->open_gaussian_file(gaussian_filename);
}