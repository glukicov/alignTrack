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
    dispX = 0.01;  // mnaul displacemnt by 0.01 cm [0.1]
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

   // Track parameters for rand-generated line MC [start and end positions outside of detectors]
    rand_num = (( RandomBuffer::instance()->get_uniform_number()+ RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max()));
    float x_0 = beamPositionLength * (rand_num); //uniform vertex
    rand_num = ((RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max()));
    float x_1 = beamPositionLength * (rand_num); //uniform exit point: so fitting a line to these two points
    //float x_slope=(x_1-x_0)/distance[layerN];
    float x_slope=(x_1-x_0)/(beamStop-beamStart); //beams starts at z=0, beamsStop=25; 
    //XXX Wouldn't the correct detention be: float x_slope=(beamStop-beamStart)/(x_1-x_0); m=dz/dx? 

    line.x_gen.push_back(x_0);
    line.z_gen.push_back(beamStart);
    line.x_gen.push_back(x_1);
    line.z_gen.push_back(beamStop);
    
    float x = x_0;
    float dx = x_slope;
    float previousPosition = 0.0;  // previous position in "z"
   
   //The main loop is for modules [which are misaligned]
    //TODO within modules we will have loops over views, layers, straws 
    for(int i=0; i<moduleN; i++){
        //distance between consecutive layers 
        float dPostion = distance[i] - previousPosition;
        previousPosition = distance[i];

        //positions on the detector plane  
        float xs=x_0 + distance[i] * x_slope;
        line.x_det.push_back(xs);
        
        //true track position [for MS]
        x=x+dx*dPostion;
        line.x_true.push_back(x);

        //multiple scattering
        rand_gaus= RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
        dx = dx+ rand_gaus * scatterError;
                
        // which straw was hit [0,5 Y; 0,10 X] MC rejection if beyond plane geometry
         //Rejection of hits
        //TODO rewrite this ugly part 
        // For now, NO hit rejection 
        int imx=1;
        if (imx < 0 || imx >= strawXN){ //[between (0,10]       
            continue;
        } 
        
        //Module number
        int ihit= i; // which module are we on
        line.i_hits.push_back(ihit); // vector of planes that were actually hit

        float xl=x-sdevX[layer[i]]; //misalign. 
        line.x_mis.push_back(xl);
       
        // we seem to now redefine the coordinates so that x is now the distance and y is a measure of the residual
        line.z_hits.push_back(distance[i]);
        // the residual looks to be deltaX + deltaY rather than the magnitude of the distance... worth noting?
        // projection Y is always 0 for non-stereo modules?? what is the motivation? 
        rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
        float hit = (xl-xs)*projectionX[i] + rand_gaus *resolution; 
        line.x_hits.push_back(hit);
        line.hit_sigmas.push_back(resolution);
        line.hit_count++;
   
    }// end of looping over detector layers
    return line; // Return data from simulated track
    
} // end of MC

//Geometry of detector arrangement 
void Tracker::setGeometry(ofstream& debug_geom, bool debugBool){
    float seperation=startingDistanceModule0;
    // Geometry of detecor arrangement 
    for (int layer_i=0; layer_i<moduleN; layer_i++){
        layer.push_back(layer_i);  // layer [starting from 1st layer]
        distance.push_back(seperation);  //distance between planes  [14]
        resolutionLayer.push_back(resolution); //resolution
        projectionX.push_back(float(1.0));  // x
        seperation+=moduleSpacing;  // incrementing distance between detectors
    }  // end of looping over layers
} // end of geom

// MC misalignment of detectors 
void Tracker::misalign(ofstream& debug_mis, bool debugBool){
    float rand_gaus;
   //Now misaligning detectors
    for (int i=0; i<moduleN; i++){   
                rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
                sdevX[i] = dispX * rand_gaus;
    } // end of modules
}//end of misalign


/**
   Write a constraint file to the supplied file-stream.

    @param constraint_file Reference to ofstream to write constraint file to. 
 */
void Tracker::write_constraint_file(ofstream& constraint_file, ofstream& debug_con, bool debugBool) {
    // constraints: fix center straws in first/last layer
    // Check constraints file is open, then write. 
    if (constraint_file.is_open()) {
       
        //Evaluation of constraints TODO ??? 
        float one = 1.0;
        
        for (int i = 0; i < moduleN; i++){ 
            constraint_file << "Constraint 0.0" << endl;
                int labelt=i;
                constraint_file << labelt << " " << fixed << setprecision(5) << one<< endl;
                //sdevX[i]=0.0;     // fix center modules at 0. XXX
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