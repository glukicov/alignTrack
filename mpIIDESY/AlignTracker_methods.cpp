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
    resolution=0.01;  // 100um = 0.01 cm setting this in the constructor //TODO decide if this is a constant 
    dispX = 0.15;  // manual displacement by 0.15 cm 
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
  Distance (parallel) of closest approach between a line (track) and a point (straw/wire).
  //BIG TODO: Implement this properly 

  Equation of a line: z=mx+c
  The (x,z) coordinate of a straw are known. 

    1) Find x_hit [=(z-c)/m] for the line at given z
    2) Find straw closest to x_hit
    3) Find x_straw-x_hit 
   @return dca
*/
float Tracker::DCA(std::vector<float> layer_x, float z_straw, float slope, float intercept){  // XXX HACK 

    //Find the x coordinate on that line at given z
    float x_hit = (z_straw - intercept) / (slope);

    //Find the closest x straw coordinate given x on the line [the vector of straw x is naturally in an ascending order]
    float lastDx=LONG_MAX; //ensure failure of the first comparison
    float thisDx=0; 
    float dca = 0;
    for (int i=0; i<layer_x.size(); i++){
        thisDx=layer_x[i]-x_hit;
        if (abs(lastDx)>abs(thisDx)){
            lastDx=thisDx;
        }
        if (abs(thisDx)>abs(lastDx)){
            dca=lastDx; 
        }
    }

    return dca; 
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
    // [0,0.5]
    rand_num = (( RandomBuffer::instance()->get_uniform_number()+ RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max()));
    float x_0 = beamPositionLength * (rand_num); //uniform vertex
    rand_num = ((RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max()));
    float x_1 = beamPositionLength * (rand_num); //uniform exit point: so fitting a line to these two points
    //float x_slope=(x_1-x_0)/distance[layerN];
    //float x_slope=(x_1-x_0)/(beamStop-beamStart); //beams starts at z=0, beamsStop=25; 
    //XXX Wouldn't the correct detention be: 
    float x_slope=(beamStop-beamStart)/(x_1-x_0); //m=dz/dx?

    float x_intercept = beamStart - x_slope*x_1; 

    float x = x_0;
    float dx = x_slope;
    float previousPosition = 0.0;  // previous position in "z"
   
   //The main loop is for modules [which are misaligned]
    //TODO within modules we will have loops over views, layers, straws 
    for(int i=0; i<layerTotalN; i++){
        
        //Only want to fill for track not all hits
        if (i==0){
        line.x0_gen.push_back(x_0);
        line.z0_gen.push_back(float(beamStart));
        line.x1_gen.push_back(x_1);
        line.z1_gen.push_back(float(beamStop));
        line.x_m.push_back(x_slope);
        line.x_c.push_back(x_intercept);
        }
        
        //distance between consecutive modules 
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
                
        //loops over views, layers, straws 
       // for (int i_view=0; i_view<viewN; i_view++){
            for (int i_layer=0; i_layer<layerN; i_layer++){
                for (int i_straw=0; i_straw<strawN; i_straw++){
                
                } //end of Straws loop
            }//end of Layers lopp
       // }// end of View loop

        // BIG TODO: implement this in a different place allowing to dynamically name and reserve size of vectors/vector # based on geom. constants
        // then used in all for loops with i, n = 0, 1, 2 .....
        for (int i=0; i<mod_0_lyr_0.size(); i++){
            //TODO DCA function goes here
        float hit_distance_1 = DCA(mod_0_lyr_0, distance[i], dx, x_intercept);
        float hit_distance_2 = DCA(mod_0_lyr_1, distance[i], dx, x_intercept);
        float hit_distance_3 = DCA(mod_1_lyr_0, distance[i], dx, x_intercept);
        float hit_distance_4 = DCA(mod_1_lyr_1, distance[i], dx, x_intercept);
        line.dca1.push_back(hit_distance_1);
        line.dca2.push_back(hit_distance_2);
        line.dca3.push_back(hit_distance_3);
        line.dca4.push_back(hit_distance_4);
        }

        // which straw was hit [0,5 Y; 0,10 X] MC rejection if beyond plane geometry
         //Rejection of hits
        //TODO rewrite this ugly part 
        // For now, NO hit rejection 
        int imx=1; //XXX HACK
        if (imx < 0 || imx >= strawN){ //[between (0,10]       
            continue;
        } 
        
        
        //Module number
        line.i_hits.push_back(i); // vector of planes that were actually hit

        float xl=x-sdevX[i]; //misalignment: same for all layers in a moduels 
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
   
    }// end of looping over modules 
    return line; // Return data from simulated track
    
} // end of MC

//Geometry of detector arrangement 
void Tracker::setGeometry(ofstream& debug_geom, bool debugBool){
    float dZ=startingZDistanceModule0; // to add distance in z
    
    //Z position of layers
     //first layer has a specified starting point
    //starting from the second layer there is a pattern 
    // XXX conditions will be modified when multiple views are added 
    for (int layer_n=0; layer_n<layerTotalN; layer_n++){
        distance.push_back(dZ);
        // for layers in the same view (previous position + 0.515)
        if (layer_n%2 == 0){
            dZ = distance[layer_n]+layerSpacing; 
        }

        // for layers in the next module (previous position of first straw + 0.515)
        if (layer_n%2 != 0){
           dZ = distance[layer_n-1]+moduleSpacing; 
        }
    layer.push_back(layer_n);  // layer label array [starting from 0th layer] 
    }// total # layers 

    if (debugBool){
        cout << "Plane Dis. Z: ";
        for (int i = 0; i<distance.size(); i++){
            cout << distance[i] << " ";
        }
        cout << endl;
    }

    // Geometry of detector arrangement in X 
    float dX1 = startingXDistanceModule0; // starting on the x-axis (z, 0)
    float dX2 = layerDisplacement; 

    for (int straw_i=0; straw_i<strawN; straw_i++){
            //TODO 0,1 labels will be more efficiently used in a for loop.... 
            mod_0_lyr_0.push_back(dX1);
            mod_0_lyr_1.push_back(dX2);
            mod_1_lyr_0.push_back(dX1);
            mod_1_lyr_1.push_back(dX2);
            resolutionLayer.push_back(resolution); //resolution
            projectionX.push_back(float(1.0));  // x
            //Constant displacement of straws thereafter 
            dX1=strawSpacing+dX1;
            dX2=strawSpacing+dX2;
    }  // end of looping over straws  
  
    if (debugBool){
        cout << "Mod0 Lyr0 X: ";
        for (int i = 0; i<mod_0_lyr_0.size(); i++){
            cout <<  mod_0_lyr_0[i] << " ";
        }
        cout << endl;
        cout << "Mod0 Lyr1 X: ";
        
        for (int i = 0; i<mod_0_lyr_1.size(); i++){
            cout <<  mod_0_lyr_1[i] << " ";
        }
         cout << endl;
        cout << "Mod1 Lyr0 X: ";
        for (int i = 0; i<mod_1_lyr_0.size(); i++){
            cout <<  mod_1_lyr_0[i] << " ";
        }
         cout << endl;
        cout << "Mod1 Lyr1 X: ";
        
        for (int i = 0; i<mod_1_lyr_1.size(); i++){
            cout <<  mod_1_lyr_1[i] << " ";
        }
        cout << endl;
    }

} // end of geom

// MC misalignment of detectors 
void Tracker::misalign(ofstream& debug_mis, bool debugBool){
    float rand_gaus;
    float sign = 1.0; //for +x or -x direction of misalignment
   //Now misaligning detectors
    for (int i=0; i<moduleN; i++){   
                rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev()); 
                // Before misalignment was smeared XXX
                //sdevX[i] = dispX * rand_gaus;
                sdevX.push_back(dispX * sign);
                sign = -sign;  //next module will be displaced in the opposite direction
    } // end of modules
}//end of misalign


/**
   Write a constraint file to the supplied file-stream.

    @param constraint_file Reference to ofstream to write constraint file to. 
 */
void Tracker::write_constraint_file(ofstream& constraint_file, ofstream& debug_con, bool debugBool) {
    // constraints: fix centre straws in first/last layer
    // Check constraints file is open, then write. 
    if (constraint_file.is_open()) {
       
        //Evaluation of constraints TODO ??? 
        float one = 1.0;
        
        for (int i = 0; i < moduleN; i++){ 
            constraint_file << "Constraint 0.0" << endl;
                int labelt=i+1; // Millepede doesn't like 0 as a label...
                constraint_file << labelt << " " << fixed << setprecision(5) << one<< endl;
                //sdevX[i]=0.0;     // fix centre modules at 0. XXX
           // constraint_file << "Constraint 0.0" << endl;
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