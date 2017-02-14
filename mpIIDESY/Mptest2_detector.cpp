/*

This source file contains definitions of various functions 
in Detector class. 

*/ 

#include "mptest1_port_detector.h"

using namespace std;

// Set up empty pointer for instance of class.
Detector* Detector::s_instance = NULL; 

/** 
	Empty constructor for detector class.
 */
Detector::Detector() {
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
LineData Detector::genlin2() {

	// Set up new container for track data, with hit count set to zero`
	Line_data line;
    line.hit_count = 0;

    // Track parameters for rand-generated line
    float x_0 = layerSize * (rand_gen->Rndm()-0.5); //uniform vertex
    float y_0 = layerSize * (rand_gen->Rndm()-0.5); //uniform vertex 
    float x_1 = layerSize * (rand_gen->Rndm()-0.5); //uniform exit point: so fitting a line to these two points
    float y_1 = layerSize * (rand_gen->Rndm()-0.5); //uniform exit point: 
    float x_slope=(x_1-x_1)/arcLength[layerN];
    float y_slope=(y_1-y_0)/arcLength[layerN];
     

    if (ip != 0){
        cout << "" << endl;
        cout << "Track: " << "x0= " << x_0 << " y0= " << y_0 << "x_slope = " << x_slope << "y_slope = " << y_slope << endl;
   } 
    
    float x = x_0;
    float dx = x_slope;
    float y = y_0;
    float dy = y_slope;
    float sold = 0.0;  // XXX ??? 
    
    for(int i=0; i<layerN; i++){
        float ds = arcLength[i] - sold;
        sold = arcLength[i];

        //position with parameters 1. hit
        float xs=x_0 + arcLength[i] * x_slope;
        float ys=y_0 + arcLength[i] * y_slope;
        
        //true track position
        x=x+dx*ds;
        y=y+dy*ds;

        //multiple scattering
        dx = dx+ rand_gen-> Gaus(0,1) * scatterError;
        dy = dy+ rand_gen -> Gaus(0,1) * scatterError;

        // TODO understand purpose of this part properly 
        float imx=int(x+layerSize*0.5)/layerSize*float(moduleXN);
        if (imx < 0. || imx >= moduleXN) continue;
        float imy=int(y+layerSize*0.5)/layerSize*float(moduleYN);
        if (imy < 0. || imy >= moduleYN) continue;      

        
        int ihit= ((i)*moduleYN+imy)*moduleXN; // XXX i from 0 to 14
        int ioff=((layer[i]-1)*moduleYN+imy)*moduleXN*imx+1;
        //line.i_hits.push_back(ihit);
        line.i_hits.push_back(i); //XXX which one do we need?
        float xl=x-sdevX[ioff];
        float yl=y-sdevY[ioff];
        line.x_hits.push_back(arcLength[i]);
        line.y_hits.push_back((xl-xs)*projection[1][i]+(yl-ys)*projection[2][i]+ rand_gen -> Gaus(0,1)*resolutionLayer[i]);
        line.hit_sigmas.push_back(resolutionLayer[i]);
        line.hit_count++;

        h_2 -> Fill(line.x_hits[i], line.y_hits[i]);

        if (ip != 0){
            cout << "" << endl;
            cout << "Generated Line data: " <<  "Sigma= " << line.hit_sigmas[i] << " X hit= "<< line.x_hits[i] << " Y hit =" << line.y_hits[i] << endl;
        }
    }// end of looping over detector layers
    
    return line; // Return data from simulated track

} // end of genlin2



/**
   Write a constraint file to the supplied file-stream.

    @param constraint_file Reference to ofstream to write constraint file to. 
 */
void Detector::write_constraint_file(ofstream& constraint_file) {

	// Check constraints file is open, then write. [Note - Don't yet understand these]
	if (constraint_file.is_open()) {
		
		cout << "" << endl;
		cout << "Writing constraint file..." << endl;
		cout << "" << endl;
		
		//Evaluation of constraints
        moduleXYN = moduleXN*moduleYN; 

		constraint_file << "Constraint 0.0" << endl;
		for (int i=0; i<PLANE_COUNT; i++) {
			int labelt = 10 + (i + 1) * 2;
			constraint_file << labelt << " " << fixed << setprecision(7) << 1.0 << endl;
		}


		float d_bar = 0.5 * (PLANE_COUNT - 1) * PLANE_X_SEP; 
		float x_bar = PLANE_X_BEGIN + (0.5 * (PLANE_COUNT - 1) * PLANE_X_SEP);
		constraint_file << "Constraint 0.0" << endl;
		for (int i=0; i<PLANE_COUNT; i++) {
			
			int labelt = 10 + (i + 1) * 2;
			
			float x = PLANE_X_BEGIN + (i * PLANE_X_SEP);
			float ww = (x - x_bar) / d_bar;

			constraint_file << labelt << " " << fixed << setprecision(7) << ww << endl;
		} //plabe count	
	} // constraing file open 
} // end of wrting cons file 









