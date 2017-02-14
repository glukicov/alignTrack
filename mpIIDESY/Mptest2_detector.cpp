/*

This source file contains definitions of various functions 
in Detector class. 

*/ 

#include "Mptest2_detector.h"

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
	LineData line;
    line.hit_count = 0;

    // Track parameters for rand-generated line
    float x_0 = layerSize * (RandomBuffer::instance()->get_uniform_number()-0.5); //uniform vertex
    float y_0 = layerSize * (RandomBuffer::instance()->get_uniform_number()-0.5); //uniform vertex 
    float x_1 = layerSize * (RandomBuffer::instance()->get_uniform_number()-0.5); //uniform exit point: so fitting a line to these two points
    float y_1 = layerSize * (RandomBuffer::instance()->get_uniform_number()-0.5); //uniform exit point: 
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
        dx = dx+ RandomBuffer::instance()->get_gaussian_number() * scatterError;
        dy = dy+ RandomBuffer::instance()->get_gaussian_number() * scatterError;

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
        line.y_hits.push_back((xl-xs)*projection[1][i]+(yl-ys)*projection[2][i]+ RandomBuffer::instance()->get_gaussian_number()*resolutionLayer[i]);
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

//Geometry of detecor arrangement 
void setGeometry(){
	float s=arcLength_Plane1;
    int i_counter = 0;
    float sign = 1.0;

    // Geometry of detecor arrangement 
    for (int layer_i=1; layer_i<=10; layer_i++){
        i_counter++;
        layer[i_counter] = layer_i;  // layer
        arcLength[i_counter] = s;  //arclength
        resolutionLayer[i_counter] = resolution; //resolution
        projection[1][i_counter]=1.0;  // x
        projection[2][i_counter]=0.0;  // y
        //taking care of stereo modules 
        if ((layer_i % 3) == 1){
            i_counter++;
            layer[i_counter] = layer_i;  // layer
            arcLength[i_counter] = s+offset;  //arclength
            resolutionLayer[i_counter] = resolution; //resolution
            projection[1][i_counter]=sqrt(1.0-pow(stereoTheta,2));  // x
            projection[2][i_counter]=stereoTheta*sign;  // y
            sign=-sign;
        }
        s=s+planeDistance;  // incrimenting distance between detecors 
    }  // end of looping over layers


} // end of geom

// MC misalignment of detecors 
void misalign(){
	//Now misaligning detecors
    float dispX = 0.01; // module displacement in Y .05 mm * N(0,1)
    float dispY = 0.01; // module displacement in Y .05 mm * N(0,1)

    //so we are only displacing 9/10 detectors? XXX
    for (int i=0; i<detectorN-1; i++){   /// XXX
        for(int k=0; k<=moduleYN-1; k++){
            for(int l=1; l<=moduleXN; l++){
                sdevX[(i*moduleYN+k)*moduleXN+l] = dispX * RandomBuffer::instance()->get_gaussian_number();  //XXX where is that used in imodel=0?? 
                sdevY[(i*moduleYN+k)*moduleXN+l] = dispY * RandomBuffer::instance()->get_gaussian_number();         
            } // // end of number of modules in x
        } // end of number of modules in y 
    } // end of layers

}//end of misalign


/**
   Write a constraint file to the supplied file-stream.

    @param constraint_file Reference to ofstream to write constraint file to. 
 */
void Detector::write_constraint_file(ofstream& constraint_file) {
//TODO fix constraint file 

	// Check constraints file is open, then write. [Note - Don't yet understand these]
	if (constraint_file.is_open()) {
		
		cout << "" << endl;
		cout << "Writing constraint file..." << endl;
		cout << "" << endl;

		//Evaluation of constraints
    	int ncx = (moduleXN+1)/2; 
    	int lunt = 9;
    	float one = 1.0;
         

		for (int i = 1; i <= detectorN; i=i+(detectorN-1)){  //XXX coorect implimentation of DO i=1,nlyr,nlyr-1
            constraint_file << "Constraint 0.0" << endl;
            for (int k=0; k<=moduleYN-1; k++){
                int labelt=(i*moduleYN+k)*moduleXN+ncx-1;
                constraint_file << labelt << " " << fixed << setprecision(7) << one<< endl;
                sdevX[((i-1)*moduleYN+k)*moduleXN+ncx]=0.0;      // fix center modules at 0.
            } // end of y loop
            constraint_file << "Constraint 0.0" << endl;
            for(int k=0; k<=moduleYN-1; k++){
                int labelt=(i*moduleYN+k)*moduleXN+ncx+1000-1;
                constraint_file << labelt << " " << fixed << setprecision(7) << one<< endl;
                sdevY[((i-1)*moduleYN+k)*moduleXN+ncx]=0.0; // fix center modules at 0.
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







