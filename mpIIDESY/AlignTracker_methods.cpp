/*

This source file contains definitions of various functions in the method class. 

*/  

#include "AlignTracker_methods.h"

using namespace std;

// Set up empty pointer for instance of class.
Tracker* Tracker::s_instance = NULL; 

/** 
	Empty constructor for detector class.
 */
Tracker::Tracker() {
  //XXX non-static variables definition here 

    resolution=0.002;  // 20um = 0.002 cm setting this in the constructor //TODO decide if this is a constant 

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
  
	// Set up new container for track data, with hit count set to zero`
	LineData line;
    line.hit_count = 0;

    float rand_num;
    float rand_gaus;

   if (debugBool){
           debug_calc << "Track # (C)        " << line.hit_count << endl;
           debug_calc << "–––––––––––––––––––––––––––––––––––––––––––––––" <<  endl;
           debug_calc << endl; 
       }

    // Track parameters for rand-generated line MC 
    rand_num = ( RandomBuffer::instance()->get_uniform_number()+ RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max());
    if (debugBool){debug_mc << "rand_num= " <<rand_num << " uniform_ran_max= "<<RandomBuffer::instance()->get_uniform_ran_max()<<" two= "<<twoR<< endl;
    debug_mc << "Rand= " << rand_num;}
    float x_0 = layerSize * (rand_num-0.5); //uniform vertex
    if (debugBool){debug_mc << " x0= " << x_0 << endl;}
    rand_num = (RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max());
    if (debugBool){ debug_mc << " Rand= " << rand_num;}
    float y_0 = layerSize * (rand_num-0.5); //uniform vertex 
    if (debugBool){ debug_mc << " y0= " << y_0 << endl;}
    rand_num = (RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max());
    if (debugBool){ debug_mc << "Rand= " << rand_num; }
    float x_1 = layerSize * (rand_num-0.5); //uniform exit point: so fitting a line to these two points
    if (debugBool){debug_mc << " x1= " << x_1 << endl; }
    rand_num = (RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max());
    if (debugBool){debug_mc << " Rand= " << rand_num; }
    float y_1 = layerSize * (rand_num-0.5); //uniform exit point: 
    if (debugBool){ debug_mc << " y1= " << y_1 << endl; }
    float x_slope=(x_1-x_0)/distance[layerN-1];
    if (debugBool){ debug_mc << " x_slope= " << x_slope << " distance[layerN-1]= " << distance[layerN-1] << endl; }
    float y_slope=(y_1-y_0)/distance[layerN-1];
    if (debugBool){ debug_mc << " y_slope= " << y_slope <<  " distance[layerN-1]= " << distance[layerN-1] << endl;}
    

    float x = x_0;
    float dx = x_slope;
    float y = y_0;
    float dy = y_slope;
    float s_old = 0.0;  // previous position in "z"

    if (debugBool){debug_mc << "x= " << x << " dx= " << dx << " y= " << y << " dy= " << dy << " s_old= " << s_old << endl; 
    debug_mc << endl; }
    if (debugBool){
              debug_calc << "x_0= "<< x_0<< " y_0= "<< y_0 << " x_1= "<< x_1<< " y_1= "<< y_1 << endl; 
               debug_calc << "x_slope= "<< x_slope<< " y_slope= "<< y_slope << endl;
               debug_calc << endl;
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
        if (debugBool){debug_mc << "Rand= " << rand_gaus << " dx= " << dx << " scatterError= " << scatterError << endl;}
        rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
        dy = dy+ rand_gaus * scatterError;
        if (debugBool){debug_mc << "Rand= " << rand_gaus << " dy= " << dy << " scatterError= " << scatterError << endl;}

        if (debugBool){
        debug_mc << "ds= " << ds << " distance[i]= " << " s_old= " << s_old <<endl; 
        debug_mc << "xs= " << xs << " x_0= " << x_0 << " x_slope= " << x_slope << endl;
        debug_mc << "ys= " << ys << " y_0= " << y_0 << " y_slope= " << y_slope << endl;
        debug_mc << "x= " << x << " dx= " << dx << " ds= " << ds << endl;
        debug_mc << "y= " << y << " dy= " << dy << " ds= " << ds << endl;\
        }

        
        // which pixel was hit [0,5 Y; 0,10 X] MC rejection if beyond plane geometry
        int imx=int((x+layerSize*0.5)/layerSize*float(pixelXN));
        if (debugBool){debug_off << (x+layerSize*0.5)/layerSize*float(pixelXN) << endl;}
        if (imx < 0 || imx >= pixelXN){ //[between (0,10]       
            continue;
        } 
        int imy=int((y+layerSize*0.5)/layerSize*float(pixelYN));
        if (debugBool){ debug_off << (y+layerSize*0.5)/layerSize*float(pixelYN) << endl;}
         if (imy < 0 || imy >= pixelYN){  //[between (0,5] 
            continue;
        }
       
        int ihit= ((i)*pixelYN+imy)*pixelXN+imx; // i from 0 to 13 (incl.)
        line.i_hits.push_back(ihit); // vector of planes that were actually hit

        float xl=x-sdevX[layer[i]-1][imy][imx]; //misalign. 
        float yl=y-sdevY[layer[i]-1][imy][imx];
        line.x_mis.push_back(xl);
        line.y_mis.push_back(yl);

        if (debugBool){debug_mc << "xl= " << xl << " yl= " << yl << " sdevX= " << sdevX[layer[i]-1][imy][imx] << " sdevY= " << sdevY[layer[i]-1][imy][imx] << endl;}
       
        // we seem to now redefine the coordinates so that x is now the distance and y is a measure of the residual
        line.x_hits.push_back(distance[i]);
        // the residual looks to be deltaX + deltaY rather than the magnitude of the distance... worth noting?
        // projection Y is always 0 for non-stero modules?? what is the motivation? 
        rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
        float yhit = (xl-xs)*projectionX[i]+(yl-ys)*projectionY[i]+ rand_gaus *resolution; 
        if (debugBool){debug_mc << "yhit = " << yhit << "rand_num= " << rand_gaus << " projectionY[i]= " << projectionY[i] << " projectionX[i]= " << projectionX[i] << " resolution= " << resolution << endl;}
        line.y_hits.push_back(yhit);
        line.hit_sigmas.push_back(resolution);
        line.hit_count++;

       if (debugBool){
               debug_calc << "xs= " << xs << "  ys= " << ys << "  x= " << x << "  y= " << y << endl;
               debug_calc << "imx= " << imx << "  imy= " << imy << endl;
               debug_calc << "ihit= " << ihit << "  xl= " << xl << "  yl= " << yl << "  xhit= " << distance[i]  << "  yhit= " << yhit << endl;
               debug_calc << "sdevX[layer[i]-1][imy][imx]= " << sdevX[layer[i]-1][imy][imx] << " sdevY[layer[i]-1][imy][imx]= " << sdevY[layer[i]-1][imy][imx] << endl; 
               debug_calc << "projectionX[i]= " << projectionX[i] << " projectionY[i]= " << projectionY[i] << endl; 
               debug_calc << "nhits= " << line.hit_count << endl; 
               debug_calc << endl; 
            }      
   
    }// end of looping over detector layers
    if (debugBool){
    debug_mc << "------------------------------------------------------------------------------" << endl; 
    debug_mc << endl; 
    }
    return line; // Return data from simulated track
    
} // end of MC

//Geometry of detecor arrangement 
void Tracker::setGeometry(ofstream& debug_geom, bool debugBool){
	float s=startingDistancePlane1;
    int i_counter = 0;
    float sign = 1.0;

    // Geometry of detecor arrangement 
    for (int layer_i=1; layer_i<=10; layer_i++){
        
        layer.push_back(layer_i);  // layer [starting from 1st layer]
        distance.push_back(s);  //distance between planes  [14]
        resolutionLayer.push_back(resolution); //resolution
        projectionX.push_back(1.0);  // x
        projectionY.push_back(0.0);  // y
        if (debugBool){
               debug_geom << "layer_i= " << layer_i << " layer[i_counter]= " << layer[i_counter]  << endl;
               debug_geom << "i_counter= " << i_counter << " distance[i_counter] " << distance[i_counter] << endl;
               debug_geom << "projectionX= " << projectionX[i_counter] << " projectionY= " << projectionY[i_counter] << endl;
               debug_geom << endl; 
        } 
        i_counter++; 
        //taking care of stereo planes [have no pixels] 1, 4, 7, 10
        if (((layer_i) % 3) == 1){
            layer.push_back(layer_i);  // layer
            distance.push_back(s+offset);  //distance between planes  [14]
            resolutionLayer.push_back(resolution); //resolution
            projectionX.push_back(std::sqrt(1.0-std::pow(stereoTheta,2)));  // x
            projectionY.push_back(stereoTheta*sign);  // y
            if (debugBool){
               debug_geom << "S_layer_i= " << layer_i << " layer[i_counter]= " << layer[i_counter]  << endl;
               debug_geom << "i_counter= " << i_counter << " distance[i_counter] " << distance[i_counter] << endl;
               debug_geom << "projectionX= " << projectionX[i_counter] << " projectionY= " << projectionY[i_counter] << endl;
               debug_geom << endl; 
            }  
            sign=-sign;
            i_counter++;    
        }

        s=s+planeDistance;  // incrimenting distance between detecors

    }  // end of looping over layers


} // end of geom

// MC misalignment of detecors 
void Tracker::misalign(ofstream& debug_mis, bool debugBool){

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
                rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
                sdevY[i][k][l] = dispY * rand_gaus;
                counterMis++;

                if (debugBool){
                    debug_mis << "i= " << i << " k= " << k << " l= " << l << endl;
                    debug_mis << "sdevX[i][k][l]= " << sdevX[i][k][l] << " sdevY[i][k][l]= " << sdevY[i][k][l] << endl;  
                    debug_mis << endl; 
                }  
            
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
void Tracker::write_constraint_file(ofstream& constraint_file) {

    // constraints: fix center pixels in first/last layer
	// Check constraints file is open, then write. 
	if (constraint_file.is_open()) {
		

		//Evaluation of constraints
    	int ncx = (pixelXN+1)/2; 
    	int lunt = 9;
    	float one = 1.0;
         

		for (int i = 1; i <= detectorN; i=i+(detectorN-1)){ 
            constraint_file << "Constraint 0.0" << endl;
            for (int k=0; k<=pixelYN-1; k++){
                int labelt=(i*pixelYN+k)*pixelXN+ncx-1;
                constraint_file << labelt << " " << fixed << setprecision(5) << one<< endl;
                sdevX[detectorN][pixelYN][pixelXN]=0.0;      // fix center pixels at 0.
            } // end of y loop
            constraint_file << "Constraint 0.0" << endl;
            for(int k=0; k<=pixelYN-1; k++){
                int labelt=(i*pixelYN+k)*pixelXN+ncx+1000-1;
                constraint_file << labelt << " " << fixed << setprecision(5) << one<< endl;
                sdevY[detectorN][pixelYN][pixelXN]=0.0;

            } // end of x loop
        } // end of detecors loop 

	} // constraing file open 
} // end of wrting cons file 

/**
   Set filename to read uniform random numbers from.

   @param uniform_filename String for filename of file of uniform random numbers.
 */
void Tracker::set_uniform_file(string uniform_filename) {
	RandomBuffer::instance()->open_uniform_file(uniform_filename);
}


/**
   Set filename to read gaussian random numbers from.

   @param uniform_filename String for filename of file of gaussian random numbers.
 */
void Tracker::set_gaussian_file(string gaussian_filename) {
	RandomBuffer::instance()->open_gaussian_file(gaussian_filename);
}







