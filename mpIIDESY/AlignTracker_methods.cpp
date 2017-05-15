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
   Simple function to calculate the Distance of Closest Approach
   between a point (x0, y0) and a line (defined by two points: (x1, y1), (x2, y2) ).
  
   @return dca
*/
float Tracker::DCA(float x0,float y0,float x1,float y1,float x2,float y2){

    float dca = abs((y2-y1)*x0-(x2-x1)*y0+x2*y1-y2*x1)/sqrt((y2-y1)*(y2-y1)+(x2-x1)*(x2-x1));

    return dca;
}

float Tracker::DCA_simple(float x_0,float x_straw){

    //accounting for negative coordinates 
    float dca = sqrt((x_0-x_straw)*(x_0-x_straw));

    return dca;
}



/** Uses DCA function to find the shortest dca between straws (i.e. which straw was hit in that layer)

    @Inputs (see DCA method) + vector of straws' x coordinates in a layer, and a return type: "dca_hit"  or "x_line"
  
    @return hit_distance - dca of the straw that was hit, or x coordinate of the line on the same z as a straw
*/
float Tracker::DCAHit(std::vector<float> layer_x, float z_straw, float x1,float z1,float x2, float z2, float c, float m, bool debugBool){ 

    bool StrongDebugBool=false; //quick hack XXX for even more debug output

    //Find the two closest x straw coordinates given x on the line - with same z [the vector of straw x is naturally in an ascending order]
   //TODO ensure this also true for >= and <= for hits going through the actual wire!! 
   // if hit bigger that all straws all smaller  

    float x_line = (z_straw-c)/m;
    
    float lower, upper, hit_distance;

    vector<float>::iterator it;
    it = lower_bound(layer_x.begin(), layer_x.end(), x_line);
    if (it == layer_x.begin()) upper = *it; // no smaller value  than val in vector
    else if (it == layer_x.end()) lower = *(it-1); // no bigger value than val in vector
    else {
        lower = *(it-1);
        upper = *it;
    }
   
    float hit_distance_low = Tracker::DCA(lower, z_straw, x1, z1, x2, z2);
    float hit_distance_up = Tracker::DCA(upper, z_straw, x1, z1, x2, z2); 

    if (hit_distance_up>hit_distance_low){hit_distance = hit_distance_low;}
    if (hit_distance_up<hit_distance_low){hit_distance = hit_distance_up;}
    if (hit_distance_up==hit_distance_low){hit_distance = hit_distance_low; cout << "Ambiguity which straw registered hit";}  //XXX

    if (debugBool && StrongDebugBool){
        cout << "Track (line) point in-line with the layer at ( " << x_line << " , "  << z_straw << " ); belongs to the line with  generated points: ( " << x1 <<" , " << beamStart << ");" <<endl;
        cout << "( " << x2 << " , " << beamStop << " ); slope= " << m << "; intercept= " << c << endl;
        cout << "Two straws closest to that point are " << lower << ", and " << upper <<  "; with DCAs "<< hit_distance_low <<  " and " << hit_distance_up << ", respectively." << endl;
        cout << "Selected DCA as the correct hit distance is " << hit_distance << endl;
    }
    
        return hit_distance; // as dca of the closest straw
}


DCAData Tracker::DCAHit_simple(std::vector<float> layer_x, float z_straw, float x_0, bool debugBool){ 

    DCAData dca_data;

    //Only storing dca and id for the current layer
    //dca_data.strawID.clear(); 
    //dca_data.dca.clear();
    //dca_data.LRSign.clear();

    bool StrongDebugBool=false; //quick hack XXX for even more debug output

    //Find the two closest x straw coordinates given x on the line - with same z [the vector of straw x is naturally in an ascending order]
   //TODO ensure this also true for >= and <= for hits going through the actual wire!! 
   // if hit bigger that all straws all smaller  

    float x_line = x_0;

    if (x_0>beamPositionLength || x_0 < 0){
        cout << "Error: generated point out of bounds";
        exit(0);
    }
    
    float lower, upper, hit_distance;
    int vector_lastID = strawN-1; // the ID of the last straw in the vector
    float comparator; // to find the ID
    float LR; //L or R hit

    if (debugBool && StrongDebugBool && 1==0){
        //cout << "Track (line) point in-line with the layer at ( " << x_line << " , "  << z_straw << " ); belongs to the line with  generated points: ( " << x_0 <<" , " << beamStart << ");" <<endl;
        //cout << "( " << x_0 << " , " << beamStop << " );" << endl;
        cout << "Hit at " << x_0 << endl;
    }

    //Iterator to find two straws between the hit:
    vector<float>::iterator it;
    it = lower_bound(layer_x.begin(), layer_x.end(), x_line);
    // no smaller value  than val in vector
    if (it == layer_x.begin()){
        upper = *it; 
        //if (debugBool && StrongDebugBool){ cout << " upper = " << upper << endl;}
         // hit distance is the dca from the L of the straw 0
        hit_distance = Tracker::DCA_simple(upper, x_0);
        comparator = upper;
        LR=-1.0;
        if (debugBool && StrongDebugBool){cout << "Hit at " << x_0 << " The first straw is closest at " << upper <<  "; with DCA "<< hit_distance << endl;}
    }   
    // no bigger value than val in vector
    else if (it == layer_x.end()){
        lower = *(it-1); 
        //if (debugBool && StrongDebugBool){ cout << "lower = " << lower << endl;}
        // hit distance is the dca from the R of the last straw
        hit_distance = Tracker::DCA_simple(lower, x_0); // hit distance is the dca from the L of the straw 0
        comparator = lower;
        LR=1.0;
        if (debugBool && StrongDebugBool){cout << "Hit at " << x_0 << " The last straw is closest at " << lower <<  "; with DCA "<< hit_distance << endl;}
    }
    else {
        lower = *(it-1);
        upper = *it;
        //if (debugBool && StrongDebugBool){ cout << "lower = " << lower << " upper = " << upper << endl;}

            float hit_distance_low = Tracker::DCA_simple(lower, x_0);
            float hit_distance_up = Tracker::DCA_simple(upper, x_0); 

            if (hit_distance_up>hit_distance_low){
                hit_distance = hit_distance_low;
                comparator = lower;
                LR=+1.0;
            }
            if (hit_distance_up<hit_distance_low){
                hit_distance = hit_distance_up;
                comparator = upper;
                LR=-1.0;
            }
            //XXX resolve this [won't be an issue once physical geometry and hit rejection is added]
            if (hit_distance_up==hit_distance_low){
                hit_distance = hit_distance_low; cout << "Ambiguity which straw registered hit" << endl;
                exit(0);
            } 
            if (debugBool && StrongDebugBool && 1==0){
                cout << "Two straws closest to that point are " << lower << ", and " << upper <<  "; with DCAs "<< hit_distance_low <<  " and " << hit_distance_up << ", respectively." << endl; 
            }
    } // end of iterator to find straws between hits
   
    //Iterator to find straw ID [implemented separately, as now we look up once only]:
    int index;
    it = find(layer_x.begin(), layer_x.end(), comparator);
    if (it == layer_x.end()){   
        cout << "Error: straw x not found in vector!" << endl;
        exit(0);
      //TODO error happens here: wrong side of the straw + crazy ID number. 
        // XXX Good test for the algorithm - investigate here :) 
    } else {
        index = std::distance(layer_x.begin(), it);
        if (index > vector_lastID || index < 0){
            cout << "Error: ID is out of range" << endl;
            exit(0);
        }
        dca_data.strawID = index;
    }

    dca_data.dca = hit_distance;
    dca_data.LRSign = LR;

    if (debugBool && StrongDebugBool && 1==0){
       
        cout << "Selected DCA as the correct hit distance is " << hit_distance << ". Straw ID: " << dca_data.strawID;
        if (LR  > 0){
            cout << ". The straw was hit from the right" << endl;
        }
        if (LR < 0) {
            cout << ". The straw was hit from the left" << endl;
        }
    }//debug
        return dca_data; // as dca and id of the closest straw
}



float Tracker::GetIdealLineX(int x_det_ID, float x_det_dca, float LRSign, vector<float> layer_x){
    
    float x_IdealStraw=layer_x[x_det_ID];

    float x_Fitline = x_IdealStraw + (LRSign * x_det_dca); //add dca for R 

    return x_Fitline;
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
    //rand_num = ((RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max()));
    //float x_1 = beamPositionLength * (rand_num); //uniform exit point: so fitting a line to these two points
    
    //float x_slope=(x_1-x_0)/(beamStop-beamStart); //beams starts at z=0, beamsStop=25; 
    //XXX Wouldn't the correct detention be: 
    //float x_slope=(beamStop-beamStart)/(x_1-x_0); //m=dz/dx? Slope for a line with coordinates (x_0, beamStart), (x_1, beamStop)
    //float x_intercept = beamStart - x_slope*x_0; // for line z=mx+c -> c= z - mx

    float x_track = x_0; //true track position x=(z_straw-c)/m;
    //float dx = x_slope;
    float previousPosition = 0.0;  // previous position in "z"
    float x_det; // this will the dca from the hit_layer

    //debug_calc << x_0 << " " << x_1 << endl;
    
    //Sanity Plots: Tracks
    line.x0_gen.push_back(x_0);
    line.z0_gen.push_back(float(beamStart));
    //line.x1_gen.push_back(x_1);
    line.z1_gen.push_back(float(beamStop));
    //line.x_m.push_back(x_slope);
    //line.x_c.push_back(x_intercept);

   
    //The main loop is for modules [they produce label of Global Parameters]: 
    // Then looping over layers [we know there will be at most 1 hit per layer] and views
    // within layers we will have loops over straws in that layer [to find out which straw was hit and dca]
    int z_counter=0; // distance in z is incremented from 0 to TotalLayerN 
    for(int i_module=0; i_module<moduleN; i_module++){  
        
    //This was left over from MP2 Toy model: we are implementing dca instead...some useful methods still might be there... XXX
    #if 0  
        //distance between consecutive modules 
        float dPostion = distance[i] - previousPosition; // in Z
        previousPosition = distance[i];
        //positions on the detector plane  
        float x_det=x_0 + distance[i] * x_slope;
        line.x_det.push_back(x_det);
        //true track position
        x_track=x_track+dx*dPostion;  
        line.x_true.push_back(x_track);
        //multiple scattering
        rand_gaus= RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
        dx = dx+ rand_gaus * scatterError;
    #endif

        //loops over views, layers, straws 
       // for (int i_view=0; i_view<viewN; i_view++){
        for (int i_layer=0; i_layer<layerN; i_layer++){ //per layer
            //for (int i_straw=0; i_straw<strawN; i_straw++){

                //x_track = (distance[z_counter]-x_intercept)/x_slope;   // true track position [from line] WRONGG.... XXX need dca here
                x_track = x_0; // the x of the parallel line  
                //the dca between x_0 and the x_straw [misaligned]
                
                DCAData x_det= Tracker::DCAHit_simple(mod_lyr_strawMisPosition[i_module][i_layer], distance[z_counter], x_0, debugBool); // position on the detector [from dca]

                cout << x_det.dca << endl;
                float x_det_dca =  x_det.dca;  //dca of the hit straw
                cout << x_det_dca << endl;
 
                int x_det_ID =  x_det.strawID; // if of the hit straw [to identify the correct straw in the Fit function]

                float x_det_LRSign = x_det.LRSign;

                float x_ideal = Tracker::GetIdealLineX(x_det_ID, x_det_dca, x_det_LRSign, mod_lyr_strawIdealPosition[i_module][i_layer]);

                x_fitted.push_back(x_ideal); // vector to store x coordinates of the track as seen from the ideal detector

                //Rejection of hits due to geometry (i.e. missed hits) TODO for later
                //No signal in a straw = no signal in the layer
                 // if (x_det > strawRadius){ //[between (0,0.25]       
                //     //continue;
                // } 
                
                //Module number [for labelling]
                line.i_hits.push_back(i_module); // vector of planes that were actually hit [after passing rejection test]

                //float x_mis=x_track-sdevX[i_module]; //x position due to misalignment: same for all layers in a modules [true position - misalignment]
                
                
                //Z-coordinate of hits [it was always known from geometry - no z-misalignment for now...]
                line.z_hits.push_back(distance[z_counter]);
                // the residual looks to be deltaX + deltaY rather than the magnitude of the distance... worth noting?
                // projection Y is always 0 for non-stereo modules?? what is the motivation? 
                rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
                //float hit = (x_mis-x_det)*projectionX[i_module] + rand_gaus *resolution; //RESIDUAL between hit due to misalignment [] and detected position * smeared resolution 
                float dca_misStraw = x_det_dca*projectionX[i_module] + rand_gaus *resolution; //RESIDUAL between hit due to misalignment [] and detected position * smeared resolution 
                line.x_hits.push_back(dca_misStraw);
                line.hit_sigmas.push_back(resolution);
                line.hit_count++;

                //Sanity Plots 
                //Hits
                line.x_det.push_back(x_det_dca);
                //line.x_mis.push_back(x_mis);
                line.x_true.push_back(x_track);

           // } //end of Straws loop
        z_counter++;
        }//end of Layers loop
   // }// end of View loop
   
    }// end of looping over modules 
    return line; // Return data from simulated track
    
} // end of MC

//Geometry of detector arrangement (Ideal Geometry)
void Tracker::setGeometry(ofstream& debug_geom, bool debugBool){
    float dZ=startingZDistanceStraw0; // to add distance in z
    
     //Z position of layers
     //first layer has a specified starting point
    //starting from the second layer there is a pattern 
    // XXX conditions will be modified when multiple views are added 
    //We have layerTotalN elements in z to be matched with potion in z
    for (int layer_n=0; layer_n<layerTotalN; layer_n++){
        distance.push_back(dZ);
        // for layers in the same view (previous position + 0.515)
        if (layer_n%2 == 0){
            dZ = distance[layer_n]+layerSpacing; 
        }

        // for layers in the next module (previous position of first straw + 18.735)
        // This will be modified with addition of views XXX
        if (layer_n%2 != 0){
           dZ = distance[layer_n-1]+moduleSpacing; 
        }
    layer.push_back(layer_n);  // layer label array [starting from 0th layer]
    resolutionLayer.push_back(resolution); //resolution in each layer
    projectionX.push_back(float(1.0));  // x projection of hits in each layer
    }// total # layers 


    // Geometry of detector arrangement in X 
    float dX = startingXDistanceStraw0; // starting on the x-axis (z, 0)

    // XXX conditions will be modified when multiple views are added 
    for (int i_module=0; i_module<moduleN; i_module++){
        mod_lyr_strawIdealPosition.push_back(vector<vector<float> >()); //initialize the first index with a 2D vector
         //for (int i_view=0; i_view<viewN; i_view++){
                for (int i_layer=0; i_layer<layerN; i_layer++){ //per module
                    mod_lyr_strawIdealPosition[i_module].push_back(vector<float> ()); //initialize the first index with a 2D vector
                    for (int i_straw=0; i_straw<strawN; i_straw++){
                        mod_lyr_strawIdealPosition[i_module][i_layer].push_back(dX); 
                        dX=strawSpacing+dX; //while we are in the same layer: increment straw spacing in x
                    } //end of Straws loop
                    dX=layerDisplacement; //set displacement in x for the next layer in the view
                }//end of Layers loop
                dX=startingXDistanceStraw0;
           // }// end of View loop
        }//Modules  
  
    if (debugBool){
        int Zcounter=0;
        for (int i_module=0; i_module<moduleN; i_module++){
          //for (int i_view=0; i_view<viewN; i_view++){
                for (int i_layer=0; i_layer<layerN; i_layer++){ //per module
                    cout << "IDEAL Mod " << i_module  << " Layer " << i_layer << " X : ";
                    for (int i_straw=0; i_straw<strawN; i_straw++){
                    cout << mod_lyr_strawIdealPosition[i_module][i_layer][i_straw]<<" ";
                    } //end of Straws loop
                    cout << " | Z= " << distance[Zcounter] << " [cm]" << endl;  // XXX fix this 
                    Zcounter++;
                } // end of Layers
                
           // }// end of View loop
        }//Modules 
    }


} // end of geom

// MC misalignment of detectors 
//TODO 
void Tracker::misalign(ofstream& debug_mis, bool debugBool){
    //float rand_gaus;
    float sign = 1.0; //for +x or -x direction of misalignment
    //Now misaligning detectors in x

    float dX = startingXDistanceStraw0+(dispX * sign); // starting on the x-axis (z, 0+dispX)

    for (int i_module=0; i_module<moduleN; i_module++){
        //rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev()); 
        // Before misalignment was also smeared XXX
        //sdevX[i] = dispX * rand_gaus;
        sdevX.push_back(dispX * sign);  // vector of misalignment 
        mod_lyr_strawMisPosition.push_back(vector<vector<float> >()); //initialize the first index with a 2D vector
         //for (int i_view=0; i_view<viewN; i_view++){
                for (int i_layer=0; i_layer<layerN; i_layer++){ //per module
                    mod_lyr_strawMisPosition[i_module].push_back(vector<float> ()); //initialize the first index with a 2D vector
                    for (int i_straw=0; i_straw<strawN; i_straw++){
                        mod_lyr_strawMisPosition[i_module][i_layer].push_back(dX); 
                        dX=strawSpacing+dX; //while we are in the same layer: increment straw spacing in x
                    } //end of Straws loop
                    dX=layerDisplacement+(dispX * sign); //set displacement in x for the next layer in the view
                }//end of Layers loop
                sign = -sign;  //next module will be displaced in the opposite direction
                dX=startingXDistanceStraw0+(dispX * sign);
           // }// end of View loop
        
        }//Modules 

    if (debugBool){
        int Zcounter=0;
        for (int i_module=0; i_module<moduleN; i_module++){
          //for (int i_view=0; i_view<viewN; i_view++){
                for (int i_layer=0; i_layer<layerN; i_layer++){ //per module
                    cout << "MIS Mod " << i_module  << " Layer " << i_layer << " X : ";
                    for (int i_straw=0; i_straw<strawN; i_straw++){
                    cout << mod_lyr_strawMisPosition[i_module][i_layer][i_straw]<<" ";
                    } //end of Straws loop
                    cout << " | Z= " << distance[Zcounter] << " [cm]" << endl;  // XXX fix this 
                    Zcounter++;
                } // end of Layers
                
           // }// end of View loop
        }//Modules 
    }

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