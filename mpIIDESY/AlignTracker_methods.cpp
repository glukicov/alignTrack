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
  // non-static variables definition here 
    resolution=0.01;  // 100um = 0.01 cm setting this in the constructor //TODO decide if this is a constant 
    dispX = 0.02;  // manual displacement by 0.15 cm

    rejectedHitsDCA = 0;
    multipleHitsLayer = 0;
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
float Tracker::DCA(float xPoint,float yPoint,float x1,float y1,float x2,float y2){

    float dca = abs((y2-y1)*xPoint-(x2-x1)*yPoint+x2*y1-y2*x1)/sqrt((y2-y1)*(y2-y1)+(x2-x1)*(x2-x1));

    return dca;
}

//DCA in 1D for parallel tracks
float Tracker::DCA_simple(float xTrack,float xStraw){

    //accounting for negative coordinates 
    float dca = sqrt((xTrack-xStraw)*(xTrack-xStraw));

    return dca;
}



/** Uses DCA function to find the shortest dca between straws (i.e. which straw was hit in that layer)

    @Inputs (see DCA method) + vector of straws' x coordinates in a layer, and a return type: "dca_hit"  or "x_line"
  
    @return hit_distance - dca of the straw that was hit, or x coordinate of the line on the same z as a straw
*/
float Tracker::DCAHit(std::vector<float> xLayer, float zStraw, float x1,float z1,float x2, float z2, float c, float m, bool debugBool){ 

    // Logic out-of-date XXXX

   //  bool StrongDebugBool=false; //quick hack XXX for even more debug output

   //  //Find the two closest x straw coordinates given x on the line - with same z [the vector of straw x is naturally in an ascending order]
   // //TODO ensure this also true for >= and <= for hits going through the actual wire!! 
   // // if hit bigger that all straws all smaller  

   //  float x_line = (z_straw-c)/m;
    
   //  float lower, upper, hit_distance;

   //  vector<float>::iterator it;
   //  it = lower_bound(layer_x.begin(), layer_x.end(), x_line);
   //  if (it == layer_x.begin()) upper = *it; // no smaller value  than val in vector
   //  else if (it == layer_x.end()) lower = *(it-1); // no bigger value than val in vector
   //  else {
   //      lower = *(it-1);
   //      upper = *it;
   //  }
   
   //  float hit_distance_low = Tracker::DCA(lower, z_straw, x1, z1, x2, z2);
   //  float hit_distance_up = Tracker::DCA(upper, z_straw, x1, z1, x2, z2); 

   //  if (hit_distance_up>hit_distance_low){hit_distance = hit_distance_low;}
   //  if (hit_distance_up<hit_distance_low){hit_distance = hit_distance_up;}
   //  if (hit_distance_up==hit_distance_low){hit_distance = hit_distance_low; cout << "Ambiguity which straw registered hit";}  //XXX

   //  if (debugBool && StrongDebugBool){
   //      cout << "Track (line) point in-line with the layer at ( " << x_line << " , "  << z_straw << " ); belongs to the line with  generated points: ( " << x1 <<" , " << beamStart << ");" <<endl;
   //      cout << "( " << x2 << " , " << beamStop << " ); slope= " << m << "; intercept= " << c << endl;
   //      cout << "Two straws closest to that point are " << lower << ", and " << upper <<  "; with DCAs "<< hit_distance_low <<  " and " << hit_distance_up << ", respectively." << endl;
   //      cout << "Selected DCA as the correct hit distance is " << hit_distance << endl;
   //  }
    
   //      return hit_distance; // as dca of the closest straw
    return 0;
}


DCAData Tracker::DCAHit_simple(std::vector<float> layer_x, float zStraw, float xTrack, bool debugBool){ 

    DCAData dca_data;

    bool StrongDebugBool=false; //quick hack XXX for even more debug output

    //Find the two closest x straw coordinates given x on the line - with same z [the vector of straw x is naturally in an ascending order]
   //TODO ensure this also true for >= and <= for hits going through the actual wire!! 
   // if hit bigger that all straws all smaller  
    
    float lower, upper, hitDistance;
    int lastID = strawN-1; // the ID of the very last straw in the vector
    float comparator; // to find the ID
    float LR; //L or R hit

    if (debugBool && StrongDebugBool){
        //cout << "Track (line) point in-line with the layer at ( " << x_line << " , "  << z_straw << " ); belongs to the line with  generated points: ( " << x_0 <<" , " << beamStart << ");" <<endl;
        //cout << "( " << x_0 << " , " << beamStop << " );" << endl;
    }

    //Iterator to find two straws between the hit:
    vector<float>::iterator it;
    it = lower_bound(layer_x.begin(), layer_x.end(), xTrack);
    // no smaller value  than val in vector
    if (it == layer_x.begin()){
        upper = *it; 
        //if (debugBool && StrongDebugBool){ cout << " upper = " << upper << endl;}
         // hit distance is the dca from the L of the straw 0
        hitDistance = Tracker::DCA_simple(upper, xTrack);
        comparator = upper;
        LR=-1.0;
        if (debugBool && StrongDebugBool && (hitDistance > strawSpacing/2)){
            cout << "Hit at " << xTrack << " The first straw is closest at " << upper <<  "; with DCA "<< hitDistance << endl;
        
        }
    }   
    // no bigger value than val in vector
    else if (it == layer_x.end()){
        lower = *(it-1); 
        //if (debugBool && StrongDebugBool){ cout << "lower = " << lower << endl;}
        // hit distance is the dca from the R of the last straw
        hitDistance = Tracker::DCA_simple(lower, xTrack); // hit distance is the dca from the L of the straw 0
        comparator = lower;
        LR=1.0;
        if (debugBool && StrongDebugBool && (hitDistance > strawSpacing/2)){
            cout << "Hit at " << xTrack << " The last straw is closest at " << lower <<  "; with DCA "<< hitDistance << endl;
            
        }
    }
    else {
        lower = *(it-1);
        upper = *it;
        //if (debugBool && StrongDebugBool){ cout << "lower = " << lower << " upper = " << upper << endl;}

            float hit_distance_low = Tracker::DCA_simple(lower, xTrack);
            float hit_distance_up = Tracker::DCA_simple(upper, xTrack); 

            if (hit_distance_low < strawRadius && hit_distance_up < strawRadius){
                if (debugBool){cout << "multiple Straw in layer were hit!" << endl;}
                multipleHitsLayer++;
            }

            if (hit_distance_up>hit_distance_low){
                hitDistance = hit_distance_low;
                comparator = lower;
                LR=+1.0;
            }
            if (hit_distance_up<hit_distance_low){
                hitDistance = hit_distance_up;
                comparator = upper;
                LR=-1.0;
            }
            //XXX resolve this [won't be an issue once physical geometry and hit rejection is added]
            if (hit_distance_up==hit_distance_low){
                hitDistance = hit_distance_low; cout << "Ambiguity which straw registered hit" << endl;
                exit(0);
            } 
            if (debugBool && StrongDebugBool && 1==0){
                cout <<  "Hit at " << xTrack << "Two straws closest to that point are " << lower << ", and " << upper <<  "; with DCAs "<< hit_distance_low <<  " and " << hit_distance_up << ", respectively." << endl; 
            }
    } // end of iterator to find straws between hits
   
    //Iterator to find straw ID [implemented separately, as now we look up once only]:
    int index;
    it = find(layer_x.begin(), layer_x.end(), comparator);
    if (it == layer_x.end()){   
        cout << "Error: straw x not found in vector!" << endl;
        exit(0);
      
    } else {
        index = std::distance(layer_x.begin(), it);
        if (index > lastID || index < 0){
            cout << "Error: ID is out of range" << endl;
            exit(0);
        }
        dca_data.strawID = index;
    }

    dca_data.dca = hitDistance;
    dca_data.LRSign = LR;

    if (debugBool && StrongDebugBool && 1==0){
       
        cout << "Selected DCA as the correct hit distance is " << hitDistance << ". Straw ID: " << dca_data.strawID;
        if (LR  > 0){
            cout << ". The straw was hit from the right" << endl;
        }
        if (LR < 0) {
            cout << ". The straw was hit from the left" << endl;
        }
    }//debug
        return dca_data; // as dca and id of the closest straw
}

// Adds dca to the ideal geometry 
float Tracker::GetIdealPoint(int x_det_ID, float x_det_dca, float LRSign, vector<float> layer_x){
    
    float x_IdealStraw=layer_x[x_det_ID];

    float x_IdealPoint = x_IdealStraw + (LRSign * x_det_dca); //add dca for R 

    //cout << "x_IdealPoint= " << x_IdealPoint << "; x_IdealStraw= " << x_IdealStraw << " LRSign = " << LRSign  <<  " x_det_dca= " << x_det_dca << endl;

    return x_IdealPoint;
}

// Function to return residuals to the fitted line (due to dca point scatter + resolution of the detector)
// @ Inputs: ideal points (x,z)
// Algorithm based on SLR: https://en.wikipedia.org/wiki/Simple_linear_regression#Fitting_the_regression_line
vector<float> Tracker::GetResiduals(vector<float> IdealPoints, vector<float> distanceZ){

    //XXX Hack: for now this is just the smearing
    //TODO implement best fit line
    vector<float> residuals; //to store residuals after fit
    for (int i_size=0; i_size<IdealPoints.size(); i_size++){
        //RESIDUAL between a point (dca on a straw due to misalignment + ideal position) and detected position *
        float residual_i=IdealPoints[i_size]; 
    
        float rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
        residual_i = residual_i + rand_gaus * resolution; // smeared resolution 
        residuals.push_back(residual_i);
    }
    return residuals;
}

// Function to return residuals to the fitted line to a parallel tracks: SIMPLE: sum(x)/N(x)
// @ Inputs: ideal points (x,z)
// Algorithm based on SLR: https://en.wikipedia.org/wiki/Simple_linear_regression#Fitting_the_regression_line
ResidualData Tracker::GetResiduals_simple(vector<float> IdealPoints, ofstream& plot_fit){

    ResidualData resData;
    float SumX=0;
    for (int i_size=0; i_size<IdealPoints.size(); i_size++){
        //RESIDUAL between a point (dca on a straw due to misalignment + ideal position) and detected position *
        SumX+=IdealPoints[i_size];
        
    }
   
    float x_fit=SumX/IdealPoints.size(); 
    plot_fit << x_fit << endl;
    
    //cout << "x_ftt =" << x_fit; 
    for (int i_size=0; i_size<IdealPoints.size(); i_size++){
        //RESIDUAL between a point (dca on a straw due to misalignment + ideal position) and detected position *
        float residual_i = abs(x_fit-IdealPoints[i_size]);
        resData.residuals.push_back(residual_i);
        resData.x_fitted.push_back(x_fit);
    }

    return resData;
}
    


/**
   Simulate the passage of a randomly generated linear track through the detector, calculating properties of hits by the track on detector planes.
  
   @return LineData struct containing data about detector hits.
*/
MCData Tracker::MC(float scatterError, ofstream& debug_calc, ofstream& debug_off, ofstream& debug_mc, ofstream& plot_fit, bool debugBool) {
  
    // Set up new container for track data, with hit count set to zero
    MCData MC;
    MC.hit_count = 0; //count only counts through detector
    //vector to store x coordinates of the track as seen from the ideal detector 
    vector<float> x_idealPoints;  

    float rand_num; 
    float rand_gaus;

   // Track parameters for rand-generated line MC [start and end positions outside of detectors]
    rand_num = (( RandomBuffer::instance()->get_uniform_number()+ RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max()));
    float x0 = beamPositionLength * (rand_num)+beamOffset; //uniform vertex
    //rand_num = ((RandomBuffer::instance()->get_uniform_number() + RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max()));
    //float x1 = beamPositionLength * (rand_num)+beamOffset; //uniform exit point: so fitting a line to these two points

    //float x_slope=(beamStop-beamStart)/(x_1-x_0); //m=dz/dx? Slope for a line with coordinates (x_0, beamStart), (x_1, beamStop)
    //float x_intercept = beamStart - x_slope*x_0; // for line z=mx+c -> c= z - mx

    //float xTrack; //true track position x=(z_straw-c)/m; XXX always x_0 for parallel track 
    //float dx = x_slope;
    float previousPosition = 0.0;  // previous position in "z"
    
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
        xTrack=xTrack+dx*dPostion;  
        line.x_true.push_back(xTrack);
        //multiple scattering
        rand_gaus= RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
        dx = dx+ rand_gaus * scatterError;
    #endif

        //loops over views, layers, straws 
       // for (int i_view=0; i_view<viewN; i_view++){
        for (int i_layer=0; i_layer<layerN; i_layer++){ //per layer
            //for (int i_straw=0; i_straw<strawN; i_straw++){

                //xTrack = (distance[z_counter]-x_intercept)/x_slope;   // true track position [from line] WRONGG.... XXX need dca here
                float xTrack = x0; // always the same for parallel line 
                //the dca between xTrack and the xStraw [misaligned]
                DCAData x_MisDetector= Tracker::DCAHit_simple(mod_lyr_strawMisPosition[i_module][i_layer], distance[z_counter], xTrack, debugBool); // position on the detector [from dca]           
                float x_mis_dca =  x_MisDetector.dca;  //dca of the hit straw
                int x_mis_ID =  x_MisDetector.strawID; // ID of the hit straw [to identify the correct straw in the Fit function]
                float x_mis_LRSign = x_MisDetector.LRSign;
                MC.strawID.push_back(x_mis_ID);
                MC.LR.push_back(x_mis_LRSign);

                if(debugBool){cout << "DCA is= " << x_mis_dca << " for straw ID= " << x_mis_ID << " was hit from " << x_mis_LRSign << endl;}

                //Finding the hit as seen from the ideal detector
                float xIdeal = Tracker::GetIdealPoint(x_mis_ID, x_mis_dca, x_mis_LRSign, mod_lyr_strawIdealPosition[i_module][i_layer]);
                x_idealPoints.push_back(xIdeal); // vector to store x coordinates of the track as seen from the ideal detector
                
                if(debugBool){cout << "Ideal x= " << xIdeal << endl;}

                //Rejection of hits due to geometry (i.e. missed hits)  
                //No signal in a straw = no signal in the layer
                if (x_mis_dca > strawRadius){ //[between (0,0.25]       
                    continue;
                    rejectedHitsDCA++;
                    if (debugBool){cout << "Hit Rejected: outside of straw!" << endl;}
                } 
                
                //Module number [for labelling] - after (if) passing the rejection...
                MC.i_hits.push_back(i_module); // vector of modules that were actually hit [after passing rejection test]
                
                //Z-coordinate of hits [it was always known from geometry - no z-misalignment for now...]
                MC.z_hits.push_back(distance[z_counter]);
                MC.hit_sigmas.push_back(resolution); 
            
            //This was left over from MP2 Toy model: we are implementing dca instead...some useful methods still might be there... XXX
            #if 0 
                
                // the residual looks to be deltaX + deltaY rather than the magnitude of the distance... worth noting?
                // projection Y is always 0 for non-stereo modules?? what is the motivation? 
                //rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
                //float hit = (x_mis-x_det)*projectionX[i_module] + rand_gaus *resolution; //RESIDUAL between hit due to misalignment [] and detected position * smeared resolution 
               // float dca_misStraw = x_det_dca*projectionX[i_module] + rand_gaus *resolution; //RESIDUAL between hit due to misalignment [] and detected position * smeared resolution 
               // line.x_hits.push_back(dca_misStraw);
            #endif
                
                //Sanity Plots: Hits
                MC.x_mis_dca.push_back(x_mis_dca);
                MC.x_track.push_back(xTrack);
                MC.x_ideal.push_back(xIdeal);
                
                //Sanity Plots: Tracks
                if (MC.hit_count == 0){
                    //cout << "Track plot push back" << endl;
                    //MC.x_m.push_back(x_slope);
                    //MC.x_c.push_back(x_intercept);
                    MC.x0_gen.push_back(xTrack);
                    MC.x1_gen.push_back(xTrack);
                    MC.z0_gen.push_back(float(beamStart)); //XXX 
                    MC.z1_gen.push_back(float(beamStop));  //XXX 
                } 

                MC.hit_count++;

           // } //end of Straws loop
        z_counter++;
        }//end of Layers loop
   // }// end of View loop
   
    }// end of looping over modules

    //This happens once per MC function call [as we now accumulated x coordinates of "ideal" point for all hits
    // and need to do a simultaneous fit once - to return #hits of residuals]
    ResidualData res_Data = Tracker::GetResiduals_simple(x_idealPoints, plot_fit);
    vector<float> residuals = res_Data.residuals;
    vector<float> x_fitted=res_Data.x_fitted; 

    ///Here we will push back into "x_hits" for loop for 0-> hit_count to unpack in main in the hit loop
    for (int i_HitCounter=0; i_HitCounter<MC.hit_count; i_HitCounter++){
        MC.x_hits.push_back(residuals[i_HitCounter]); 
        MC.x_fitted.push_back(x_fitted[i_HitCounter]); //Sanity Plot 
    }

    return MC; // Return data from simulated track
    
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
   
   // Print out:
    if (debugBool){
        int Zcounter=0;
        for (int i_module=0; i_module<moduleN; i_module++){
          //for (int i_view=0; i_view<viewN; i_view++){
                for (int i_layer=0; i_layer<layerN; i_layer++){ //per module
                    cout << "IDEAL Mod " << i_module  << " Layer " << i_layer << " X : ";
                    for (int i_straw=0; i_straw<strawN; i_straw++){
                    cout << mod_lyr_strawIdealPosition[i_module][i_layer][i_straw]<<" ";
                    debug_geom << mod_lyr_strawIdealPosition[i_module][i_layer][i_straw] << " ";
                    } // straws 
                    cout << " | Z= " << distance[Zcounter] << " [cm]" << endl;  // XXX fix this 
                    debug_geom << distance[Zcounter] << endl;
                    Zcounter++;
                } // end of Layers
                
           // }// end of View loop
        }//Modules 
    }


} // end of geom

// MC misalignment of detectors 
void Tracker::misalign(ofstream& debug_mis, bool debugBool){
    //float rand_gaus;
    //float sign = 1.0; //for +x or -x direction of misalignment
    float factor = 1.0; // to enlarge the effect of misalignment
    //Now misaligning detectors in x

    float dX = startingXDistanceStraw0+(dispX * factor); // starting on the x-axis (z, 0+dispX)

    for (int i_module=0; i_module<moduleN; i_module++){
        //rand_gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev()); 
        // Before misalignment was also smeared XXX
        //sdevX[i] = dispX * rand_gaus;
        sdevX.push_back(dispX * factor);  // vector of misalignment 
        mod_lyr_strawMisPosition.push_back(vector<vector<float> >()); //initialize the first index with a 2D vector
         //for (int i_view=0; i_view<viewN; i_view++){
                for (int i_layer=0; i_layer<layerN; i_layer++){ //per module
                    mod_lyr_strawMisPosition[i_module].push_back(vector<float> ()); //initialize the first index with a 2D vector
                    for (int i_straw=0; i_straw<strawN; i_straw++){
                        mod_lyr_strawMisPosition[i_module][i_layer].push_back(dX); 
                        dX=strawSpacing+dX; //while we are in the same layer: increment straw spacing in x
                    } //end of Straws loop
                    dX=layerDisplacement+(dispX * factor); //set displacement in x for the next layer in the view
                }//end of Layers loop
                //sign = -sign;  //next module will be displaced in the opposite direction
                factor = 1.25;  //next module will be displaced 1.25 more in the same direction 
                dX=startingXDistanceStraw0+(dispX * factor);
           // }// end of View loop
        
        }//Modules 

    // Print out:
    if (debugBool){
        int Zcounter=0;
        for (int i_module=0; i_module<moduleN; i_module++){
          //for (int i_view=0; i_view<viewN; i_view++){
                for (int i_layer=0; i_layer<layerN; i_layer++){ //per module
                    cout << "MIS Mod " << i_module  << " Layer " << i_layer << " X : ";
                    for (int i_straw=0; i_straw<strawN; i_straw++){
                    cout << mod_lyr_strawMisPosition[i_module][i_layer][i_straw]<<" ";
                    debug_mis << mod_lyr_strawMisPosition[i_module][i_layer][i_straw] << " ";
                    
                    } //end of Straws loop
                    cout << " | Z= " << distance[Zcounter] << " [cm]" << endl;  // XXX fix this 
                    debug_mis << distance[Zcounter] << endl;
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

/**
 * Returns the peak (maximum so far) resident set size [RSS] (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS. By Dr. David R. Nadeau:
 * http://nadeausoftware.com/articles/2012/07/c_c_tip_how_get_process_resident_set_size_physical_memory_use
 */
size_t Tracker::getPeakRSS( ){
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
        return (size_t)0L;      /* Can't open? */
    if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
    {
        close( fd );
        return (size_t)0L;      /* Can't read? */
    }
    close( fd );
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;          /* Unsupported. */
#endif
}