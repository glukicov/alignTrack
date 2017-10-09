/*
This source file contains definitions of various functions in the method class. 
*/  

#include "AlignTracker_methods.h"

using namespace std;

// Set up empty pointer for instance of the class.
Tracker* Tracker::s_instance = NULL;

/** 
Constructor for tracker class.
 */
Tracker::Tracker() {
    // non-static variables definition here

    // Set mapping for U0...V1
    string tempMapping[4] = {"U0", "U1", "V0", "V1"};
   	for (int i_view=0; i_view<viewN; i_view++){
        UVmapping.push_back(vector<string> ()); //initialize the first index with a 2D vector
          for (int i_layer=0; i_layer<layerN; i_layer++){
             if (i_view==0){UVmapping[i_view].push_back(tempMapping[i_layer]);}
             if (i_view==1){UVmapping[i_view].push_back(tempMapping[i_layer+2]);}
         } //layers
    } // views      
} // constructor 

/** 
    Empty destructor for tracker class.
*/
Tracker::~Tracker() {
}


/**
   Get pointer to only instance of tracker class, creating this instance if it doesn't already exist.
   instance is a static member function of Tracker (singleton)
   @return Pointer to tracker instance.
 */
Tracker* Tracker::instance() {

    // Create pointer to class instance if one doesn't exist already, then return that pointer.
    if (s_instance == NULL) s_instance = new Tracker();
    return s_instance;
}

void Tracker::write_steering_file(ofstream& steering_file) {
    if(steering_file.is_open()){
        steering_file <<  "* g-2 Tracker Alignment: PEDE Steering File" << endl
        << " "  << endl
        << "Tracker_con.txt   ! constraints text file " << endl
        << "Cfiles ! following bin files are Cfiles" << endl 
        << "Tracker_data.bin   ! binary data file" << endl
        << "method inversion 1 0.01" << endl
        << "printrecord  -1 -1      ! debug printout for bad data records" << endl
        << "printrecord 1 -1 ! produces mpdebug.txt"<< endl     //  
        << " "  << endl
        << "end ! optional for end-of-data"<< endl;
    } // steering file open 
} // steering function
/**
   Write a constraint file to the supplied file-stream.
    @param constraint_file Reference to ofstream to write constraint file to. 
 */
// XXX constraints are ignored with HIP method 
void Tracker::write_constraint_file(ofstream& constraint_file, ofstream& debug_con, bool debugBool) {
    // Check constraints file is open, then write. 
    if (constraint_file.is_open()) {
        //Evaluation of constraints
        float one = 1.0;
        //Fixing module 0 and the last module
        for (int i_module = 0; i_module < moduleN; i_module++){ 
            if (i_module==0 || i_module==moduleN-1){
                constraint_file << "Constraint 0.0" << endl;
                int labelt=i_module+1; // Millepede accepts +ive labels only
                constraint_file << labelt << " " << fixed << setprecision(5) << one<< endl;
            } // end of fixed modules
        } // end of detectors loop 
    } // constrain file open 
} // end of writing cons file 

/**
   DCA for a point to a line in 2D space
*/
float Tracker::pointToLineDCA(float z_straw, float x_straw, float x_slope, float x_intercept){

	//http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html

	//converting from x=mz+c -> mz+x+c=0
	float a = -x_slope;
    float b = 1.0; 
    float c = -x_intercept;
    float x0 = z_straw;
    float y0 = x_straw;

    float dca = abs( a*x0 + b*y0 + c ) / ( a*a + b*b  ) ;

    return dca;
}


/** Uses DCA function to find the shortest dca between straws (i.e. which straw was hit in that layer)
    @Inputs (see DCA method) + vector of straws' x coordinates in a layer, and a return type: "dca_hit"  or "x_line"
    @return hit_distance - dca of the straw that was hit, or x coordinate of the line on the same z as a straw
*/
DCAData Tracker::DCAHit(vector<float> xLayer, float zStraw, float xTrack, float xSlpoe, float xIntercept, bool debugBool){

    DCAData dca_data;

    bool StrongDebugBool=false; //quick hack XXX for even more debug output

    //Find the two closest x straw coordinates given x on the line - with same z [the vector of straw x is naturally in an ascending order]
      
    float lower, upper, hitDistance;
    int lastID = strawN-1; // the ID of the very last straw in the vector
    float comparator; // to find the ID
    float LR; //L or R hit looking downstream

    //Iterator to find two straws between the hit:
    vector<float>::iterator it;
    it = lower_bound(xLayer.begin(), xLayer.end(), xTrack);
    // no smaller value  than val in vector
    if (it == xLayer.begin()){
        upper = *it; 
        if (debugBool && StrongDebugBool){ cout << " upper = " << upper << endl;}
         // hit distance is the dca from the L of the straw 0
        hitDistance = pointToLineDCA(zStraw, upper, xSlpoe, xIntercept);
        comparator = upper;
        LR=-1.0;
        if (debugBool && StrongDebugBool && (hitDistance > strawSpacing/2)){
            cout << "Track in line at " << xTrack << " The first straw is closest at " << upper <<  "; with DCA "<< hitDistance << endl;
        }
    }   
    // no bigger value than val in vector
    else if (it == xLayer.end()){
        lower = *(it-1); 
        if (debugBool && StrongDebugBool){ cout << "lower = " << lower << endl;}
        // hit distance is the dca from the R of the last straw
        hitDistance = pointToLineDCA(zStraw, lower, xSlpoe, xIntercept); // hit distance is the dca from the L of the straw 0
        comparator = lower;
        LR=1.0;
        if (debugBool && StrongDebugBool && (hitDistance > strawSpacing/2)){
            cout << "Track in line at " << xTrack << " The last straw is closest at " << lower <<  "; with DCA "<< hitDistance << endl;
        }
    }
    else {
        lower = *(it-1);
        upper = *it;
        if (debugBool && StrongDebugBool){ cout << "lower = " << lower << " upper = " << upper << endl;}

            float hit_distance_low = pointToLineDCA(zStraw, lower, xSlpoe, xIntercept);
            float hit_distance_up = pointToLineDCA(zStraw, upper, xSlpoe, xIntercept); 

            if (hit_distance_low < strawRadius && hit_distance_up < strawRadius){
                if (debugBool){cout << "Multiple straws in layer were hit!" << endl;}
                incMultipleHitsLayer();
            }

            // if DCA in higher straw ID is bigger, select straw ID with smaller DCA
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
            
            if (hit_distance_up==hit_distance_low){
                cout << "BIAS" << endl;
                float random = Tracker::generate_uniform();
                if (random<0.5){
                    hitDistance = hit_distance_low;
                }
                if (random>0.5){
                    hitDistance = hit_distance_up;
                }
                if (1==1){cout << "Ambiguity which straw registered hit" << endl;}
                cout << "Ambiguity which straw registered hit" << endl;
                incAmbiguityHit();
                //exit(0);
            } 
            if (debugBool && StrongDebugBool && 1==0){
                cout <<  "Track in Line " << xTrack << "Two straws closest to that point are " << lower << ", and " << upper <<  "; with DCAs "<< hit_distance_low <<  " and " << hit_distance_up << ", respectively." << endl; 
            }
    } // end of iterator to find straws between hits
   
    //Iterator to find straw ID [implemented separately, as now we look up once only]:
    int index;
    it = find(xLayer.begin(), xLayer.end(), comparator);
    if (it == xLayer.end()){   
        cout << "Error: straw x not found in vector!" << endl;
        exit(0);
    }
    else {
        index = std::distance(xLayer.begin(), it);
        if (index > lastID || index < 0){
            cout << "Error: ID is out of range" << endl;
            exit(0);
        }
    dca_data.strawID = index;
    }

    dca_data.dcaUnsmeared = hitDistance;
    float hitDistanceSmeared=hitDistance+Tracker::instance()->getResolution()*Tracker::generate_gaus();
    dca_data.dca = hitDistanceSmeared;
    float residualTruth= hitDistanceSmeared - hitDistance;
	dca_data.residualTruth = residualTruth;
    dca_data.LRSign = LR;

    if (debugBool && StrongDebugBool && 1==0){
        cout << "Selected DCA as the correct hit distance is " << hitDistance << ". Straw ID: " << dca_data.strawID;
        if (LR  > 0){cout << ". The straw was hit from the right" << endl;}
        if (LR < 0) {cout << ". The straw was hit from the left" << endl;}
    }//debug
    return dca_data; // as dca and id of the closest straw
}


// Adds dca to the ideal geometry 
ReconData Tracker::HitRecon(int det_ID, float det_dca, vector<float> xLayer, float z_distance){

	ReconData recon_data;
    
    bool StrongDebugBool=false; //quick hack XXX for even more debug output
    float x_IdealStraw=xLayer[det_ID];
    float recon_dca = det_dca; //adds dca to the assumed straw x coordinate  
    if (StrongDebugBool){ cout << "recon_dca= " << recon_dca << "; x_IdealStraw= " << x_IdealStraw <<  " det_dca= " << det_dca << endl;}
    recon_data.z = z_distance; 
    recon_data.x = x_IdealStraw;
    recon_data.dcaRecon = recon_dca; 
    return recon_data;
}

// Function to return residuals to the fitted line (due to dca point scatter + resolution of the detector)
// @ Inputs: ideal points (x,z)

//Truth can be used as an optional input 
  
ResidualData Tracker::GetResiduals(vector<float> zRecon, vector<float> xRecon, vector<float> radRecon, int dataSize, ofstream& plot_fit, bool debugBool, bool useTruth){

    ResidualData resData;

    bool StrongDebugBool = false; 

    //Circle Fit goes here; TODO

           
    
    //resData.residuals.push_back(Res);   // residual between the (centre of the straw and the fitted line [pointToLineDCA]) and radius of the fit circle; 
    //resData.slope_recon=slope;
    //resData.intercept_recon=intercept;
    // resData.meanXReconTrack=AVGy;
    // resData.meanZReconTrack=AVGx;

    //if (debugBool){ plot_fit <<  slope*beamStart+y_intercept << " "  << slope*beamStop+y_intercept  <<   " " <<  beamStart  << " " << beamStop << endl; }
    return resData;
}


/**
   Simulate the passage of a randomly generated linear track through the detector, calculating properties of hits by the track on detector planes.
   @return MCData struct containing data about detector hits.
*/
MCData Tracker::MC_launch(float scatterError, ofstream& debug_calc, ofstream& debug_off, ofstream& debug_mc, ofstream& plot_fit, ofstream& plot_gen, ofstream& plot_hits_gen, ofstream& plot_hits_fit, bool debugBool) {
    
    bool StrongDebugBool = false; // HACK XXX set by hand for now
    // Set up new container for track data, with hit count set to zero
    MCData MC;
    MC.hit_count = 0; //count only counts through detector
    hitLayerCounter=-1; // start counting from 0, increment before hit-rejection
    
    //clearing containers for next track 
    xRecon.clear();  // ideal straw x
    zRecon.clear();  // ideal straw z
    radRecon.clear();
    
    //Track parameters for rand-generated line MC [start and end positions outside of detectors]
    // redefining the track as x=ym+c
    float x0 = beamPositionLength * Tracker::generate_uniform()-1.0; //uniform vertex
    float signXSlope;
    if (Tracker::generate_uniform() >= 0.5){
    	signXSlope=1.0;
    }
    else{
    	signXSlope=-1.0;
    }  
    float xSlope = (Tracker::generate_uniform()*signXSlope) * (0.5*beamPositionLength/beamStop); 
    float xIntercept =x0; // by definition 
    float x1 = xSlope*beamStop + xIntercept;

    float xTrack; //true track position x=zm+c 

    MC.x0 = x0;
    MC.x1 = x1;
    MC.slope_truth = xSlope;
    MC.intercept_truth = xIntercept;
        
    //The main loop is for modules [they produce label of Global Parameters]: 
    // Then looping over layers and views
    // within layers we will have loops over straws in that layer [to find out which straw was hit and dca]
    int z_counter=0; // distance in z is incremented from 0 to TotalLayerN 
    if (debugBool && StrongDebugBool){cout << "Calculating hits/dca:" << endl;}
    //loops over modules, views, layers
    for(int i_module=0; i_module<moduleN; i_module++){    
    	for (int i_view=0; i_view<viewN; i_view++){
        	for (int i_layer=0; i_layer<layerN; i_layer++){ //per layer
                hitLayerCounter++;
         
                // TRUTH PARAMETERS
                //The registered hit position on the misaligned detector is smeared by its resolution
                // zTrack is just the distance[z_counter] - in-line with the layer
                xTrack = xSlope*distance[z_counter]+xIntercept;  // true track position in-line with the layer [from line x=ym+c]
                MC.x_track_true.push_back(xTrack);  // True in-line (gen.) track position
                

                //float xHit=xTrack+Tracker::instance()->getResolution()*Tracker::generate_gaus();
                //float residualTrack= xHit - xTrack;
                   
	            //the DCA between the track and the (zStraw, xStraw) [misaligned]
                // DCA will return the pointToLineDCA (smeared by resolution) and strawID [+truth parameters: LR sign etc. xHit and zHit]

	            DCAData MisDetector= Tracker::DCAHit(mod_lyr_strawMisPosition[i_module][i_view][i_layer], distance[z_counter], xTrack, xSlope, xIntercept, debugBool); // position on the detector [from dca]  
	            float mis_dca =  MisDetector.dca;  //dca (smeared) of the hit straw
	      		      
	            //Rejection of hits due to geometry (i.e. missed hits)  
	            //No signal in a straw = no signal in the layer
	            // if (mis_dca > 0.5*strawSpacing){ //[between (0,0.25] 
	            //     MC.hit_bool.push_back(0);       
	            //     incRejectedHitsDCA();
	            //     if (debugBool){cout << "Hit Rejected: outside of straw layer with dca =" << mis_dca << endl;}
	            //     stringstream absolute_hit;
	            //     absolute_hit <<"#"; 
	            //     MC.absolute_straw_hit.push_back(absolute_hit.str().c_str());
	            //	   z_counter++;  // incrementing distance of planes   
	            //     continue;
	            // } 

	            //Find the truth ID, and LR hit for a straw    
	            float dcaUnsmeared = MisDetector.dcaUnsmeared;
				float residualTruth = MisDetector.residualTruth;
	            int mis_ID =  MisDetector.strawID; // ID of the hit straw [to identify the correct straw in the Fit function]
	            float mis_LRSign = MisDetector.LRSign;
	            MC.strawID.push_back(mis_ID);
	            MC.LR.push_back(mis_LRSign);
                //Recording hit information
	            MC.hit_list.push_back(hitLayerCounter); //push back the absolute layer # into the hit list
	            MC.hit_bool.push_back(1);  // 1 = hit
	            stringstream absolute_hit;
	            //Making absolute strawID for hit for plotting:
	            absolute_hit << mis_ID;
	            MC.absolute_straw_hit.push_back(absolute_hit.str().c_str());         
	            if(debugBool && StrongDebugBool){cout << "DCA is= " << mis_dca << " for straw ID= " << mis_ID << " was hit from " << mis_LRSign << endl;}
	            if (debugBool){plot_hits_gen << mod_lyr_strawMisPosition[i_module][i_view][i_layer][mis_ID] << " " << distance[z_counter] << " "  << mis_dca << endl;}
	            // RECONSTRUCTED PARAMETERS
                //Reconstructing the hit as seen from the ideal detector [shift circle x coordinate by misalignment]

	            //TODO: remove LRSign, reconstruction returns 

	            ReconData ReconDetector = Tracker::HitRecon(mis_ID, mis_dca, mod_lyr_strawIdealPosition[i_module][i_view][i_layer], distance[z_counter]);
	            float zCircle = ReconDetector.z;
	            float xCircle = ReconDetector.x;
	            float radCircle = ReconDetector.dcaRecon;
	            xRecon.push_back(xCircle); // vector to store x coordinates of circles as seen from the ideal detector
	            zRecon.push_back(zCircle); // vector to store z coordinates of circles as seen from the ideal detector
	            radRecon.push_back(radCircle); // vector to store radius of circles as seen from the ideal detector

	            if (debugBool){plot_hits_fit << xCircle << " " << zCircle << " "  << radCircle << endl;}
	            //Module number [for labelling] - after (if) passing the rejection.
                // Millepede accepts only positive non-zero integers as labels
	            ostringstream oss; oss << i_module + 1; 
	            istringstream iss(oss.str()); int label_int; iss >> label_int;
	            MC.i_hits.push_back(label_int); // vector of modules that were actually hit [after passing rejection test: for MP2 labelling]
	      
	            //Z-coordinate of hits [it was always known from geometry - no z-misalignment for now...]
	            MC.z_hits.push_back(distance[z_counter]);
	            MC.hit_sigmas.push_back(Tracker::instance()->getResolution()); 
	                
	            //Sanity Plots: Hits
	            MC.mis_dca.push_back(radCircle); // DCA
	            MC.residuals_gen.push_back(residualTruth); // Truth residual
                MC.Module_i.push_back(i_module);
                MC.View_i.push_back(i_view);
                MC.Layer_i.push_back(i_layer);
                MC.Straw_i.push_back(mis_ID);
	           
                MC.hit_count++;
	            z_counter++;  // incrementing distance of planes
	        }//end of Layers loop
        }// end of View loop
    }// end of looping over modules

    if (debugBool && StrongDebugBool){cout << "Calculating residuals..." << endl;}
    //This happens once per MC function call [as we now accumulated x coordinates of "ideal" points for all hits
    // and need to do a simultaneous fit once - to return #hits worth of residuals]
    //ResidualData res_Data = Tracker::GetResiduals(xReconPoints, distance, plot_fit, debugBool);
    
   

    // XXX HACK for now
    bool useTruth = false;



    ResidualData res_Data = GetResiduals(zRecon, xRecon, radRecon, MC.hit_count, plot_fit, debugBool, useTruth);
    MC.residuals = res_Data.residuals;
    MC.slope_recon = res_Data.slope_recon;
    MC.intercept_recon = res_Data.intercept_recon;
    MC.meanXReconTrack=res_Data.meanXReconTrack;
    MC.meanZReconTrack=res_Data.meanZReconTrack;
    
    if(debugBool){
                        plot_gen << x0 << " " << x1 << " " << beamStart << " " << beamStop << " " << endl;
                   
                    }

    return MC; // Return data from simulated track
    
} // end of MC

//Geometry of detector arrangement (Ideal Geometry)
void Tracker::setGeometry(ofstream& debug_geom,  bool debugBool){
   // float dZ=startingZDistanceStraw0; // to add distance in z
    
    float dZ=startingZDistanceStraw0; // the increment in z for consecutive layers
    int layer_n = 0;  // layer label 
    for (int i_module=0; i_module<moduleN; i_module++){
    	for (int i_view=0; i_view<viewN; i_view++){
    		for (int i_layer=0; i_layer<layerN; i_layer++){
    			distance.push_back(dZ); // vector will contain all z coordinates of layers 
    			layer.push_back(layer_n);  // layer label array [starting from 0th layer]
				//resolutionLayer.push_back(Tracker::instance()->getResolution()); //resolution in each layer
				projectionX.push_back(float(1.0));  // x projection of hits in each layer
				layer_n++; 
    			if (i_layer==0){ dZ+=layerSpacing; } // increment spacing between layers in a view once only 
    		} // layers
    		if (i_view==0) { dZ+=viewSpacing; } // increment spacing between views in a module once only 
    	} // views
    	dZ+=moduleSpacing;
    } // modules

    // Geometry of detector arrangement in X 
    float dX = startingXDistanceStraw0; // starting on the x-axis (z, 0)
    for (int i_module=0; i_module<moduleN; i_module++){
        mod_lyr_strawIdealPosition.push_back(vector<vector<vector<float> > >()); //initialize the first index with a 2D vector
         for (int i_view=0; i_view<viewN; i_view++){
         	 mod_lyr_strawIdealPosition[i_module].push_back(vector<vector<float> >()); //initialize the first index with a 2D vector
                for (int i_layer=0; i_layer<layerN; i_layer++){ 
                    mod_lyr_strawIdealPosition[i_module][i_view].push_back(vector<float> ()); //initialize the first index with a 2D vector
                    for (int i_straw=0; i_straw<strawN; i_straw++){
                        mod_lyr_strawIdealPosition[i_module][i_view][i_layer].push_back(dX); 
                        dX=strawSpacing+dX; //while we are in the same layer: increment straw spacing in x
                    } //end of Straws loop
            if (i_view==0){ dX=layerDisplacement+startingXDistanceStraw0; } //set displacement in x for the next layer in the view
            if (i_view==1){ dX=startingXDistanceStraw0; } //set displacement in x for the next layer in the view
            }//end of Layers loop
        }// end of View loop
    }//Modules  
 
   // Print out:
    if (debugBool){
        stringstream gm1;
        gm1 << endl << "Ideal Detector Position:";
        Logger::Instance()->write(Logger::WARNING, gm1.str());
        int Zcounter=0;
        for (int i_module=0; i_module<moduleN; i_module++){
          for (int i_view=0; i_view<viewN; i_view++){
                for (int i_layer=0; i_layer<layerN; i_layer++){ //per module
                    cout << "IDEAL Mod " << i_module  << " " << UVmapping[i_view][i_layer] << " X : ";
                    for (int i_straw=0; i_straw<strawN; i_straw++){
                    cout << mod_lyr_strawIdealPosition[i_module][i_view][i_layer][i_straw]<<" ";
                    debug_geom << mod_lyr_strawIdealPosition[i_module][i_view][i_layer][i_straw] << " ";
                    } // straws 
                    cout << "  | Z= " << distance[Zcounter] << " [cm]" << endl;  // TODO align the cout better 
                    debug_geom << distance[Zcounter] << endl;
                    Zcounter++;
                } // end of Layers
                
            }// end of View loop
        cout << endl;
        }//Modules 
    }

} // end of geom

// MC misalignment of detectors 
void Tracker::misalign(ofstream& debug_mis, ofstream& pede_mis, bool debugBool){
     
    //Now misaligning detectors in x
    float misDispX; // effective misalignment 
    for (int i_module=0; i_module<moduleN; i_module++){
    	misDispX=dispX[i_module];
               
        float dX = startingXDistanceStraw0+(misDispX); // starting on the x-axis (z, 0+disp)
    	    	sdevX.push_back(misDispX);  // vector to store the actual of misalignment
        mod_lyr_strawMisPosition.push_back(vector<vector<vector<float> > >()); //initialize the first index with a 2D vector
         for (int i_view=0; i_view<viewN; i_view++){
         	 mod_lyr_strawMisPosition[i_module].push_back(vector<vector<float> >()); //initialize the first index with a 2D vector
                for (int i_layer=0; i_layer<layerN; i_layer++){ 
                    mod_lyr_strawMisPosition[i_module][i_view].push_back(vector<float> ()); //initialize the first index with a 2D vector
                    for (int i_straw=0; i_straw<strawN; i_straw++){
                        mod_lyr_strawMisPosition[i_module][i_view][i_layer].push_back(dX); 
                        dX=strawSpacing+dX; //while we are in the same layer: increment straw spacing in x
                    } //end of Straws loop
            if (i_view==0){ dX=layerDisplacement+startingXDistanceStraw0+(misDispX); } //set displacement in x for the next layer in the view
            if (i_view==1){ dX=startingXDistanceStraw0+(misDispX); } //set displacement in x for the next layer in the view
            }//end of Layers loop
        }// end of View loop 
        
    }//Modules  

    // Print out:
    if (debugBool){
        stringstream gm2;
        gm2 << "Misaligned Detector Position:";
        Logger::Instance()->write(Logger::WARNING, gm2.str());
        int Zcounter=0;
        for (int i_module=0; i_module<moduleN; i_module++){
          for (int i_view=0; i_view<viewN; i_view++){
                for (int i_layer=0; i_layer<layerN; i_layer++){ //per module
                    cout << "MIS Mod " << i_module  << " " << UVmapping[i_view][i_layer] << " X : ";
                    for (int i_straw=0; i_straw<strawN; i_straw++){
                    cout << mod_lyr_strawMisPosition[i_module][i_view][i_layer][i_straw]<<" ";
                    debug_mis << mod_lyr_strawMisPosition[i_module][i_view][i_layer][i_straw] << " ";
                    
                    } //end of Straws loop
                    cout << " | Z= " << distance[Zcounter] << " [cm]" << endl;  
                    debug_mis << distance[Zcounter] << endl;
                    Zcounter++;
                } // end of Layers
                
            }// end of View loop
        cout << endl;
        }//Modules 
    }

// Estimating misalignment parameters from geometry and assumed constants:
   
    //--- For Modules only [MC misalignment set per module -> transfer to each layer] -- //

    //Overall Misalignment [w.r.t to other modules]
    float sum_of_elems=0.0;
    for(vector<float>::iterator it = sdevX.begin(); it != sdevX.end(); ++it){
        sum_of_elems += *it;
    }
    overallMis=sum_of_elems/moduleN;
    cout << "Manual Misalignment: " << endl;
    for (int i_module=0; i_module<moduleN; i_module++){
        cout << showpos << "Module " << i_module <<" :: Characteristic:  " << sdevX[i_module] << " cm. ";  // absolute misalignment [as set by MC]   
        if (debugBool){pede_mis <<sdevX[i_module] << " "; }
        float relMisTmp = sdevX[i_module]-overallMis;
        // now push these misalignment parameters for all layers in the module [for use later]
        for (int i_view=0; i_view<viewN; i_view++){
            for (int i_layer=0; i_layer<layerN; i_layer++){
            relMis.push_back(relMisTmp);           //relative misalignment per layer [w.r.t to other modules]
            charMis.push_back(sdevX[i_module]);  // absolute misalignment per layer
            } // view
        } // layer
        cout << showpos << "Relative: " << relMisTmp << " cm." << endl;
    } // modules
    cout << noshowpos; 
    cout << "The overall misalignment was " << overallMis << " cm" <<  endl << endl;

    //--- Calculations for Each Layer -- //
    // XXX Verbose calculation via separate loops on purpose - to demonstrate the alignment methods 
    
    // constants
    float N = float(layerTotalN); // total N of layers
    float D = Tracker::instance()->getResolution(); // Original detector resolution
    int i_totalLayers; // dummy layer counter 
           
    // Pivot Point
    i_totalLayers=0;
    for (int i_module=0; i_module<moduleN; i_module++){
        for (int i_view=0; i_view<viewN; i_view++){
            for (int i_layer=0; i_layer<layerN; i_layer++){
                pivotPoint_estimated+=distance[i_totalLayers];
                i_totalLayers++;
            }// layer
        } // view 
    } // modules
    pivotPoint_estimated = pivotPoint_estimated/N;
    
    // Centred Distances and squared distances 
    i_totalLayers=0;
    for (int i_module=0; i_module<moduleN; i_module++){
        for (int i_view=0; i_view<viewN; i_view++){
            for (int i_layer=0; i_layer<layerN; i_layer++){
            float tmpZDistance= distance[i_totalLayers]-pivotPoint_estimated;
            zDistance_centered.push_back(tmpZDistance);
            squaredZSum += pow(tmpZDistance,2);
            i_totalLayers++;
            }// layer
        } // view 
    } // modules

    // Estimated SD of residuals, Sum of mis., and squared sum of mis.   
    i_totalLayers=0;
    for (int i_module=0; i_module<moduleN; i_module++){
        sigma_recon_estimated.push_back(vector< vector< float > > ()); 
        for (int i_view=0; i_view<viewN; i_view++){
            sigma_recon_estimated[i_module].push_back( vector<float>  () ); 
            for (int i_layer=0; i_layer<layerN; i_layer++){
                float simga_est = D * sqrt( (N-1)/N - (pow(zDistance_centered[i_totalLayers],2)/squaredZSum) );
                sigma_recon_estimated[i_module][i_view].push_back(simga_est);
                MisZdistanceSum += charMis[i_totalLayers]*zDistance_centered[i_totalLayers];
                MisSum += charMis[i_totalLayers];
                i_totalLayers++;
            }// layer
        } // view 
    } // modules

    // Predicted means of residuals [a.k.a "shear misalignment"]
    i_totalLayers=0;
    for (int i_module=0; i_module<moduleN; i_module++){
        for (int i_view=0; i_view<viewN; i_view++){
            for (int i_layer=0; i_layer<layerN; i_layer++){
                float charM = charMis[i_totalLayers];
                float z = zDistance_centered[i_totalLayers];
                float shearMistmp = charM - MisSum/N - (z*MisZdistanceSum)/squaredZSum;
                shearMis.push_back(shearMistmp);
                i_totalLayers++;
            }// layer
        } // view 
    } // modules
   
   // The Chi2 estimation requires misalignment parameters per layer
   Chi2_recon_estimated=N-2.0;  // 2 parameter fit
   i_totalLayers=0;
    for (int i_module=0; i_module<moduleN; i_module++){
        for (int i_view=0; i_view<viewN; i_view++){
            for (int i_layer=0; i_layer<layerN; i_layer++){
                float M = shearMis[i_totalLayers];
                float z = zDistance_centered[i_totalLayers];
                Chi2_recon_estimated += ( (M*M) - (M)/(N) - (2.0*M*z)/(squaredZSum) ) / (D*D) ;
                i_totalLayers++;
            }// layer
        } // view 
    } // modules

}//end of misalign

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

float Tracker::generate_gaus(){
    float gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
    return gaus; 
}

float Tracker::generate_uniform(){
    float uniform = (( RandomBuffer::instance()->get_uniform_number()+ RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max()));
    return uniform; 
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