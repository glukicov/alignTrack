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

}

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

float Tracker::generate_gaus(){
	float gaus = RandomBuffer::instance()->get_gaussian_number() / float(RandomBuffer::instance()->get_gaussian_ran_stdev());
	return gaus; 
}

float Tracker::generate_uniform(){
	float uniform = (( RandomBuffer::instance()->get_uniform_number()+ RandomBuffer::instance()->get_uniform_ran_max()) / (twoR * RandomBuffer::instance()->get_uniform_ran_max()));
	return uniform; 
}

/**
   Simple function to calculate the Distance of Closest Approach
   between a point (x0, y0) and a line (defined by two points: (x1, y1), (x2, y2) ).
  
   @return dca
*/
float Tracker::DCA(float xHit,float xStraw){

    float dca = sqrt((xHit-xStraw)*(xHit-xStraw));

    return dca;
}


/** Uses DCA function to find the shortest dca between straws (i.e. which straw was hit in that layer)

    @Inputs (see DCA method) + vector of straws' x coordinates in a layer, and a return type: "dca_hit"  or "x_line"
  
    @return hit_distance - dca of the straw that was hit, or x coordinate of the line on the same z as a straw
*/
DCAData Tracker::DCAHit(std::vector<float> xLayer, float zStraw, float xHit, bool debugBool){

    DCAData dca_data;

    bool StrongDebugBool=false; //quick hack XXX for even more debug output

    //Find the two closest x straw coordinates given x on the line - with same z [the vector of straw x is naturally in an ascending order]
      
    float lower, upper, hitDistance;
    int lastID = strawN-1; // the ID of the very last straw in the vector
    float comparator; // to find the ID
    float LR; //L or R hit

    //Iterator to find two straws between the hit:
    vector<float>::iterator it;
    it = lower_bound(xLayer.begin(), xLayer.end(), xHit);
    // no smaller value  than val in vector
    if (it == xLayer.begin()){
        upper = *it; 
        if (debugBool && StrongDebugBool){ cout << " upper = " << upper << endl;}
         // hit distance is the dca from the L of the straw 0
        hitDistance = Tracker::DCA(upper, xHit);
        comparator = upper;
        LR=-1.0;
        if (debugBool && StrongDebugBool && (hitDistance > strawSpacing/2)){
            cout << "Hit at " << xHit << " The first straw is closest at " << upper <<  "; with DCA "<< hitDistance << endl;
        
        }
    }   
    // no bigger value than val in vector
    else if (it == xLayer.end()){
        lower = *(it-1); 
        if (debugBool && StrongDebugBool){ cout << "lower = " << lower << endl;}
        // hit distance is the dca from the R of the last straw
        hitDistance = Tracker::DCA(lower, xHit); // hit distance is the dca from the L of the straw 0
        comparator = lower;
        LR=1.0;
        if (debugBool && StrongDebugBool && (hitDistance > strawSpacing/2)){
            cout << "Hit at " << xHit << " The last straw is closest at " << lower <<  "; with DCA "<< hitDistance << endl;
            
        }
    }
    else {
        lower = *(it-1);
        upper = *it;
        if (debugBool && StrongDebugBool){ cout << "lower = " << lower << " upper = " << upper << endl;}

            float hit_distance_low = Tracker::DCA(lower, xHit);
            float hit_distance_up = Tracker::DCA(upper, xHit); 

            if (hit_distance_low < strawRadius && hit_distance_up < strawRadius){
                if (debugBool){cout << "Multiple straws in layer were hit!" << endl;}
                incMultipleHitsLayer();
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
            
            if (hit_distance_up==hit_distance_low){
                hitDistance = hit_distance_low; cout << "Ambiguity which straw registered hit" << endl;
                exit(0);
            } 
            if (debugBool && StrongDebugBool && 1==0){
                cout <<  "Hit at " << xHit << "Two straws closest to that point are " << lower << ", and " << upper <<  "; with DCAs "<< hit_distance_low <<  " and " << hit_distance_up << ", respectively." << endl; 
            }
    } // end of iterator to find straws between hits
   
    //Iterator to find straw ID [implemented separately, as now we look up once only]:
    int index;
    it = find(xLayer.begin(), xLayer.end(), comparator);
    if (it == xLayer.end()){   
        cout << "Error: straw x not found in vector!" << endl;
        exit(0);
      
    } else {
        index = std::distance(xLayer.begin(), it);
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
float Tracker::HitRecon(int x_det_ID, float x_det_dca, float LRSign, vector<float> xLayer){
    
    bool StrongDebugBool=false; //quick hack XXX for even more debug output

    float x_IdealStraw=xLayer[x_det_ID];

    float xRecon = x_IdealStraw + (LRSign * x_det_dca); //add dca for R  

    if (StrongDebugBool){ cout << "HitRecon= " << xRecon << "; x_IdealStraw= " << x_IdealStraw << " LRSign = " << LRSign  <<  " x_det_dca= " << x_det_dca << endl;}

    return xRecon;
}

// Function to return residuals to the fitted line (due to dca point scatter + resolution of the detector)
// @ Inputs: ideal points (x,z)

ResidualData Tracker::GetResiduals(vector<float> ReconPoints, ofstream& plot_fit){

    bool StrongDebugBool=false;
    ResidualData resData;
    float SumX=0;
    for (int i_size=0; i_size<ReconPoints.size(); i_size++){
        //RESIDUAL between a point (dca on a straw due to misalignment + ideal position) and detected position *
        SumX+=ReconPoints[i_size];
    if (StrongDebugBool){cout << "ReconPoints[i_size]= " << ReconPoints[i_size] << endl;}
    }
     if (StrongDebugBool){cout << "SumX= " << SumX <<endl;}

    float x_fit=SumX/ReconPoints.size(); 
   
    plot_fit << x_fit << " "  << x_fit << " " << beamStart << " " << beamStop << endl;
    
    for (int i_size=0; i_size<ReconPoints.size(); i_size++){
        //RESIDUAL between a point (dca on a straw due to misalignment + ideal position) and detected position *
        float residual_i = x_fit-ReconPoints[i_size];
        resData.residuals.push_back(residual_i);
        resData.x_fitted.push_back(x_fit);
        if (StrongDebugBool){cout << "residual_i= " << residual_i << " x_fit= " << x_fit << " ReconPoints= " << ReconPoints[i_size] << endl;}
    }
    
    return resData;
}
   


/**
   Simulate the passage of a randomly generated linear track through the detector, calculating properties of hits by the track on detector planes.
  
   @return MCData struct containing data about detector hits.
*/
MCData Tracker::MC_launch(float scatterError, ofstream& debug_calc, ofstream& debug_off, ofstream& debug_mc, ofstream& plot_fit, ofstream& plot_gen, ofstream& plot_hits_gen, bool debugBool) {
  
    // Set up new container for track data, with hit count set to zero
    MCData MC;
    MC.hit_count = 0; //count only counts through detector
    hitLayerCounter=-1; // start counting from 0, increment before hit-rejection
    
    //vector to store x coordinates of the track as seen from the ideal detector 
    vector<float> xReconPoints;  
    xReconPoints.clear();
    
   //Track parameters for rand-generated line MC [start and end positions outside of detectors]
    float x0 = 2.0 * Tracker::generate_uniform()-1.0; //uniform vertex
    //float x0 = Tracker::generate_uniform(); //XXX

    float xTrack=x0; //true track position always the same for parallel track
        
    //The main loop is for modules [they produce label of Global Parameters]: 
    // Then looping over layers [we know there will be at most 1 hit per layer] and views
    // within layers we will have loops over straws in that layer [to find out which straw was hit and dca]
    int z_counter=0; // distance in z is incremented from 0 to TotalLayerN 
    if (debugBool){cout << "Calculating hits/dca:" << endl;}
    //loops over modules, views, layers
    for(int i_module=0; i_module<moduleN; i_module++){    
    	for (int i_view=0; i_view<viewN; i_view++){
        	for (int i_layer=0; i_layer<layerN; i_layer++){ //per layer
	            hitLayerCounter++;
         
               //The registered hit position on the misaligned detector is smeared by its resolution 
               float xHit=xTrack+resolution*Tracker::generate_gaus();
               float residualTrack= xHit - xTrack;
               if (debugBool){plot_hits_gen << xHit << " " << distance[z_counter] << " ";}

	            //the dca between xHit and the xStraw [misaligned]
	            DCAData x_MisDetector= Tracker::DCAHit(mod_lyr_strawMisPosition[i_module][i_view][i_layer], distance[z_counter], xHit, debugBool); // position on the detector [from dca]  
	            float x_mis_dca =  x_MisDetector.dca;  //dca of the hit straw

	            //Rejection of hits due to geometry (i.e. missed hits)  
	            //No signal in a straw = no signal in the layer
	            // if (x_mis_dca > 0.5*strawSpacing){ //[between (0,0.25] 
	            //     MC.hit_bool.push_back(0);       
	            //     incRejectedHitsDCA();
	            //     if (debugBool){cout << "Hit Rejected: outside of straw layer with dca =" << x_mis_dca << endl;}
	            //     stringstream absolute_hit;
	            //     absolute_hit <<"#"; 
	            //     MC.absolute_straw_hit.push_back(absolute_hit.str().c_str());    
	            //     continue;
	            // } 

	            //Find the ID, and LR hit for a straw
	            int x_mis_ID =  x_MisDetector.strawID; // ID of the hit straw [to identify the correct straw in the Fit function]
	            float x_mis_LRSign = x_MisDetector.LRSign;
	            MC.strawID.push_back(x_mis_ID);
	            MC.LR.push_back(x_mis_LRSign);

                //DEBUG
                

                //Recording hit information
	            MC.hit_list.push_back(hitLayerCounter); //push back the absolute layer # into the hit list
	            MC.hit_bool.push_back(1);  // 1 = hit
	            stringstream absolute_hit;
	            //Making 5 chars consistent strings:
	            absolute_hit << x_mis_ID;
	            MC.absolute_straw_hit.push_back(absolute_hit.str().c_str());         

	            if(debugBool){cout << "DCA is= " << x_mis_dca << " for straw ID= " << x_mis_ID << " was hit from " << x_mis_LRSign << endl;}

	            //Reconstructing the hit as seen from the ideal detector
	            float xRecon = Tracker::HitRecon(x_mis_ID, x_mis_dca, x_mis_LRSign, mod_lyr_strawIdealPosition[i_module][i_view][i_layer]);
	            xReconPoints.push_back(xRecon); // vector to store x coordinates of the track as seen from the ideal detector
	            if (debugBool){plot_hits_gen << xRecon << " " << distance[z_counter] << endl;}
	            if(debugBool){cout << "Recon x= " << xRecon << endl;}
	            
	            //Module number [for labelling] - after (if) passing the rejection.
	            ostringstream oss;
	            oss << i_module + 1; // Millepede accepts only positive non-zero integers as labels
	            istringstream iss(oss.str());
	            int label_int;
	            iss >> label_int;
	            MC.i_hits.push_back(label_int); // vector of modules that were actually hit [after passing rejection test: for MP2 labelling]
	      
	            //Z-coordinate of hits [it was always known from geometry - no z-misalignment for now...]
	            MC.z_hits.push_back(distance[z_counter]);
	            MC.hit_sigmas.push_back(Tracker::instance()->getResolution()); 
	                
	            //Sanity Plots: Hits
	            MC.x_mis_dca.push_back(x_mis_dca);
	            MC.x_track.push_back(xTrack);
	            MC.x_recon.push_back(xRecon);
	            MC.residuals_gen.push_back(residualTrack);
	            //Sanity Plots: Tracks
	            if (MC.hit_count == 0){
	                           
	            } 

	            MC.hit_count++;

	            z_counter++;
	        }//end of Layers loop
        }// end of View loop
    }// end of looping over modules

    if (debugBool){cout << "Calculating residuals:" << endl;}
    //This happens once per MC function call [as we now accumulated x coordinates of "ideal" points for all hits
    // and need to do a simultaneous fit once - to return #hits worth of residuals]
    ResidualData res_Data = Tracker::GetResiduals(xReconPoints, plot_fit);
    MC.x_residuals = res_Data.residuals;
    MC.x_fitted=res_Data.x_fitted; //Sanity Plot 
    
    if(debugBool){
                        plot_gen << x0 << " " << x0 << " " << beamStart << " " << beamStop << " ";
                        for (int i=0; i<MC.absolute_straw_hit.size(); i++){
                        plot_gen     << MC.absolute_straw_hit[i] << " ";
                             }
                        plot_gen << endl;
                    }

    return MC; // Return data from simulated track
    
} // end of MC

//Geometry of detector arrangement (Ideal Geometry)
void Tracker::setGeometry(ofstream& debug_geom, bool debugBool){
   // float dZ=startingZDistanceStraw0; // to add distance in z
    
    float dZ=startingZDistanceStraw0; // the increment in z for consecutive layers
    int layer_n = 0;  // layer label 
    for (int i_module=0; i_module<moduleN; i_module++){
    	for (int i_view=0; i_view<viewN; i_view++){
    		for (int i_layer=0; i_layer<layerN; i_layer++){
    			distance.push_back(dZ); // vector will contain all z coordinates of layers 
    			layer.push_back(layer_n);  // layer label array [starting from 0th layer]
				//resolutionLayer.push_back(resolution); //resolution in each layer
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

    // Set mapping for U0...V1
    string tempMapping[4] = {"U0", "U1", "V0", "V1"};
   	for (int i_view=0; i_view<viewN; i_view++){
        UVmapping.push_back(vector<string> ()); //initialize the first index with a 2D vector
          for (int i_layer=0; i_layer<layerN; i_layer++){
             if (i_view==0){UVmapping[i_view].push_back(tempMapping[i_layer]);}
             if (i_view==1){UVmapping[i_view].push_back(tempMapping[i_layer+2]);}
         } //layers
    } // views

   
   // Print out:
    if (debugBool){
        int Zcounter=0;
        for (int i_module=0; i_module<moduleN; i_module++){
          for (int i_view=0; i_view<viewN; i_view++){
                for (int i_layer=0; i_layer<layerN; i_layer++){ //per module
                    cout << "IDEAL Mod " << i_module  << " " << UVmapping[i_view][i_layer] << " X : ";
                    for (int i_straw=0; i_straw<strawN; i_straw++){
                    cout << mod_lyr_strawIdealPosition[i_module][i_view][i_layer][i_straw]<<" ";
                    debug_geom << mod_lyr_strawIdealPosition[i_module][i_view][i_layer][i_straw] << " ";
                    } // straws 
                    cout << " | Z= " << distance[Zcounter] << " [cm]" << endl;  // XXX fix this 
                    debug_geom << distance[Zcounter] << endl;
                    Zcounter++;
                } // end of Layers
                
            }// end of View loop
        cout << endl;
        }//Modules 
    }


} // end of geom

// MC misalignment of detectors 
void Tracker::misalign(ofstream& debug_mis, bool debugBool){
     
    //Now misaligning detectors in x
    float misDispX; // effective misalignment 
    float sign = -1.0; // for up/down +/- direction
    for (int i_module=0; i_module<moduleN; i_module++){
    	//Fix first and last modules with no misalignment
        if (i_module==0 || i_module==moduleN-1){
            misDispX=0;
            sign=1.0;
        }
        if (i_module==1){
             misDispX=0.05;
        }
        else if(i_module==2){
            misDispX=0.1;
        }

        else if(i_module==3){
            misDispX=-0.03;
        }

        else if(i_module==4){
            misDispX=-0.045;
        }
        
        float dX = startingXDistanceStraw0+(misDispX*sign); // starting on the x-axis (z, 0+disp)
    	    	sdevX.push_back(misDispX*sign);  // vector of misalignment 
        mod_lyr_strawMisPosition.push_back(vector<vector<vector<float> > >()); //initialize the first index with a 2D vector
         for (int i_view=0; i_view<viewN; i_view++){
         	 mod_lyr_strawMisPosition[i_module].push_back(vector<vector<float> >()); //initialize the first index with a 2D vector
                for (int i_layer=0; i_layer<layerN; i_layer++){ 
                    mod_lyr_strawMisPosition[i_module][i_view].push_back(vector<float> ()); //initialize the first index with a 2D vector
                    for (int i_straw=0; i_straw<strawN; i_straw++){
                        mod_lyr_strawMisPosition[i_module][i_view][i_layer].push_back(dX); 
                        dX=strawSpacing+dX; //while we are in the same layer: increment straw spacing in x
                    } //end of Straws loop
            if (i_view==0){ dX=layerDisplacement+startingXDistanceStraw0+(misDispX*sign); } //set displacement in x for the next layer in the view
            if (i_view==1){ dX=startingXDistanceStraw0+(misDispX*sign); } //set displacement in x for the next layer in the view
            }//end of Layers loop
        }// end of View loop 
        sign=sign;
    }//Modules  


    // Print out:
    if (debugBool){
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

}//end of misalign


/**
   Write a constraint file to the supplied file-stream.

    @param constraint_file Reference to ofstream to write constraint file to. 
 */
// XXX constraints are ignored with HIP method 
void Tracker::write_constraint_file(ofstream& constraint_file, ofstream& debug_con, bool debugBool) {
    // constraints: fix centre straws in first/last layer
    // Check constraints file is open, then write. 
    if (constraint_file.is_open()) {
       
        //Evaluation of constraints
        float one = 1.0;
        
        for (int i_module = 0; i_module < moduleN; i_module++){ 
            constraint_file << "Constraint 0.0" << endl;
                int labelt=i_module+1; // Millepede doesn't like 0 as a label...
                //Fixing module 0 and the last module
                if (i_module==0 || i_module==moduleN-1){
                    constraint_file << labelt << " " << fixed << setprecision(5) << one<< endl;
                }
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