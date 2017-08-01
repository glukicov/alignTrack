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

// ResidualData Tracker::GetResiduals(vector<float> ReconPoints,  vector<float> z_distance, ofstream& plot_fit, bool debugBool){

//     bool StrongDebugBool=true;
//      ResidualData resData;
//     //from: http://www.bragitoff.com/2015/09/c-program-to-linear-fit-the-data-using-least-squares-method/
//     int i,j,k,n;
//     //the no. of data pairs to be used
//     n=ReconPoints.size(); // # points to fit = number of layers that saw a hit 

//     double slope,intercept;
        
//     double xsum=0,x2sum=0,ysum=0,xysum=0;                //variables for sums/sigma of xi,yi,xi^2,xiyi etc
//     for (i=0;i<n;i++)
//     {
//         xsum=xsum+z_distance[i];                        //calculate sigma(xi)
//         ysum=ysum+ReconPoints[i];                        //calculate sigma(yi)
//         x2sum=x2sum+pow(z_distance[i],2);                //calculate sigma(x^2i)
//         xysum=xysum+z_distance[i]*ReconPoints[i];                    //calculate sigma(xi*yi)
//     }
//     slope=(n*xysum-xsum*ysum)/(n*x2sum-xsum*xsum);            //calculate slope
//     intercept=(x2sum*ysum-xsum*xysum)/(x2sum*n-xsum*xsum);            //calculate intercept
                             
//     vector<double> x_fit;  //an array to store the new fitted values of y  
//     for (i=0;i<n;i++)
//         x_fit.push_back(slope*z_distance[i]+intercept);                    //to calculate y(fitted) at given x points
    
//     if (debugBool && StrongDebugBool){
//         cout<<"# "<<setw(5)<<"z"<<setw(19)<<"x(observed)"<<setw(19)<<"x(fitted)"<<endl;
//         cout<<"-----------------------------------------------------------------\n";
//         for (i=0;i<n;i++)
//             cout<<i+1<<". "<<setw(8)<<z_distance[i]<<setw(15)<<ReconPoints[i]<<setw(18)<<x_fit[i]<<endl;//print a table of x,y(obs.) and y(fit.)    
//         cout<<"\nThe linear fit line is of the form:\nx="<<slope<<"z + "<<intercept<<endl << endl;        //print the best fit line
//     }

        
//     for (int i_size=0; i_size<n; i_size++){
//         //RESIDUAL between a point (dca on a straw due to misalignment + ideal position) and detected position
//         float residual_i = x_fit[i_size]-ReconPoints[i_size];
//         resData.residuals.push_back(residual_i);
//         resData.x_fitted.push_back(x_fit[i_size]);
//     }
//     resData.intercept_recon = intercept;
//     resData.slope_recon = slope;

//     if (debugBool){ plot_fit <<  slope*beamStart+intercept << " "  << slope*beamStop+intercept  <<   " " <<  beamStart  << " " << beamStop << endl; }
//     return resData;
// }
   
ResidualData Tracker::GetResiduals(vector<float> ReconPoints,  vector<float> z_distance, int dataSize, ofstream& plot_fit, bool debugBool){

    ResidualData resData;

    bool StrongDebugBool = false; 

    //from: http://codesam.blogspot.com/2011/06/least-square-linear-regression-of-data.html
    double SUMx = 0;     //sum of x values
    double SUMy = 0;     //sum of y values
    double SUMxy = 0;    //sum of x * y
    double SUMxx = 0;    //sum of x^2
    double SUMres = 0;   //sum of squared residue
    double res = 0;      //residue squared
    double slope = 0;    //slope of regression line
    double y_intercept = 0; //y intercept of regression line
    double SUM_Yres = 0; //sum of squared of the discrepancies
    double AVGy = 0;     //mean of y
    double AVGx = 0;     //mean of x
    double Yres = 0;     //squared of the discrepancies
    double Rsqr = 0;     //coefficient of determination

    //calculate various sums 
    for (int i = 0; i < dataSize; i++){
        //sum of x
        SUMx = SUMx + z_distance[i];
        //sum of y
        SUMy = SUMy + ReconPoints[i];
        //sum of squared x*y
        SUMxy = SUMxy +z_distance[i] * ReconPoints[i];
        //sum of squared x
        SUMxx = SUMxx + z_distance[i] * z_distance[i];
    }

    //calculate the means of x and y
    AVGy = SUMy / dataSize;
    AVGx = SUMx / dataSize;

    //slope or a1
    slope = (dataSize * SUMxy - SUMx * SUMy) / (dataSize * SUMxx - SUMx*SUMx);

    //y itercept or a0
    y_intercept = AVGy - slope * AVGx;

    if (StrongDebugBool){
        printf("x mean(AVGx) = %0.5E\n", AVGx);
        printf("y mean(AVGy) = %0.5E\n", AVGy);
        printf ("\n");
        printf ("The linear equation that best fits the given data:\n");
        printf ("       y = %2.8lfx + %2.8f\n", slope, y_intercept);
        printf ("------------------------------------------------------------\n");
        printf ("   Original (x,y)   (y_i - y_avg)^2     (y_i - a_o - a_1*x_i)^2\n");
        printf ("------------------------------------------------------------\n");
    }

    //calculate squared residues, their sum etc.
    for (int i = 0; i < dataSize; i++) {
        //current (y_i - a0 - a1 * x_i)^2
        Yres = pow((ReconPoints[i] - y_intercept - (slope * z_distance[i])), 2);
        
        float xFit = slope*z_distance[i]+y_intercept; 
        resData.x_fitted.push_back(xFit);
        float xRes = xFit - ReconPoints[i];
        resData.residuals.push_back(xRes); 
        
        //sum of (y_i - a0 - a1 * x_i)^2
        SUM_Yres += Yres;

        //current residue squared (y_i - AVGy)^2
        res = pow(ReconPoints[i] - AVGy, 2);

        //sum of squared residues
        SUMres += res;

        if (StrongDebugBool){
            printf ("   (%0.2f %0.2f)      %0.5E         %0.5E\n", 
            z_distance[i], ReconPoints[i], res, Yres);
        }
    }

    //calculate r^2 coefficient of determination
    Rsqr = (SUMres - SUM_Yres) / SUMres;

    if (StrongDebugBool){
        printf("--------------------------------------------------\n");
        printf("Sum of (y_i - y_avg)^2 = %0.5E\t\n", SUMres);
        printf("Sum of (y_i - a_o - a_1*x_i)^2 = %0.5E\t\n", SUM_Yres);
        printf("Standard deviation(St) = %0.5E\n", sqrt(SUMres / (dataSize - 1)));
        printf("Standard error of the estimate(Sr) = %0.5E\t\n", sqrt(SUM_Yres / (dataSize-2)));
        printf("Coefficient of determination(r^2) = %0.5E\t\n", (SUMres - SUM_Yres)/SUMres);
        printf("Correlation coefficient(r) = %0.5E\t\n", sqrt(Rsqr));
    }
    
    resData.slope_recon=slope;
    resData.intercept_recon=y_intercept;

    return resData;
}






/**
   Simulate the passage of a randomly generated linear track through the detector, calculating properties of hits by the track on detector planes.
  
   @return MCData struct containing data about detector hits.
*/
MCData Tracker::MC_launch(float scatterError, ofstream& debug_calc, ofstream& debug_off, ofstream& debug_mc, ofstream& plot_fit, ofstream& plot_gen, ofstream& plot_hits_gen, ofstream& plot_hits_fit, bool debugBool) {
    
    bool StrongDebugBool = false; // HACK
    // Set up new container for track data, with hit count set to zero
    MCData MC;
    MC.hit_count = 0; //count only counts through detector
    hitLayerCounter=-1; // start counting from 0, increment before hit-rejection
    
    //vector to store x coordinates of the track as seen from the ideal detector 
    vector<float> xReconPoints;  
    xReconPoints.clear();
    
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

    //cout << "xSlope= " << xSlope << " x0= " << x0 << " beamPositionLength= " << beamPositionLength << " beamStop= " << beamStop << " x1= " << x1 << endl;  

    float xTrack; //true track position x=ym+c 

    MC.x0 = x0;
    MC.x1 = x1;
    MC.slope_truth = xSlope;
    MC.intercept_truth = xIntercept;
        
    //The main loop is for modules [they produce label of Global Parameters]: 
    // Then looping over layers [we know there will be at most 1 hit per layer] and views
    // within layers we will have loops over straws in that layer [to find out which straw was hit and dca]
    int z_counter=0; // distance in z is incremented from 0 to TotalLayerN 
    if (debugBool && StrongDebugBool){cout << "Calculating hits/dca:" << endl;}
    //loops over modules, views, layers
    for(int i_module=0; i_module<moduleN; i_module++){    
    	for (int i_view=0; i_view<viewN; i_view++){
        	for (int i_layer=0; i_layer<layerN; i_layer++){ //per layer
	            hitLayerCounter++;
         
               //The registered hit position on the misaligned detector is smeared by its resolution 
               xTrack = xSlope*distance[z_counter]+xIntercept;  // true track position [from line x=ym+c]
               MC.x_track_true.push_back(xTrack);  // True (gen.) track position

               float xHit=xTrack+resolution*Tracker::generate_gaus();
               float residualTrack= xHit - xTrack;
               if (debugBool){plot_hits_gen << xHit << " " << distance[z_counter] << endl;}

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

	            if(debugBool && StrongDebugBool){cout << "DCA is= " << x_mis_dca << " for straw ID= " << x_mis_ID << " was hit from " << x_mis_LRSign << endl;}

	            //Reconstructing the hit as seen from the ideal detector
	            float xRecon = Tracker::HitRecon(x_mis_ID, x_mis_dca, x_mis_LRSign, mod_lyr_strawIdealPosition[i_module][i_view][i_layer]);
	            xReconPoints.push_back(xRecon); // vector to store x coordinates of the track as seen from the ideal detector
	            if (debugBool){plot_hits_fit << xRecon << " " << distance[z_counter] << endl;}
	            if(debugBool && StrongDebugBool){cout << "Recon x= " << xRecon << endl;}
	            
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
	            MC.x_mis_dca.push_back(x_mis_dca); // DCA
	            MC.x_hit_recon.push_back(xRecon);   // Reconstructed hit position
	            MC.residuals_gen.push_back(residualTrack); // Residual between the XHit and the generated track
                MC.x_hit_true.push_back(xHit);  //True (smeared) hit position

                MC.Module_i.push_back(i_module);
                MC.View_i.push_back(i_view);
                MC.Layer_i.push_back(i_layer);
                MC.Straw_i.push_back(x_mis_ID);

	            //Sanity Plots: Tracks
	            if (MC.hit_count == 0){
	                           
	            } 

	            MC.hit_count++;

	            z_counter++;
	        }//end of Layers loop
        }// end of View loop
    }// end of looping over modules

    if (debugBool && StrongDebugBool){cout << "Calculating residuals:" << endl;}
    //This happens once per MC function call [as we now accumulated x coordinates of "ideal" points for all hits
    // and need to do a simultaneous fit once - to return #hits worth of residuals]
    //ResidualData res_Data = Tracker::GetResiduals(xReconPoints, distance, plot_fit, debugBool);
    
    ResidualData res_Data = Tracker::GetResiduals(xReconPoints, distance, MC.hit_count, plot_fit, debugBool);
    MC.x_residuals = res_Data.residuals;
    MC.x_track_recon=res_Data.x_fitted; //Sanity Plot: fitted (reconstructed) x of the track
    MC.slope_recon = res_Data.slope_recon;
    MC.intercept_recon = res_Data.intercept_recon;
    
    if(debugBool){
                        plot_gen << x0 << " " << x1 << " " << beamStart << " " << beamStop << " ";
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
    for (int i_module=0; i_module<moduleN; i_module++){
    	//Fix first and last modules with no misalignment
        if (i_module==0 || i_module==moduleN-1){
            misDispX=0.0;    
        }
        // the rest of the modules get misalignment parameter from the vector 
        else{
            misDispX=dispX[i_module];
        }       
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

float sum_of_elems;
for(std::vector<float>::iterator it = sdevX.begin(); it != sdevX.end(); ++it)
    sum_of_elems += *it;
overallMis=sum_of_elems/moduleN;

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
            //if (i_module==0){
            if (i_module==0 || i_module==moduleN-1){
            	constraint_file << "Constraint 0.0" << endl;
                	int labelt=i_module+1; // Millepede doesn't like 0 as a label...
                	//Fixing module 0 and the last module
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