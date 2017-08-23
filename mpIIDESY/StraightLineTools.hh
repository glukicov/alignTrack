//Helper functions for geometrical mathematical operations involving straight lines
//Tom Stuttard (1st March 2016)

#ifndef StraightLineTools_hh
#define StraightLineTools_hh

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "gm2geom/coordSystems/CoordSystem3Vector.hh"

namespace gm2geom {

   class StraightLineTools {

      public :

       //! default constructor
       StraightLineTools() {};

       //
       // Compare line to point
       //

       //Find point on line that is closest to point (e.g. the "point of closest approach")
       //The line is defined by 2 points
       //If user specified "finiteLine", function returns false if the closest approach point is not within the two line points
       //The calculated value is returned by reference as "closestPointOnLine"
       static bool pointToLineCA(const gm2geom::CoordSystem3Vector& point, 
                                 const gm2geom::CoordSystem3Vector& lineStartPoint, const gm2geom::CoordSystem3Vector& linePointB, 
                                 gm2geom::CoordSystem3Vector& closestPointOnLine,
                                 bool finiteLine=false); 

       //Closest approach distance for point to line
       //The line is defined by 2 points
       //If user specified "finiteLine", function returns false if the closest approach point is not within the two line points
       //The calculated value is returned by reference as "dca"
       static bool pointToLineDCA(const gm2geom::CoordSystem3Vector& point, 
                                  const gm2geom::CoordSystem3Vector& lineStartPoint, const gm2geom::CoordSystem3Vector& linePointB, 
                                  double& dca,
                                  bool finiteLine=false); 


       //
       // Compare line to another line
       //

       //Find points where two skew lines approach closest
       //The lines are each defined by 2 points
       //If user specified "finiteLine", function returns false if either of the closest approach points are not within the line end points
       //The calculated values is returned by reference as "closestPointOnFirstLine" and "closestPointOnSecondLine"
       //FIXME Handle parallel line case
       static bool lineToLineCA( const gm2geom::CoordSystem3Vector& firstLineStartPoint, const gm2geom::CoordSystem3Vector& firstLineEndPoint, 
                                 const gm2geom::CoordSystem3Vector& secondLineStartPoint, const gm2geom::CoordSystem3Vector& secondLinePointB, 
                                 gm2geom::CoordSystem3Vector& closestPointOnFirstLine, gm2geom::CoordSystem3Vector& closestPointOnSecondLine,
                                 bool finiteLine=false);  

       //Closest approach distance for two lines
       //The lines are each defined by 2 points
       //If user specified "finiteLine", function returns false if either of the closest approach points are not within the line end points
       //The calculated value is returned by reference as "dca"
       static bool lineToLineDCA( const gm2geom::CoordSystem3Vector& firstLineStartPoint, const gm2geom::CoordSystem3Vector& firstLineEndPoint, 
                                  const gm2geom::CoordSystem3Vector& secondLineStartPoint, const gm2geom::CoordSystem3Vector& secondLinePointB, 
                                  double& dca,
                                  bool finiteLine=false);  


       //
       // Get point on a line from one coordinate
       //

       static bool pointOnLineFromX( const gm2geom::CoordSystem3Vector& lineStartPoint, const gm2geom::CoordSystem3Vector& linePointB, 
                                     double x,
                                     gm2geom::CoordSystem3Vector& outputPoint,
                                     bool finiteLine);

       static bool pointOnLineFromY( const gm2geom::CoordSystem3Vector& lineStartPoint, const gm2geom::CoordSystem3Vector& linePointB, 
                                     double y,
                                     gm2geom::CoordSystem3Vector& outputPoint,
                                     bool finiteLine);

       static bool pointOnLineFromZ( const gm2geom::CoordSystem3Vector& lineStartPoint, const gm2geom::CoordSystem3Vector& linePointB, 
                                     double z,
                                     gm2geom::CoordSystem3Vector& outputPoint,
                                     bool finiteLine);

   }; //end class

}

#endif
