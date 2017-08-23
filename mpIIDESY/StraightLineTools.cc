//! header file
#include "gm2geom/common/StraightLineTools.hh"
 
bool gm2geom::StraightLineTools::pointToLineCA(
  const gm2geom::CoordSystem3Vector& point, 
  const gm2geom::CoordSystem3Vector& lineStartPoint, const gm2geom::CoordSystem3Vector& lineEndPoint, 
  gm2geom::CoordSystem3Vector& closestPointOnLine,
  bool finiteLine) 
{

  /*

  Some notes of formalism:

    Using 3D parametric line equations using two points r0=(x0,y0,z0),r1=(x1,y1,z1): 

      - A point on this line, r(t), is give by r = { x,y,z } = { x0+(x1-x0)t , y0+(y1-y0)t , z0+(z1-z0)t } 
                                               = { x0+m_x*t , y0+m_y*t , z0+m_z*t } = r0 + m*t
      - This is basically just "point 0" + "fraction t along line" * "gradient"
      - At t=0, have point r0, at t=1, have point r1, so for finite line range of t is [0,1]
      - Note that r(t) is a position relative to the origin, NOT a vector along the line direction

  */

  //Line is r(t), defined by end points r0,r1
  //Point is p = { x_p , y_p , z_p }

  //Distance from point to line is minimized when t = (r0-p).(r1-r0) / |r1-r0|^2
  //From http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  double tMin = - (lineStartPoint-point).dot(lineEndPoint-lineStartPoint) / pow((lineEndPoint-lineStartPoint).mag(),2);

  //If line is finite, check that this closest point lies between the two points specified
  //Do this by checking that the minimized t lies in range [0,1], as r(t=0)=r0 and r(t=1)=r1
  if(finiteLine) {
    if( (tMin<0.) || (tMin>1.) ) return false;
  }

  //Put tMin back into equation to get point along line that is closest to p
  closestPointOnLine.set( lineStartPoint.x() + (lineEndPoint.x()-lineStartPoint.x())*tMin ,
                          lineStartPoint.y() + (lineEndPoint.y()-lineStartPoint.y())*tMin ,
                          lineStartPoint.z() + (lineEndPoint.z()-lineStartPoint.z())*tMin ,
                          point.coordSystemName );

  return true;
}


bool gm2geom::StraightLineTools::pointToLineDCA(
  const gm2geom::CoordSystem3Vector& point, 
  const gm2geom::CoordSystem3Vector& lineStartPoint, const gm2geom::CoordSystem3Vector& lineEndPoint, 
  double& dca,
  bool finiteLine) 
{

  //Get closest approach point on line
  gm2geom::CoordSystem3Vector closestPointOnLine;
  bool success = pointToLineCA(point,lineStartPoint,lineEndPoint,closestPointOnLine,finiteLine);

  //Return magnitude of vector from point to point on line as DCA
  dca = success ? (closestPointOnLine-point).mag() : -1.;
  return success;

}


bool gm2geom::StraightLineTools::lineToLineCA(
  const gm2geom::CoordSystem3Vector& firstLineStartPoint, const gm2geom::CoordSystem3Vector& firstLineEndPoint, 
  const gm2geom::CoordSystem3Vector& secondLineStartPoint, const gm2geom::CoordSystem3Vector& secondLineEndPoint, 
  gm2geom::CoordSystem3Vector& closestPointOnFirstLine, gm2geom::CoordSystem3Vector& closestPointOnSecondLine,
  bool finiteLine)
{

  //Using http://mipp-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=415&filename=dca_vert.pdf
  //See section "Distance and Points of Closest Approach for Straight Tracks"
  //TODO Also code up the curved tracks version and make available here

  //Have two lines here, P and Q, defined parametrically:
  //  p(t) = p0 + tu (p(t) is point on P, p0 = p(t=0) is known point on P, u is unit vector along P (from p(t=0) to p(t=1)), t is parameteric scalar
  //  q(t) = q0 + sv (q(t) is point on Q, q0 = q(s=0) is known point on Q, v is unit vector along Q (from q(s=0) to q(s=1)), s is parameteric scalar
  //t and s are measured in units of "unit vectors" (note difference from line-point equation above, where t is expressed in range [0,1])
  //Also need difference between the two known points: d = p0 - q0

  const gm2geom::CoordSystem3Vector & p0 = firstLineStartPoint;
  const gm2geom::CoordSystem3Vector u = ( firstLineEndPoint - firstLineStartPoint ).unit();

  const gm2geom::CoordSystem3Vector & q0 = secondLineStartPoint;
  const gm2geom::CoordSystem3Vector v = ( secondLineEndPoint - secondLineStartPoint ).unit();

  const gm2geom::CoordSystem3Vector d = p0 - q0;

  //Closest approach is when magnitude of d is minimized w.r.t. t and s
  //These are given by:
  //  t_min = [-d.u + (d.v)(u.v)] / [1 - (u.v)^2]
  //  s_min = [ d.v - (d.u)(u.v)] / [1 - (u.v)^2]

  //Pre-compute repeated terms for speed
  double u_dot_v = u.dot(v);
  double d_dot_u = d.dot(u);
  double d_dot_v = d.dot(v);
  double denominator = 1 - (u_dot_v*u_dot_v);

  //Calculate t and s values at closest approach 
  double tMin = ( -d_dot_u + (d_dot_v*u_dot_v) ) / denominator;
  double sMin = (  d_dot_v - (d_dot_u*u_dot_v) ) / denominator;

  //If lines are finite, check the closest approach points are within range
  //This means that length of line is <= t or s (since t or s are in units of [unit vectors]) 
  // or if get -ve t or s
  if(finiteLine) {
    if( tMin > (firstLineEndPoint-firstLineStartPoint).mag() ) return false;
    if( tMin < 0. ) return false;
    if( sMin > (secondLineEndPoint-secondLineStartPoint).mag() ) return false;
    if( sMin < 0. ) return false;
  }

  //Get the closest approach point on each line using the minimised parameteric scalars
  closestPointOnFirstLine = firstLineStartPoint + ( tMin * u );
  closestPointOnSecondLine = secondLineStartPoint + ( sMin * v );

  return true;

}


bool gm2geom::StraightLineTools::lineToLineDCA(
  const gm2geom::CoordSystem3Vector& firstLineStartPoint, const gm2geom::CoordSystem3Vector& firstLineEndPoint, 
  const gm2geom::CoordSystem3Vector& secondLineStartPoint, const gm2geom::CoordSystem3Vector& secondLineEndPoint, 
  double& dca,
  bool finiteLine)
{

   //Get closest approach points on the lines
   gm2geom::CoordSystem3Vector closestPointOnFirstLine, closestPointOnSecondLine;
   bool success = gm2geom::StraightLineTools::lineToLineCA( firstLineStartPoint, firstLineEndPoint,
                                                            secondLineStartPoint, secondLineEndPoint,
                                                            closestPointOnFirstLine, closestPointOnSecondLine,
                                                            finiteLine);

  //Return magnitude of vector between the two points as DCA
  dca = success ? (closestPointOnSecondLine-closestPointOnFirstLine).mag() : -1.;
  return success;

}  


bool gm2geom::StraightLineTools::pointOnLineFromX(
  const gm2geom::CoordSystem3Vector& lineStartPoint, const gm2geom::CoordSystem3Vector& lineEndPoint, 
  double x,
  gm2geom::CoordSystem3Vector& outputPoint,
  bool finiteLine)
{

  //Get t for given x coord (line eqn: r = r0 + mt)
  gm2geom::CoordSystem3Vector gradient = lineEndPoint - lineStartPoint;  
  double t = ( x - lineStartPoint.x() ) / gradient.x();

  //If line is finite, check that this closest point lies between the two points specified
  //Do this by checking that the minimized t lies in range [0,1], as r(t=0)=r0 and r(t=1)=r1
  if(finiteLine) {
    if( (t<0.) || (t>1.) ) return false;
  }

  //Use t to get all point coords
  outputPoint.set( x , lineStartPoint.y()+(gradient.y()*t) , lineStartPoint.z()+(gradient.z()*t) , lineStartPoint.coordSystemName );

  return true;

}


bool gm2geom::StraightLineTools::pointOnLineFromY(
  const gm2geom::CoordSystem3Vector& lineStartPoint, const gm2geom::CoordSystem3Vector& lineEndPoint, 
  double y,
  gm2geom::CoordSystem3Vector& outputPoint,
  bool finiteLine)
{

  //Get t for given y coord (line eqn: r = r0 + mt)
  gm2geom::CoordSystem3Vector gradient = lineEndPoint - lineStartPoint;  
  double t = ( y - lineStartPoint.y() ) / gradient.y();

  //If line is finite, check that this closest point lies between the two points specified
  //Do this by checking that the minimized t lies in range [0,1], as r(t=0)=r0 and r(t=1)=r1
  if(finiteLine) {
    if( (t<0.) || (t>1.) ) return false;
  }

  //Use t to get all point coords
  outputPoint.set( lineStartPoint.x()+(gradient.x()*t) , y , lineStartPoint.z()+(gradient.z()*t) , lineStartPoint.coordSystemName );

  return true;

}


bool gm2geom::StraightLineTools::pointOnLineFromZ(
  const gm2geom::CoordSystem3Vector& lineStartPoint, const gm2geom::CoordSystem3Vector& lineEndPoint, 
  double z,
  gm2geom::CoordSystem3Vector& outputPoint,
  bool finiteLine)
{

  //Get t for given x coord (line eqn: r = r0 + mt)
  gm2geom::CoordSystem3Vector gradient = lineEndPoint - lineStartPoint;  
  double t = ( z - lineStartPoint.z() ) / gradient.z();

  //If line is finite, check that this closest point lies between the two points specified
  //Do this by checking that the minimized t lies in range [0,1], as r(t=0)=r0 and r(t=1)=r1
  if(finiteLine) {
    if( (t<0.) || (t>1.) ) return false;
  }

  //Use t to get all point coords
  outputPoint.set( lineStartPoint.x()+(gradient.x()*t) , lineStartPoint.y()+(gradient.y()*t) , z , lineStartPoint.coordSystemName );

  return true;

}

