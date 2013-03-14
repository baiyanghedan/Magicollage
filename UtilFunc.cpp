#include "UtilFunc.h"

//
double SquaredDistance(const Vector2D vec1, const Vector2D vec2)
{
	return ((vec1.x-vec2.x)*(vec1.x-vec2.x) + (vec1.y-vec2.y)*(vec1.y-vec2.y));
}
//
double SqrtDistance(const Vector2D vec1, const Vector2D vec2)
{
	double res = sqrt((vec1.x-vec2.x)*(vec1.x-vec2.x) + (vec1.y-vec2.y)*(vec1.y-vec2.y));
	return res;
}
//
bool isInCircle(double x, double y)
{
	double sqrdist = (x-center_x)*(x-center_x) + (y-center_y)*(y-center_y);
	if(sqrdist < radius_r*radius_r)
		return true;
	return false;
}
//
bool isInBox(double x, double y)
{
	if(x>bound_x_min && x<bound_x_max && y>bound_y_min && y<bound_y_max)
		return true;
	return false;
}
//
bool isCircleInCircle(double x, double y, double r, double cx, double cy, double cr, double& delta)
{
	double dist1 = (x-cx)*(x-cx) + (y-cy)*(y-cy);
	double dist2 = (cr-r)*(cr-r);
	if (dist1 < dist2)
		return true;
	delta = sqrt(dist1)+r-cr;
	return false;
}
//
bool is2CirclesOverlapped(double x1, double y1, double r1, double x2, double y2, double r2)
{
	double dist1 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
	double dist2 = (r1+r2)*(r1+r2);

	//The tolerance is very critical here!!!
	//Observation: 1e-10 plays better than 1e-12 in DS1
	//The smaller the tolerance is, the easier it is classified as overlapped, and be dropped off	
	if ((dist1 + TOLER) < dist2)
		return true;
	return false;
}
//
bool IsSamePoint(Point p1, Point p2)
{
	if ( (fabs(p1.x() - p2.x()) < TOLER ) && ( fabs(p2.y() - p1.y()) < TOLER ) )
		return true;
	return false;
}

//
bool  IsSameWeightedPoint(WeightedPoint wp1, WeightedPoint wp2)
{
	if (IsSamePoint(wp1.point(), wp2.point()) && (fabs(wp1.weight() - wp2.weight()) < TOLER ) )
		return true;
	return false;
}

// ------
// |    |
// |    |
// ------
// 0: y=bound_y_min; 1: x=bound_x_max; 2: y=bound_y_max; 3 x=bound_x_min
double intersect_rect_bound(double x, double y, double dx, double dy, int& index)
{
	double t1=10e10,t2=10e10,t3=10e10,t4=10e10;
	if(dx!=0)
	{
		t1 = (bound_x_max-x)/dx;
		if(t1<0)
			t1 = 10e10;
		t2 = (bound_x_min-x)/dx;
		if(t2<0)
			t2 = 10e10;
	}
	if(dy!=0)
	{
		t3 = (bound_y_max-y)/dy;
		if(t3<0)
			t3 = 10e10;
		t4 = (bound_y_min-y)/dy;
		if(t4<0)
			t4 = 10e10;
	}

	double t;
	index = t1<t2?1:3;		t = t1<t2?t1:t2;	
	index = t<t3?index:2;	t = t<t3?t:t3;
	index = t<t4?index:0;	t = t<t4?t:t4;
	return t;
}

//
bool intersect_circle_bound(double x, double y, double dx, double dy, double& t)
{
	double dxdy = dx*dx + dy*dy;
	double deter = 2 * x * dx * y * dy - 2 * x * dx * dy * center_y - 2 * y * dy * dx * center_x + 2 * dy * center_y * dx * center_x - dx * dx * center_y * center_y + 2 * dx * dx * y * center_y - dx * dx * y * y  - dy * dy * x * x - dy * dy * center_x * center_x + 2 * dy * dy * x * center_x + dxdy * radius_r * radius_r;
	if (deter < 0)
	{//no intersection
		return false;
	}
	double tmp = -x*dx - y*dy + dy*center_y + dx*center_x;
	double sqrtdeter = sqrt(deter);
	double t1 = (tmp + sqrtdeter) /dxdy;
	double t2 = (tmp - sqrtdeter) /dxdy;
	t = t1>t2?t1:t2;
	return true;
}
//
bool intersect_circle_bound(double x, double y, double dx, double dy, double& t1, double& t2)
{
	double dxdy = dx*dx + dy*dy;
	if (dxdy < TOLER)
	{
		return false;
	}
	if ((abs(x)>1e15) ||(abs(y)>1e15) ||(abs(dx)>1e15)||(abs(dy)>1e15)){
		double denom = sqrt(dxdy);
		dx = dx/denom;
		dy = dy/denom;
		double tmp = 1-(dx*y-dy*x)*(dx*y-dy*x);
		if (tmp<0)
		{
			return false;
		}
		double t = -(x*dx + y*dy+sqrt(tmp));
		x = x + t*dx;
		y = y + t*dy;
		dxdy = 1;
	}
	double deter = 2 * x * dx * y * dy - 2 * x * dx * dy * center_y - 2 * y * dy * dx * center_x + 2 * dy * center_y * dx * center_x - dx * dx * center_y * center_y + 2 * dx * dx * y * center_y - dx * dx * y * y  - dy * dy * x * x - dy * dy * center_x * center_x + 2 * dy * dy * x * center_x + dxdy * radius_r * radius_r;
	if (deter < 0)
	{//no intersection
		return false;
	}
	double tmp = -x*dx - y*dy + dy*center_y + dx*center_x;
	double sqrtdeter = sqrt(deter);
	t1 = (tmp + sqrtdeter) /dxdy;
	t2 = (tmp - sqrtdeter) /dxdy;
	assert(t1>t2);
	return true;
}
//
void GetBoundary(int idx, Vector2D& E, Vector2D& F)
{
	switch(idx) {
		case 0:
			E = Vector2D(bound_x_min, bound_y_min);
			F = Vector2D(bound_x_max, bound_y_min);
			break;
		case 1:
			E = Vector2D(bound_x_max, bound_y_min);
			F = Vector2D(bound_x_max, bound_y_max);
			break;
		case 2:
			E = Vector2D(bound_x_max, bound_y_max);
			F = Vector2D(bound_x_min, bound_y_max);
			break;
		case 3:
			E = Vector2D(bound_x_min, bound_y_max);
			F = Vector2D(bound_x_min, bound_y_min);
			break;
		default:
			break;
	}	
}
//
int LinesegLineInter(Point2D linesp, Vector2D linedir, Point2D spt, Point2D ept, Point2D& pnt2d, double tol)
{
	double delta, deltat;
	Vector2D vn1 = ept -spt;
	Vector2D vn2 = linedir;
	vn2.SetUnit();
	Vector2D veca, vecb;
	veca = spt - linesp;
	vecb = ept - linesp;

	double d1, d2;
	d1 = veca.x*vn2.y - veca.y*vn2.x;
	d2 = vecb.x*vn2.y - vecb.y*vn2.x;

	double len1, len2;
	len1 = Length(veca);
	len2 = Length(vecb);

	double tol1, tol2;
	tol1 = tol*len1;
	tol2 = tol*len2;

	if(len1<tol)
	{
		if(fabs(d2)>tol2)
		{
			pnt2d = spt;
			return 1;   
		}
		else   
		{
			return 2;
		}		
	}

	if(len2<tol)
	{
		if(fabs(d1)>tol1)
		{
			pnt2d = ept;
			return 1;    
		}
		else    
		{
			return 2;
		}		
	}


	if(d1<-tol1) 
	{
		if(d2<-tol2)
			return 0;   
		else if(d2<tol2)
		{
			pnt2d = ept;
			return 1;
		}
		else 
		{
			delta = vn1.x*vn2.y - vn1.y*vn2.x;
			deltat = vn1.x*veca.y - vn1.y*veca.x;
			pnt2d = linesp + vn2*deltat/delta;
			return 1;
		}
	}

	else if(d1<tol1)
	{
		if(fabs(d2)>tol2)
		{
			pnt2d = spt;
			return 1;
		}
		else 
		{
			return 2;
		}
	}
	else  
	{
		if(d2<-tol2)  
		{
			delta = vn1.x*vn2.y - vn1.y*vn2.x;
			deltat = vn1.x*veca.y - vn1.y*veca.x;
			pnt2d = linesp + vn2*deltat/delta;
			return 1;
		}
		else if(d2<tol2) 
		{
			pnt2d = ept;
			return 1;
		}
		else  
		{
			return 0;
		}
	}
}

//Rotate the point (x,y,0) around the vector (u,v,0) by 45 degrees
//pt and newvec must be the unit vector
void RotatePointAroundVector(Vector2D pt, Vector2D vec, Vector3D& newvec)
{
	double ux = pt.x*vec.x; //ux#=u*x
	double uy = pt.y*vec.x; //uy#=u*y
	
	double vx = pt.x*vec.y; //vx#=v*x
	double vy = pt.y*vec.y; //vy#=v*y	
	
	double sa = 0.5*sqrtf(2); //sa#=Sin(a)
	double ca = -sa; //ca#=Cos(a)
	double x = vec.x*(ux+vy) + (pt.x*(vec.y*vec.y)-vec.x*vy)*ca;
	double y = vec.y*(ux+vy) + (pt.y*(vec.x*vec.x)-vec.y*ux)*ca;
	double z = (-vx+uy)*sa;
	newvec = Vector3D(x, y, z);

	newvec.SetUnit();
}

//
bool ThreePlanesInter(Plane pl1, Plane pl2, Plane pl3, Point3D& pt)
{
	Vector3D pl23 = pl2.normal^pl3.normal;
	if (fabs(Length(pl23)) < TOLER )
		return false;
	pt =  (pl1.d*(pl2.normal^pl3.normal) + pl2.d*(pl3.normal^pl1.normal)+  pl3.d*(pl1.normal^pl2.normal))/(pl1.normal*pl23);
	return true;
}

//
bool TwoPlanesConeInter(Plane pl1, Plane pl2, Point p1, Point p2, bool& b2res, Point3D& pt1, Point3D& pt2)
{	
	Vector3D pl12 = pl1.normal^pl2.normal;
	double tmp = pl2.normal.x*pl1.normal.y - pl1.normal.x*pl2.normal.y;
	double x = (pl2.d*pl1.normal.y - pl1.d*pl2.normal.y)/tmp;
	double y = (pl1.d*pl2.normal.x - pl2.d*pl1.normal.x)/tmp;
	double z = 0;

	double nx = pl12.x;
	double ny = pl12.y;
	double nz = pl12.z;

	double r3 = pow(radius_r, 3);
	double r4 = pow(radius_r, 4);
	double deter = -4 * nz * nz * radius_r * radius_r * y * center_y - 4 * nz * r3 * center_x * nx + 2 * center_x * nx * center_y * ny + 2 * center_x * nx * radius_r * nz + 2 * center_y * ny * radius_r * nz + 2 * nz * nz * radius_r * radius_r * center_x * center_x + 2 * nz * nz * radius_r * radius_r * y * y + 2 * nz * nz * radius_r * radius_r * center_y * center_y + 2 * nz * nz * radius_r * radius_r * x * x + 4 * nx * x * nz * r3 - 4 * nz * nz * radius_r * radius_r * x * center_x + 2 * nx * nx * r4 + 2 * ny * ny * r4 + 4 * ny * y * nz * r3 + 2 * ny * y * nx * x - nx * nx * center_y * center_y - nx * nx * y * y - ny * ny * center_x * center_x - ny * ny * x * x - nx * nx * radius_r * radius_r - nz * nz * center_y * center_y - nz * nz * y * y - nz * nz * center_x * center_x - nz * nz * x * x - ny * ny * radius_r * radius_r - 2 * ny * y * radius_r * nz - 2 * ny * y * center_x * nx + 2 * ny * ny * x * center_x + 2 * nx * nx * y * center_y + 2 * nz * nz * y * center_y - 2 * nx * x * radius_r * nz - 2 * nx * x * center_y * ny + 2 * nz * nz * x * center_x - 4 * nz * r3 * center_y * ny;

	if (deter < 0)
	{//no intersection
		return false;
	}
	double sqrt_deter = sqrt(deter);
	double a = nx * x - center_x * nx - center_y * ny + ny * y + 2 * nz * r3 - nz * radius_r;

	double denom = -nx * nx - ny * ny + 2 * nz * nz * radius_r * radius_r - nz * nz;

	double t1 = (a + sqrt_deter)/denom;
	double t2 = (a - sqrt_deter)/denom;

	Point2D pos1 = Point2D(x+t1*nx, y+t1*ny);
	Point2D pos2 = Point2D(x+t2*nx, y+t2*ny);
	
	bool bIncircle1 = isInCircle(pos1.x, pos1.y);
	bool bIncircle2 = isInCircle(pos2.x, pos2.y);

	bool bPositveZ1 = (t1*nz > TOLER) && (t1*nz < radius_r);
	bool bPositveZ2 = (t2*nz > TOLER) && (t2*nz < radius_r);

	if (!bIncircle1 && !bIncircle2)
	{
		return false;
	}else if (!bIncircle1)
	{
		if (bPositveZ2)
		{
			pt1 = Point3D(pos2.x, pos2.y, t2*nz );
		}else{
			return false;
		}
	}else if (!bIncircle2)
	{
		if (bPositveZ1)
		{
			pt1 = Point3D(pos1.x, pos1.y, t1*nz );
		}else{
			return false;
		}
	}else{// two pts are both inside the circle
		if (bPositveZ1 && bPositveZ2)
		{
			Vector2D vvt1 = pos1 - Point2D(center_x, center_y);
			Vector2D vvt2 = pos2 - Point2D(center_x, center_y);
			Vector2D vv1 = Vector2D(p1.x(), p1.y())- Vector2D(center_x, center_y);
			Vector2D vv2 = Vector2D(p2.x(), p2.y())- Vector2D(center_x, center_y);

			vvt1.SetUnit(); vvt2.SetUnit(); vv1.SetUnit(); vv2.SetUnit();
			double res1 = (vv1^vvt1) * (vvt1^vv2);
			double res2 = (vv1^vvt2) * (vvt2^vv2);

			if (res1 < TOLER && res2 < TOLER)
			{
				return false;
			}else if (res1 < TOLER)
			{
				pt1 = Point3D(pos2.x, pos2.y, t2*nz );
			}else if (res2 < TOLER)
			{
				pt1 = Point3D(pos1.x, pos1.y, t1*nz );
			}else{
				pt1 = Point3D(pos1.x, pos1.y, t1*nz );
				pt2 = Point3D(pos2.x, pos2.y, t2*nz );
				b2res = true;
			}		
		}else if (bPositveZ1)
		{
			pt1 = Point3D(pos1.x, pos1.y, t1*nz );
		}else if (bPositveZ2)
		{
			pt1 = Point3D(pos2.x, pos2.y, t2*nz );
		}else{
			return false;
		}
	}

	return true;
}

bool isInsideCone(Point3D pt)
{	
	double r = radius_r - pt.z;
	double dist = (pt.x - center_x)*(pt.x - center_x) + (pt.y - center_y)*(pt.y - center_y);
	if ((dist - r*r)>TOLER)
	{
		return false;
	}
	return true;
}

//the rect boundary is (x0, y0), (x1, y0), (x1, y1), (x0, y1)
//to find the circle (x,y),r inscribed with the rect boundary and circumscribed with (xp, yp) ,rp
bool tangentto2Circles_inscribe_rect(const double x0, const double y0, const double x1, const double y1, 
									 const double xp, const double yp, const double rp, const double r, 
									 double* x, double* y)
{
	bool found = false;
//case 1: x-x0=r
//case 2: x1-x=r
	for (int i=0; i<2; i++)
	{
		double xnew;
		if (i==0)	xnew = x0+r;
		else xnew = x1-r;

		double xxp = xnew*xnew -2*xnew*xp + xp*xp;
		double c = yp*yp+xxp-(r+rp)*(r+rp);
		double delta = 4*yp*yp-4*c;
		if (delta>0){
			double y_res1 = (2*yp+sqrt(delta))*0.5;
			double y_res2 = (2*yp-sqrt(delta))*0.5;
			if ((y_res1<(y1-r)) && (y_res1>(y0+r))){
				x[0] = xnew;
				y[0] = y_res1;
				found = true;
				goto stop;
			}else if((y_res2<(y1-r)) && (y_res2>(y0+r))){
				x[0] = xnew;
				y[0] = y_res2;
				found = true;
				goto stop;
			}
		}
	}
	
//case 3: y-y0=r
//case 4: y1-y=r	
	for (int i=0; i<2; i++)
	{
		double ynew;
		if (i==0)	ynew = y0+r;
		else ynew = y1-r;

		double yyp = ynew*ynew -2*ynew*yp + yp*yp;
		double c = xp*xp+yyp-(r+rp)*(r+rp);
		double delta = 4*xp*xp-4*c;
		if (delta>0){
			double x_res1 = (2*xp+sqrt(delta))*0.5;
			double x_res2 = (2*xp-sqrt(delta))*0.5;
			if ((x_res1<(x1-r)) && (x_res1>(x0+r))){
				x[0] = x_res1;
				y[0] = ynew;
				found = true;
				goto stop;
			}else if((x_res2<(x1-r)) && (x_res2>(x0+r))){
				x[0] = x_res2;
				y[0] = ynew;
				found = true;
				goto stop;
			}
		}
	}
stop:
	return found;
}

//(x0, y0) as the container circle with radius r0
//to find the circle (x, y),r inscribed with (x0, y0) and circumscribed with (x1, y1), r1
bool tangentto2Circles_inscribe(const double x0, const double y0, const double r0, 
								const double x1, const double y1, const double r1, const double r, 
								double* x, double* y)
{
	double x02 = x0*x0;
	double x03 = pow(x0, 3);
	double x04 = pow(x0, 4);
	double x12 = x1*x1;
	double x13 = pow(x1, 3);
	double x14 = pow(x1, 4);
	double y02 = y0*y0;
	double y03 = pow(y0, 3);
	double y04 = pow(y0, 4);
	double y12 = y1*y1;
	double y13 = pow(y1, 3);
	double y14 = pow(y1, 4);
	double r02 = r0*r0;
	double r03 = pow(r0, 3);
	double r04 = pow(r0, 4);
	double r12 = r1*r1;
	double r13 = pow(r1, 3);
	double r14 = pow(r1, 4);
	double rr2 = r*r;

	double deter = -8 * x0 * x1 * r1 * r02 * r + 0.16e2 * x0 * x1 * r1 * r0 * rr2 + 8 * x0 * x1 * r12 * r0 * r - 0.15e2 * x04 * x12 + 2 * x04 * r12 + 6 * pow(x0, 0.5e1) * x1 - pow(x0, 6) - y14 * x02 - y14 * x12 - pow(x1, 6) - 8 * y1 * y03 * x0 * x1 + 0.12e2 * y12 * y02 * x0 * x1 - 4 * y02 * r02 * x0 * x1 + 8 * y0 * r02 * x0 * x1 * y1 + 2 * y02 * r12 * x02 + 2 * y02 * r02 * x12 - 8 * y1 * y0 * x02 * rr2 - 8 * y1 * y0 * x12 * rr2 + 4 * x02 * r1 * r02 * r - 8 * x02 * r1 * r0 * rr2 - 4 * x02 * r12 * r0 * r + 4 * x12 * r1 * r02 * r - 8 * x12 * r1 * r0 * rr2 - 4 * x12 * r12 * r0 * r - 4 * y0 * r02 * y1 * x02 - 4 * y0 * r02 * x12 * y1 - 4 * y02 * r12 * x0 * x1 - 4 * y0 * r12 * y1 * x02 - 4 * y0 * r12 * x12 * y1 - 8 * y13 * x0 * x1 * y0 + 4 * r1 * y12 * r * x02 + 4 * r1 * y12 * r * x12 + 0.24e2 * x02 * x12 * y1 * y0 - 0.16e2 * x03 * x1 * y1 * y0 - 4 * x0 * x1 * y12 * r02 - 4 * x0 * x1 * y12 * r12 - 0.16e2 * x0 * x13 * y1 * y0 + 4 * y02 * x02 * r1 * r - 4 * y02 * x02 * r0 * r + 4 * y1 * y03 * x02 - 6 * y12 * y02 * x02 + 4 * y1 * y03 * x12 - 6 * y12 * y02 * x12 + 2 * y02 * r02 * x02 + 6 * x0 * pow(x1, 0.5e1) + 8 * y0 * r12 * x0 * x1 * y1 + 2 * y02 * r12 * x12 + 2 * y14 * x0 * x1 + 4 * y13 * y0 * x02 + 4 * y13 * x12 * y0 + 2 * y04 * x0 * x1 - 0.12e2 * x02 * x12 * y12 + 8 * x03 * x1 * y12 + 8 * x0 * x13 * y12 - 0.12e2 * x02 * x12 * y02 + 8 * x03 * x1 * y02 + 8 * x0 * x13 * y02 + 4 * y0 * x04 * y1 + 2 * y12 * x02 * r02 + 2 * y12 * x02 * r12 + 2 * y12 * r02 * x12 + 2 * y12 * r12 * x12 - 8 * r1 * y12 * r * x0 * x1 + 0.16e2 * r1 * y1 * r * x0 * x1 * y0 - 8 * r1 * y1 * r * y0 * x02 - 8 * r1 * y1 * r * x12 * y0 - y04 * x02 - y04 * x12 - 2 * y02 * x04 - 2 * y12 * x04 - 2 * x14 * y02 - 2 * x14 * y12 - 0.15e2 * x02 * x14 + 0.20e2 * x03 * x13 + 2 * x04 * r02 + 4 * x04 * rr2 - x02 * r14 - x02 * r04 + 2 * x14 * r12 + 2 * x14 * r02 + 4 * x14 * rr2 - x12 * r14 - x12 * r04 + 8 * x0 * x1 * y12 * r0 * r - 0.16e2 * x0 * x1 * y1 * y0 * r0 * r - 8 * x0 * x1 * y02 * r1 * r + 8 * x0 * x1 * y02 * r0 * r + 8 * y0 * x02 * y1 * r0 * r + 4 * r1 * y02 * r * x12 - 4 * y12 * r0 * r * x02 - 4 * y12 * r0 * r * x12 - 4 * y02 * r0 * r * x12 - 0.24e2 * x02 * x12 * r0 * r - 4 * x0 * x1 * r12 * r02 + 0.24e2 * x02 * x12 * r1 * r + 0.16e2 * x03 * x1 * r0 * r - 0.16e2 * x0 * x13 * r1 * r - 0.16e2 * x03 * x1 * r1 * r + 0.16e2 * x0 * x13 * r0 * r + 8 * x0 * x1 * r12 * rr2 + 8 * x0 * x1 * r02 * rr2 - 8 * x0 * x1 * r03 * r + 8 * x0 * x1 * r13 * r + 8 * y1 * r0 * r * x12 * y0 + 4 * x14 * y0 * y1 - 8 * x03 * x1 * r12 - 8 * x03 * x1 * r02 - 8 * x0 * x13 * r12 - 8 * x0 * x13 * r02 + 0.12e2 * x02 * x12 * r02 + 0.12e2 * x02 * x12 * r12 + 0.24e2 * x02 * x12 * rr2 - 0.16e2 * x03 * x1 * rr2 - 0.16e2 * x0 * x13 * rr2 + 2 * x0 * x1 * r14 + 2 * x0 * x1 * r04 + 2 * x02 * r12 * r02 - 4 * x04 * r0 * r + 4 * x04 * r1 * r - 4 * x02 * r12 * rr2 - 4 * x02 * r02 * rr2 + 4 * x02 * r03 * r - 4 * x02 * r13 * r + 2 * x12 * r12 * r02 + 0.16e2 * y1 * y0 * x0 * x1 * rr2 + 4 * x14 * r1 * r - 4 * x14 * r0 * r - 4 * x12 * r12 * rr2 - 4 * x12 * r02 * rr2 + 4 * x12 * r03 * r - 4 * x12 * r13 * r + 4 * y12 * x02 * rr2 + 4 * y12 * x12 * rr2 + 4 * y02 * x02 * rr2 + 4 * y02 * x12 * rr2 - 8 * y12 * x0 * x1 * rr2 - 8 * y02 * x0 * x1 * rr2;
	double deter_rot = -1;
	if(deter < 0)
	{
		return false;
	}else{
		deter_rot = -8 * y0 * y1 * r1 * r02 * r + 0.16e2 * y0 * y1 * r1 * r0 * rr2 + 8 * y0 * y1 * r12 * r0 * r - 0.15e2 * y04 * y12 + 2 * y04 * r12 + 6 * pow(y0, 0.5e1) * y1 - pow(y0, 6) - x14 * y02 - x14 * y12 - pow(y1, 6) - 8 * x1 * x03 * y0 * y1 + 0.12e2 * x12 * x02 * y0 * y1 - 4 * x02 * r02 * y0 * y1 + 8 * x0 * r02 * y0 * y1 * x1 + 2 * x02 * r12 * y02 + 2 * x02 * r02 * y12 - 8 * x1 * x0 * y02 * rr2 - 8 * x1 * x0 * y12 * rr2 + 4 * y02 * r1 * r02 * r - 8 * y02 * r1 * r0 * rr2 - 4 * y02 * r12 * r0 * r + 4 * y12 * r1 * r02 * r - 8 * y12 * r1 * r0 * rr2 - 4 * y12 * r12 * r0 * r - 4 * x0 * r02 * x1 * y02 - 4 * x0 * r02 * y12 * x1 - 4 * x02 * r12 * y0 * y1 - 4 * x0 * r12 * x1 * y02 - 4 * x0 * r12 * y12 * x1 - 8 * x13 * y0 * y1 * x0 + 4 * r1 * x12 * r * y02 + 4 * r1 * x12 * r * y12 + 0.24e2 * y02 * y12 * x1 * x0 - 0.16e2 * y03 * y1 * x1 * x0 - 4 * y0 * y1 * x12 * r02 - 4 * y0 * y1 * x12 * r12 - 0.16e2 * y0 * y13 * x1 * x0 + 4 * x02 * y02 * r1 * r - 4 * x02 * y02 * r0 * r + 4 * x1 * x03 * y02 - 6 * x12 * x02 * y02 + 4 * x1 * x03 * y12 - 6 * x12 * x02 * y12 + 2 * x02 * r02 * y02 + 6 * y0 * pow(y1, 0.5e1) + 8 * x0 * r12 * y0 * y1 * x1 + 2 * x02 * r12 * y12 + 2 * x14 * y0 * y1 + 4 * x13 * x0 * y02 + 4 * x13 * y12 * x0 + 2 * x04 * y0 * y1 - 0.12e2 * y02 * y12 * x12 + 8 * y03 * y1 * x12 + 8 * y0 * y13 * x12 - 0.12e2 * y02 * y12 * x02 + 8 * y03 * y1 * x02 + 8 * y0 * y13 * x02 + 4 * x0 * y04 * x1 + 2 * x12 * y02 * r02 + 2 * x12 * y02 * r12 + 2 * x12 * r02 * y12 + 2 * x12 * r12 * y12 - 8 * r1 * x12 * r * y0 * y1 + 0.16e2 * r1 * x1 * r * y0 * y1 * x0 - 8 * r1 * x1 * r * x0 * y02 - 8 * r1 * x1 * r * y12 * x0 - x04 * y02 - x04 * y12 - 2 * x02 * y04 - 2 * x12 * y04 - 2 * y14 * x02 - 2 * y14 * x12 - 0.15e2 * y02 * y14 + 0.20e2 * y03 * y13 + 2 * y04 * r02 + 4 * y04 * rr2 - y02 * r14 - y02 * r04 + 2 * y14 * r12 + 2 * y14 * r02 + 4 * y14 * rr2 - y12 * r14 - y12 * r04 + 8 * y0 * y1 * x12 * r0 * r - 0.16e2 * y0 * y1 * x1 * x0 * r0 * r - 8 * y0 * y1 * x02 * r1 * r + 8 * y0 * y1 * x02 * r0 * r + 8 * x0 * y02 * x1 * r0 * r + 4 * r1 * x02 * r * y12 - 4 * x12 * r0 * r * y02 - 4 * x12 * r0 * r * y12 - 4 * x02 * r0 * r * y12 - 0.24e2 * y02 * y12 * r0 * r - 4 * y0 * y1 * r12 * r02 + 0.24e2 * y02 * y12 * r1 * r + 0.16e2 * y03 * y1 * r0 * r - 0.16e2 * y0 * y13 * r1 * r - 0.16e2 * y03 * y1 * r1 * r + 0.16e2 * y0 * y13 * r0 * r + 8 * y0 * y1 * r12 * rr2 + 8 * y0 * y1 * r02 * rr2 - 8 * y0 * y1 * r03 * r + 8 * y0 * y1 * r13 * r + 8 * x1 * r0 * r * y12 * x0 + 4 * y14 * x0 * x1 - 8 * y03 * y1 * r12 - 8 * y03 * y1 * r02 - 8 * y0 * y13 * r12 - 8 * y0 * y13 * r02 + 0.12e2 * y02 * y12 * r02 + 0.12e2 * y02 * y12 * r12 + 0.24e2 * y02 * y12 * rr2 - 0.16e2 * y03 * y1 * rr2 - 0.16e2 * y0 * y13 * rr2 + 2 * y0 * y1 * r14 + 2 * y0 * y1 * r04 + 2 * y02 * r12 * r02 - 4 * y04 * r0 * r + 4 * y04 * r1 * r - 4 * y02 * r12 * rr2 - 4 * y02 * r02 * rr2 + 4 * y02 * r03 * r - 4 * y02 * r13 * r + 2 * y12 * r12 * r02 + 0.16e2 * x1 * x0 * y0 * y1 * rr2 + 4 * y14 * r1 * r - 4 * y14 * r0 * r - 4 * y12 * r12 * rr2 - 4 * y12 * r02 * rr2 + 4 * y12 * r03 * r - 4 * y12 * r13 * r + 4 * x12 * y02 * rr2 + 4 * x12 * y12 * rr2 + 4 * x02 * y02 * rr2 + 4 * x02 * y12 * rr2 - 8 * x12 * y0 * y1 * rr2 - 8 * x02 * y0 * y1 * rr2;
		if(deter_rot < 0){
			return false;
		}
	}
	
	double C = -x12 - y12 + 2 * r * r1 + r12 + x02 + y02 - r02 + 2 * r0 * r;
	double denom = -2 * x0 * x1 - 2 * y1 * y0 + x02 + x12 + y12 + y02;
	double denom1 = 1.0/denom;

	double sqrtdeter, tmp, x_res1, y_res1, x_res2, y_res2;

	if (abs(x1-x0) < 1e-15)
	{
		if (deter_rot < 0)
		{
			deter_rot = -8 * y0 * y1 * r1 * r02 * r + 0.16e2 * y0 * y1 * r1 * r0 * rr2 + 8 * y0 * y1 * r12 * r0 * r - 0.15e2 * y04 * y12 + 2 * y04 * r12 + 6 * pow(y0, 0.5e1) * y1 - pow(y0, 6) - x14 * y02 - x14 * y12 - pow(y1, 6) - 8 * x1 * x03 * y0 * y1 + 0.12e2 * x12 * x02 * y0 * y1 - 4 * x02 * r02 * y0 * y1 + 8 * x0 * r02 * y0 * y1 * x1 + 2 * x02 * r12 * y02 + 2 * x02 * r02 * y12 - 8 * x1 * x0 * y02 * rr2 - 8 * x1 * x0 * y12 * rr2 + 4 * y02 * r1 * r02 * r - 8 * y02 * r1 * r0 * rr2 - 4 * y02 * r12 * r0 * r + 4 * y12 * r1 * r02 * r - 8 * y12 * r1 * r0 * rr2 - 4 * y12 * r12 * r0 * r - 4 * x0 * r02 * x1 * y02 - 4 * x0 * r02 * y12 * x1 - 4 * x02 * r12 * y0 * y1 - 4 * x0 * r12 * x1 * y02 - 4 * x0 * r12 * y12 * x1 - 8 * x13 * y0 * y1 * x0 + 4 * r1 * x12 * r * y02 + 4 * r1 * x12 * r * y12 + 0.24e2 * y02 * y12 * x1 * x0 - 0.16e2 * y03 * y1 * x1 * x0 - 4 * y0 * y1 * x12 * r02 - 4 * y0 * y1 * x12 * r12 - 0.16e2 * y0 * y13 * x1 * x0 + 4 * x02 * y02 * r1 * r - 4 * x02 * y02 * r0 * r + 4 * x1 * x03 * y02 - 6 * x12 * x02 * y02 + 4 * x1 * x03 * y12 - 6 * x12 * x02 * y12 + 2 * x02 * r02 * y02 + 6 * y0 * pow(y1, 0.5e1) + 8 * x0 * r12 * y0 * y1 * x1 + 2 * x02 * r12 * y12 + 2 * x14 * y0 * y1 + 4 * x13 * x0 * y02 + 4 * x13 * y12 * x0 + 2 * x04 * y0 * y1 - 0.12e2 * y02 * y12 * x12 + 8 * y03 * y1 * x12 + 8 * y0 * y13 * x12 - 0.12e2 * y02 * y12 * x02 + 8 * y03 * y1 * x02 + 8 * y0 * y13 * x02 + 4 * x0 * y04 * x1 + 2 * x12 * y02 * r02 + 2 * x12 * y02 * r12 + 2 * x12 * r02 * y12 + 2 * x12 * r12 * y12 - 8 * r1 * x12 * r * y0 * y1 + 0.16e2 * r1 * x1 * r * y0 * y1 * x0 - 8 * r1 * x1 * r * x0 * y02 - 8 * r1 * x1 * r * y12 * x0 - x04 * y02 - x04 * y12 - 2 * x02 * y04 - 2 * x12 * y04 - 2 * y14 * x02 - 2 * y14 * x12 - 0.15e2 * y02 * y14 + 0.20e2 * y03 * y13 + 2 * y04 * r02 + 4 * y04 * rr2 - y02 * r14 - y02 * r04 + 2 * y14 * r12 + 2 * y14 * r02 + 4 * y14 * rr2 - y12 * r14 - y12 * r04 + 8 * y0 * y1 * x12 * r0 * r - 0.16e2 * y0 * y1 * x1 * x0 * r0 * r - 8 * y0 * y1 * x02 * r1 * r + 8 * y0 * y1 * x02 * r0 * r + 8 * x0 * y02 * x1 * r0 * r + 4 * r1 * x02 * r * y12 - 4 * x12 * r0 * r * y02 - 4 * x12 * r0 * r * y12 - 4 * x02 * r0 * r * y12 - 0.24e2 * y02 * y12 * r0 * r - 4 * y0 * y1 * r12 * r02 + 0.24e2 * y02 * y12 * r1 * r + 0.16e2 * y03 * y1 * r0 * r - 0.16e2 * y0 * y13 * r1 * r - 0.16e2 * y03 * y1 * r1 * r + 0.16e2 * y0 * y13 * r0 * r + 8 * y0 * y1 * r12 * rr2 + 8 * y0 * y1 * r02 * rr2 - 8 * y0 * y1 * r03 * r + 8 * y0 * y1 * r13 * r + 8 * x1 * r0 * r * y12 * x0 + 4 * y14 * x0 * x1 - 8 * y03 * y1 * r12 - 8 * y03 * y1 * r02 - 8 * y0 * y13 * r12 - 8 * y0 * y13 * r02 + 0.12e2 * y02 * y12 * r02 + 0.12e2 * y02 * y12 * r12 + 0.24e2 * y02 * y12 * rr2 - 0.16e2 * y03 * y1 * rr2 - 0.16e2 * y0 * y13 * rr2 + 2 * y0 * y1 * r14 + 2 * y0 * y1 * r04 + 2 * y02 * r12 * r02 - 4 * y04 * r0 * r + 4 * y04 * r1 * r - 4 * y02 * r12 * rr2 - 4 * y02 * r02 * rr2 + 4 * y02 * r03 * r - 4 * y02 * r13 * r + 2 * y12 * r12 * r02 + 0.16e2 * x1 * x0 * y0 * y1 * rr2 + 4 * y14 * r1 * r - 4 * y14 * r0 * r - 4 * y12 * r12 * rr2 - 4 * y12 * r02 * rr2 + 4 * y12 * r03 * r - 4 * y12 * r13 * r + 4 * x12 * y02 * rr2 + 4 * x12 * y12 * rr2 + 4 * x02 * y02 * rr2 + 4 * x02 * y12 * rr2 - 8 * x12 * y0 * y1 * rr2 - 8 * x02 * y0 * y1 * rr2;
		}
		sqrtdeter = sqrt(deter_rot);
		double A1_rot = x1 * x02 + x0 * r02 - x0 * r12 - x13 + 2 * r1 * x1 * r - x03 + 2 * y0 * y1 * x1 + 2 * y0 * y1 * x0 - x0 * y02 - 2 * r1 * x0 * r + 2 * x1 * r0 * r - 2 * x0 * r0 * r + x12 * x0 - x1 * y02 - x1 * r02 + x1 * r12 - y12 * x0 - y12 * x1;
		tmp = 1.0/(y1-y0)*0.5; 

		y_res1 = -((A1_rot-sqrtdeter)*x0*denom1 - (A1_rot-sqrtdeter)*x1*denom1 + C)*tmp;
		x_res1 = -(A1_rot-sqrtdeter)*denom1*0.5;

		y_res2 = -((A1_rot+sqrtdeter)*x0*denom1 - (A1_rot+sqrtdeter)*x1*denom1 + C)*tmp;
		x_res2 = -(A1_rot+sqrtdeter)*denom1*0.5;
	}
	else{
		sqrtdeter = sqrt(deter);
		double A1 = y1 * y02 + y0 * r02 - y0 * r12 - y13 + 2 * r1 * y1 * r - y03 + 2 * x0 * x1 * y1 + 2 * x0 * x1 * y0 - y0 * x02 - 2 * r1 * y0 * r + 2 * y1 * r0 * r - 2 * y0 * r0 * r + y12 * y0 - y1 * x02 - y1 * r02 + y1 * r12 - x12 * y0 - x12 * y1;
		tmp = 1.0/(x1-x0)*0.5; 

		x_res1 = -((A1-sqrtdeter)*y0*denom1 - (A1-sqrtdeter)*y1*denom1 + C)*tmp;
		y_res1 = -(A1-sqrtdeter)*denom1*0.5;

		x_res2 = -((A1+sqrtdeter)*y0*denom1 - (A1+sqrtdeter)*y1*denom1 + C)*tmp;
		y_res2 = -(A1+sqrtdeter)*denom1*0.5;
	}

	x[0] = x_res1;  x[1] = x_res2;
	y[0] = y_res1;	y[1] = y_res2;

	return true;	
}

//to find the circle (x, y) circumscribed with (x1, y2) and (x2, y2) 
bool tangentto2Circles(double x1, double y1, double r1, double x2, double y2, double r2, double r, double* x, double* y)
{
	double x12 = x1*x1;
	double x13 = pow(x1, 3);
	double x14 = pow(x1, 4);
	double x22 = x2*x2;
	double x23 = pow(x2, 3);
	double x24 = pow(x2, 4);
	double y12 = y1*y1;
	double y13 = pow(y1, 3);
	double y14 = pow(y1, 4);
	double y22 = y2*y2;
	double y23 = pow(y2, 3);
	double y24 = pow(y2, 4);
	double r02 = r*r;
	double r12 = r1*r1;
	double r13 = pow(r1, 3);
	double r14 = pow(r1, 4);
	double r22 = r2*r2;
	double r23 = pow(r2, 3);
	double r24 = pow(r2, 4);

	double deter = -0.16e2 * y1 * x23 * x1 * y2 + 4 * y22 * r2 * r * x22 + 0.16e2 * x2 * x1 * y1 * r1 * y2 * r + 0.16e2 * x2 * x1 * y1 * y2 * r2 * r - 8 * x2 * x1 * y12 * r2 * r - 8 * x2 * x1 * y12 * r1 * r + 8 * x2 * x1 * y1 * y2 * r12 + 8 * x2 * x1 * y1 * y2 * r22 - 8 * x2 * x1 * y22 * r1 * r - 8 * x2 * x1 * y22 * r2 * r - 4 * x2 * x1 * y22 * r22 + 4 * r1 * y22 * r * x22 + 4 * y22 * r02 * x12 + 4 * y22 * r02 * x22 - 4 * x12 * r23 * r - 4 * x12 * r13 * r + 4 * x14 * r1 * r + 4 * x14 * r2 * r - 4 * x12 * r02 * r12 - 4 * x12 * r02 * r22 + 4 * y12 * r2 * r * x12 + 4 * r1 * y12 * r * x12 - pow(x2, 6) + 0.12e2 * x22 * x12 * r22 + 0.24e2 * x22 * x12 * r02 - 8 * x23 * x1 * r12 - 0.16e2 * x2 * x13 * r02 - 0.16e2 * x23 * x1 * r02 - 8 * x2 * x13 * r12 - 8 * x2 * x13 * r22 - 8 * x23 * x1 * r22 + 4 * y12 * r02 * x12 + 4 * y12 * r02 * x22 + 0.16e2 * y1 * y2 * x2 * x1 * r02 + 2 * x2 * x1 * r14 + 2 * x2 * x1 * r24 + 0.12e2 * x22 * x12 * r12 + 0.24e2 * y1 * x22 * x12 * y2 - 8 * y22 * x2 * x1 * r02 + 4 * x12 * r12 * r2 * r + 8 * x12 * r1 * r2 * r02 + 4 * x12 * r1 * r22 * r + 4 * x22 * r12 * r2 * r + 8 * x22 * r1 * r2 * r02 + 4 * x22 * r1 * r22 * r - 2 * y12 * x24 - 8 * y1 * x22 * r1 * y2 * r - 8 * y1 * x22 * y2 * r2 * r + 4 * y12 * x22 * r2 * r - 2 * x14 * y22 + 4 * y12 * x22 * r1 * r - 4 * y1 * x22 * y2 * r12 - 4 * y1 * x22 * y2 * r22 - 0.16e2 * x13 * y2 * x2 * y1 - y14 * x22 + 2 * y12 * x22 * r22 - 0.12e2 * y12 * x22 * x12 + 8 * x13 * y22 * x2 - 0.12e2 * x12 * y22 * x22 + 4 * x12 * y2 * y13 + 2 * x12 * y22 * r12 + 4 * x12 * y23 * y1 - 6 * x12 * y22 * y12 + 4 * x14 * y2 * y1 + 4 * y13 * x22 * y2 + 8 * y12 * x23 * x1 + 4 * y1 * x24 * y2 - 6 * y12 * x22 * y22 + 2 * y12 * x22 * r12 + 4 * y1 * x22 * y23 + 2 * y22 * x22 * r22 + 2 * y12 * r12 * x12 - y14 * x12 - 2 * x14 * y12 + 6 * x2 * pow(x1, 5) + 6 * pow(x2, 5) * x1 + 0.20e2 * x23 * x13 - 0.15e2 * x22 * x14 - 0.15e2 * x24 * x12 + 2 * y12 * r22 * x12 - x12 * y24 - y24 * x22 - 2 * y22 * x24 + 4 * x12 * y22 * r1 * r + 4 * x12 * y22 * r2 * r - 4 * x12 * y2 * y1 * r12 - 4 * x12 * y2 * y1 * r22 - 8 * x2 * x1 * y1 * y23 + 0.12e2 * x2 * x1 * y12 * y22 - 4 * x2 * x1 * y12 * r12 - 4 * x2 * x1 * y12 * r22 - 8 * x2 * x1 * y13 * y2 - 4 * x2 * x1 * y22 * r12 - 8 * x12 * y2 * y1 * r2 * r - 8 * y1 * y2 * r02 * x12 - 8 * y1 * y2 * r02 * x22 + 8 * x2 * x1 * r23 * r - 8 * x12 * y2 * r1 * y1 * r + 2 * x12 * y22 * r22 + 2 * x2 * x1 * y14 + 8 * x2 * x13 * y12 + 2 * x2 * x1 * y24 + 8 * x23 * x1 * y22 + 2 * y22 * x22 * r12 - 8 * x2 * x1 * r12 * r2 * r - 0.16e2 * x2 * x1 * r1 * r2 * r02 - 8 * x2 * x1 * r1 * r22 * r - x12 * r14 - x12 * r24 + 4 * x14 * r02 + 2 * x14 * r12 + 2 * x14 * r22 - x22 * r14 - x22 * r24 + 2 * x24 * r12 + 4 * x24 * r02 + 2 * x24 * r22 - pow(x1, 6) + 8 * x2 * x1 * r13 * r - 0.16e2 * x23 * x1 * r1 * r - 0.16e2 * x2 * x13 * r2 * r - 0.16e2 * x23 * x1 * r2 * r + 0.24e2 * x22 * x12 * r1 * r + 0.24e2 * x22 * x12 * r2 * r + 8 * x2 * x1 * r02 * r12 + 8 * x2 * x1 * r02 * r22 - 4 * x2 * x1 * r12 * r22 - 8 * y12 * x2 * x1 * r02 - 0.16e2 * x2 * x13 * r1 * r + 2 * x12 * r12 * r22 - 4 * x22 * r23 * r - 4 * x22 * r13 * r + 4 * x24 * r1 * r + 4 * x24 * r2 * r - 4 * x22 * r02 * r12 - 4 * x22 * r02 * r22 + 2 * x22 * r12 * r22;
	double deter_rot = -1;

	//if deter=0, that means the two input circles touch each other
	if(deter < 0)
	{
		return false;
	}else{
		deter_rot = -0.16e2 * x1 * y23 * y1 * x2 + 4 * x22 * r2 * r * y22 + 0.16e2 * y2 * y1 * x1 * r1 * x2 * r + 0.16e2 * y2 * y1 * x1 * x2 * r2 * r - 8 * y2 * y1 * x12 * r2 * r - 8 * y2 * y1 * x12 * r1 * r + 8 * y2 * y1 * x1 * x2 * r12 + 8 * y2 * y1 * x1 * x2 * r22 - 8 * y2 * y1 * x22 * r1 * r - 8 * y2 * y1 * x22 * r2 * r - 4 * y2 * y1 * x22 * r22 + 4 * r1 * x22 * r * y22 + 4 * x22 * r02 * y12 + 4 * x22 * r02 * y22 - 4 * y12 * r23 * r - 4 * y12 * r13 * r + 4 * y14 * r1 * r + 4 * y14 * r2 * r - 4 * y12 * r02 * r12 - 4 * y12 * r02 * r22 + 4 * x12 * r2 * r * y12 + 4 * r1 * x12 * r * y12 - pow(y2, 6) + 0.12e2 * y22 * y12 * r22 + 0.24e2 * y22 * y12 * r02 - 8 * y23 * y1 * r12 - 0.16e2 * y2 * y13 * r02 - 0.16e2 * y23 * y1 * r02 - 8 * y2 * y13 * r12 - 8 * y2 * y13 * r22 - 8 * y23 * y1 * r22 + 4 * x12 * r02 * y12 + 4 * x12 * r02 * y22 + 0.16e2 * x1 * x2 * y2 * y1 * r02 + 2 * y2 * y1 * r14 + 2 * y2 * y1 * r24 + 0.12e2 * y22 * y12 * r12 + 0.24e2 * x1 * y22 * y12 * x2 - 8 * x22 * y2 * y1 * r02 + 4 * y12 * r12 * r2 * r + 8 * y12 * r1 * r2 * r02 + 4 * y12 * r1 * r22 * r + 4 * y22 * r12 * r2 * r + 8 * y22 * r1 * r2 * r02 + 4 * y22 * r1 * r22 * r - 2 * x12 * y24 - 8 * x1 * y22 * r1 * x2 * r - 8 * x1 * y22 * x2 * r2 * r + 4 * x12 * y22 * r2 * r - 2 * y14 * x22 + 4 * x12 * y22 * r1 * r - 4 * x1 * y22 * x2 * r12 - 4 * x1 * y22 * x2 * r22 - 0.16e2 * y13 * x2 * y2 * x1 - x14 * y22 + 2 * x12 * y22 * r22 - 0.12e2 * x12 * y22 * y12 + 8 * y13 * x22 * y2 - 0.12e2 * y12 * x22 * y22 + 4 * y12 * x2 * x13 + 2 * y12 * x22 * r12 + 4 * y12 * x23 * x1 - 6 * y12 * x22 * x12 + 4 * y14 * x2 * x1 + 4 * x13 * y22 * x2 + 8 * x12 * y23 * y1 + 4 * x1 * y24 * x2 - 6 * x12 * y22 * x22 + 2 * x12 * y22 * r12 + 4 * x1 * y22 * x23 + 2 * x22 * y22 * r22 + 2 * x12 * r12 * y12 - x14 * y12 - 2 * y14 * x12 + 6 * y2 * pow(y1, 5) + 6 * pow(y2, 5) * y1 + 0.20e2 * y23 * y13 - 0.15e2 * y22 * y14 - 0.15e2 * y24 * y12 + 2 * x12 * r22 * y12 - y12 * x24 - x24 * y22 - 2 * x22 * y24 + 4 * y12 * x22 * r1 * r + 4 * y12 * x22 * r2 * r - 4 * y12 * x2 * x1 * r12 - 4 * y12 * x2 * x1 * r22 - 8 * y2 * y1 * x1 * x23 + 0.12e2 * y2 * y1 * x12 * x22 - 4 * y2 * y1 * x12 * r12 - 4 * y2 * y1 * x12 * r22 - 8 * y2 * y1 * x13 * x2 - 4 * y2 * y1 * x22 * r12 - 8 * y12 * x2 * x1 * r2 * r - 8 * x1 * x2 * r02 * y12 - 8 * x1 * x2 * r02 * y22 + 8 * y2 * y1 * r23 * r - 8 * y12 * x2 * r1 * x1 * r + 2 * y12 * x22 * r22 + 2 * y2 * y1 * x14 + 8 * y2 * y13 * x12 + 2 * y2 * y1 * x24 + 8 * y23 * y1 * x22 + 2 * x22 * y22 * r12 - 8 * y2 * y1 * r12 * r2 * r - 0.16e2 * y2 * y1 * r1 * r2 * r02 - 8 * y2 * y1 * r1 * r22 * r - y12 * r14 - y12 * r24 + 4 * y14 * r02 + 2 * y14 * r12 + 2 * y14 * r22 - y22 * r14 - y22 * r24 + 2 * y24 * r12 + 4 * y24 * r02 + 2 * y24 * r22 - pow(y1, 6) + 8 * y2 * y1 * r13 * r - 0.16e2 * y23 * y1 * r1 * r - 0.16e2 * y2 * y13 * r2 * r - 0.16e2 * y23 * y1 * r2 * r + 0.24e2 * y22 * y12 * r1 * r + 0.24e2 * y22 * y12 * r2 * r + 8 * y2 * y1 * r02 * r12 + 8 * y2 * y1 * r02 * r22 - 4 * y2 * y1 * r12 * r22 - 8 * x12 * y2 * y1 * r02 - 0.16e2 * y2 * y13 * r1 * r + 2 * y12 * r12 * r22 - 4 * y22 * r23 * r - 4 * y22 * r13 * r + 4 * y24 * r1 * r + 4 * y24 * r2 * r - 4 * y22 * r02 * r12 - 4 * y22 * r02 * r22 + 2 * y22 * r12 * r22;
		if(deter_rot < 0){
			return false;
		}
	}
	if (abs(x1-x2) < 1e-15)
	{
		if ((abs(y1-y2) - r1-r2)>2*r)
		{
			return false;
		}
	}
	double denom = x12 + y12 + x22 + y22 - 2 * y1 * y2 -2* x1 * x2;
	double C = -x12 - y12 + 2 * r * r1 + r12 + x22 + y22 - 2 * r * r2 - r22;	
	double denom1 = 1.0/denom;

	double sqrtdeter, tmp, x_res1, y_res1, x_res2, y_res2;

	if (abs(x1-x2) < 1e-15)
	{
		if (deter_rot < 0)
		{
			deter_rot = -0.16e2 * x1 * y23 * y1 * x2 + 4 * x22 * r2 * r * y22 + 0.16e2 * y2 * y1 * x1 * r1 * x2 * r + 0.16e2 * y2 * y1 * x1 * x2 * r2 * r - 8 * y2 * y1 * x12 * r2 * r - 8 * y2 * y1 * x12 * r1 * r + 8 * y2 * y1 * x1 * x2 * r12 + 8 * y2 * y1 * x1 * x2 * r22 - 8 * y2 * y1 * x22 * r1 * r - 8 * y2 * y1 * x22 * r2 * r - 4 * y2 * y1 * x22 * r22 + 4 * r1 * x22 * r * y22 + 4 * x22 * r02 * y12 + 4 * x22 * r02 * y22 - 4 * y12 * r23 * r - 4 * y12 * r13 * r + 4 * y14 * r1 * r + 4 * y14 * r2 * r - 4 * y12 * r02 * r12 - 4 * y12 * r02 * r22 + 4 * x12 * r2 * r * y12 + 4 * r1 * x12 * r * y12 - pow(y2, 6) + 0.12e2 * y22 * y12 * r22 + 0.24e2 * y22 * y12 * r02 - 8 * y23 * y1 * r12 - 0.16e2 * y2 * y13 * r02 - 0.16e2 * y23 * y1 * r02 - 8 * y2 * y13 * r12 - 8 * y2 * y13 * r22 - 8 * y23 * y1 * r22 + 4 * x12 * r02 * y12 + 4 * x12 * r02 * y22 + 0.16e2 * x1 * x2 * y2 * y1 * r02 + 2 * y2 * y1 * r14 + 2 * y2 * y1 * r24 + 0.12e2 * y22 * y12 * r12 + 0.24e2 * x1 * y22 * y12 * x2 - 8 * x22 * y2 * y1 * r02 + 4 * y12 * r12 * r2 * r + 8 * y12 * r1 * r2 * r02 + 4 * y12 * r1 * r22 * r + 4 * y22 * r12 * r2 * r + 8 * y22 * r1 * r2 * r02 + 4 * y22 * r1 * r22 * r - 2 * x12 * y24 - 8 * x1 * y22 * r1 * x2 * r - 8 * x1 * y22 * x2 * r2 * r + 4 * x12 * y22 * r2 * r - 2 * y14 * x22 + 4 * x12 * y22 * r1 * r - 4 * x1 * y22 * x2 * r12 - 4 * x1 * y22 * x2 * r22 - 0.16e2 * y13 * x2 * y2 * x1 - x14 * y22 + 2 * x12 * y22 * r22 - 0.12e2 * x12 * y22 * y12 + 8 * y13 * x22 * y2 - 0.12e2 * y12 * x22 * y22 + 4 * y12 * x2 * x13 + 2 * y12 * x22 * r12 + 4 * y12 * x23 * x1 - 6 * y12 * x22 * x12 + 4 * y14 * x2 * x1 + 4 * x13 * y22 * x2 + 8 * x12 * y23 * y1 + 4 * x1 * y24 * x2 - 6 * x12 * y22 * x22 + 2 * x12 * y22 * r12 + 4 * x1 * y22 * x23 + 2 * x22 * y22 * r22 + 2 * x12 * r12 * y12 - x14 * y12 - 2 * y14 * x12 + 6 * y2 * pow(y1, 5) + 6 * pow(y2, 5) * y1 + 0.20e2 * y23 * y13 - 0.15e2 * y22 * y14 - 0.15e2 * y24 * y12 + 2 * x12 * r22 * y12 - y12 * x24 - x24 * y22 - 2 * x22 * y24 + 4 * y12 * x22 * r1 * r + 4 * y12 * x22 * r2 * r - 4 * y12 * x2 * x1 * r12 - 4 * y12 * x2 * x1 * r22 - 8 * y2 * y1 * x1 * x23 + 0.12e2 * y2 * y1 * x12 * x22 - 4 * y2 * y1 * x12 * r12 - 4 * y2 * y1 * x12 * r22 - 8 * y2 * y1 * x13 * x2 - 4 * y2 * y1 * x22 * r12 - 8 * y12 * x2 * x1 * r2 * r - 8 * x1 * x2 * r02 * y12 - 8 * x1 * x2 * r02 * y22 + 8 * y2 * y1 * r23 * r - 8 * y12 * x2 * r1 * x1 * r + 2 * y12 * x22 * r22 + 2 * y2 * y1 * x14 + 8 * y2 * y13 * x12 + 2 * y2 * y1 * x24 + 8 * y23 * y1 * x22 + 2 * x22 * y22 * r12 - 8 * y2 * y1 * r12 * r2 * r - 0.16e2 * y2 * y1 * r1 * r2 * r02 - 8 * y2 * y1 * r1 * r22 * r - y12 * r14 - y12 * r24 + 4 * y14 * r02 + 2 * y14 * r12 + 2 * y14 * r22 - y22 * r14 - y22 * r24 + 2 * y24 * r12 + 4 * y24 * r02 + 2 * y24 * r22 - pow(y1, 6) + 8 * y2 * y1 * r13 * r - 0.16e2 * y23 * y1 * r1 * r - 0.16e2 * y2 * y13 * r2 * r - 0.16e2 * y23 * y1 * r2 * r + 0.24e2 * y22 * y12 * r1 * r + 0.24e2 * y22 * y12 * r2 * r + 8 * y2 * y1 * r02 * r12 + 8 * y2 * y1 * r02 * r22 - 4 * y2 * y1 * r12 * r22 - 8 * x12 * y2 * y1 * r02 - 0.16e2 * y2 * y13 * r1 * r + 2 * y12 * r12 * r22 - 4 * y22 * r23 * r - 4 * y22 * r13 * r + 4 * y24 * r1 * r + 4 * y24 * r2 * r - 4 * y22 * r02 * r12 - 4 * y22 * r02 * r22 + 2 * y22 * r12 * r22;
		}
		sqrtdeter = sqrt(deter_rot);
		double A1_rot = -x1 * y22 - y12 * x2 + 2 * y2 * y1 * x1 + 2 * y2 * y1 * x2 - 2 * r1 * x2 * r + 2 * x2 * r2 * r - 2 * x1 * r2 * r + 2 * r1 * x1 * r - x23 - x2 * y22 - x13 - x2 * r12 + x1 * x22 + x1 * r12 - x1 * r22 + x12 * x2 - y12 * x1 + x2 * r22;
		tmp = 1.0/(y1-y2)*0.5; 

		y_res1 = -((A1_rot-sqrtdeter)*x2*denom1 - (A1_rot-sqrtdeter)*x1*denom1 + C)*tmp;
		x_res1 = -(A1_rot-sqrtdeter)*denom1*0.5;

		y_res2 = -((A1_rot+sqrtdeter)*x2*denom1 - (A1_rot+sqrtdeter)*x1*denom1 + C)*tmp;
		x_res2 = -(A1_rot+sqrtdeter)*denom1*0.5;
	}
	else{
		sqrtdeter = sqrt(deter);
		double A1 = -y1 * x22 - x12 * y2 + 2 * x2 * x1 * y1 + 2 * x2 * x1 * y2 - 2 * r1 * y2 * r + 2 * y2 * r2 * r - 2 * y1 * r2 * r + 2 * r1 * y1 * r - y23 - y2 * x22 - y13 - y2 * r12 + y1 * y22 + y1 * r12 - y1 * r22 + y12 * y2 - x12 * y1 + y2 * r22;
		tmp = 1.0/(x1-x2)*0.5; 

		x_res1 = -((A1-sqrtdeter)*y2*denom1 - (A1-sqrtdeter)*y1*denom1 + C)*tmp;
		y_res1 = -(A1-sqrtdeter)*denom1*0.5;

		x_res2 = -((A1+sqrtdeter)*y2*denom1 - (A1+sqrtdeter)*y1*denom1 + C)*tmp;
		y_res2 = -(A1+sqrtdeter)*denom1*0.5;
	}
	
	x[0] = x_res1;  x[1] = x_res2;
	y[0] = y_res1;	y[1] = y_res2;

	return true;	
}

bool isInsidePolygon(const Point* polygon, int N, Point pt)
{
	int counter = 0;
	int i;
	double xinters;
	Point pt1,pt2;

	pt1 = polygon[0];
	for (i=1;i<=N;i++) {
		pt2 = polygon[i % N];
		if (pt.x() > std::min(pt1.y(),pt2.y())) {
			if (pt.y() <= std::max(pt1.y(),pt2.y())) {
				if (pt.x() <= std::max(pt1.x(),pt2.x())) {
					if (pt1.y() != pt2.y()) {
						xinters = (pt.y()-pt1.y())*(pt2.x()-pt1.x())/(pt2.y()-pt1.y())+pt1.x();
						if (pt1.x() == pt2.x() || pt.x() <= xinters)
							counter++;
					}
				}
			}
		}
		pt1 = pt2;
	}

	if (counter % 2 == 0)
		return false;
	else
		return true;
}

bool intersec2Circles(double x1, double y1, double r1, double x2, double y2, double r2,  Point2D* pts, int& num)
{
	if ((abs(x1-x2) < 1e-15) && (abs(y1-y2) < 1e-15) && (abs(r1-r2) < 1e-15))
	{// circles are essentially the same
		num = 0;
		return false;
	}
	
	double x12 = x1 * x1;
	double x13 =  pow(x1, 3);
	double x14 =  pow(x1, 4);
	double x22 = x2 * x2;
	double x23 =  pow(x2, 3);
	double x24 =  pow(x2, 4);
	double y12 = y1 * y1;
	double y13 =  pow(y1, 3);
	double y14 =  pow(y1, 4);
	double y22 = y2 * y2;
	double y23 =  pow(y2, 3);
	double y24 =  pow(y2, 4);
	double r12 = r1 * r1;
	double r14 =  pow(r1, 4);
	double r22 = r2 * r2;
	double r24 =  pow(r2, 4);

	double deter = -0.15e2 * x24 * x12 + 6 * x2 * pow(x1, 5) + 0.20e2 * x23 * x13 + 6 * pow(x2, 5) * x1 - y24 * x12 + 4 * y23 * x12 * y1 + 8 * r12 * y2 * x2 * x1 * y1 - y24 * x22 - x12 * r14 - x12 * r24 + 2 * x24 * r12 + 2 * x24 * r22 - x22 * r14 - x22 * r24 - 2 * y12 * x24 - y14 * x22 - 2 * x24 * y22 - x12 * y14 + 8 * y12 * x23 * x1 + 2 * y12 * x22 * r22 + 2 * y12 * x22 * r12 + 4 * y1 * x24 * y2 + 2 * x2 * x1 * y14 + 2 * r22 * y22 * x22 + 8 * x23 * x1 * y22 - 8 * x23 * x1 * r12 + 0.12e2 * x22 * x12 * r12 - 8 * x2 * x13 * r22 - 8 * x2 * x13 * r12 - 8 * x23 * x1 * r22 + 2 * x2 * x1 * r14 + 2 * x2 * x1 * r24 + 0.12e2 * x22 * x12 * r22 + 8 * x2 * x1 * y1 * r22 * y2 + 2 * x12 * r12 * r22 + 2 * x22 * r12 * r22 + 2 * x14 * r22 + 2 * r12 * x14 - 2 * x14 * y12 + 2 * r12 * y22 * x12 + 2 * r12 * y22 * x22 - 6 * y12 * y22 * x12 + 4 * y13 * y2 * x12 + 4 * y13 * y2 * x22 - 6 * y12 * y22 * x22 + 4 * y23 * y1 * x22 + 2 * y24 * x2 * x1 + 4 * x14 * y2 * y1 + 2 * x12 * y22 * r22 + 8 * x13 * y22 * x2 - 0.12e2 * x12 * y22 * x22 - 0.12e2 * x12 * y12 * x22 + 8 * x13 * y12 * x2 + 2 * x12 * y12 * r22 + 2 * x12 * y12 * r12 - 2 * x14 * y22 - 0.15e2 * x14 * x22 - pow(x1, 6) - pow(x2, 6) - 4 * r12 * y2 * x12 * y1 - 4 * r12 * y2 * y1 * x22 - 4 * r12 * y22 * x2 * x1 - 8 * y13 * y2 * x2 * x1 + 0.12e2 * y12 * y22 * x2 * x1 - 8 * y23 * x2 * x1 * y1 + 0.24e2 * x12 * y2 * y1 * x22 - 0.16e2 * x13 * y2 * x2 * y1 - 4 * x12 * y2 * y1 * r22 - 4 * y1 * x22 * r22 * y2 - 0.16e2 * y1 * x23 * x1 * y2 - 4 * x2 * x1 * y12 * r22 - 4 * x2 * x1 * y12 * r12 - 4 * r22 * y22 * x2 * x1 - 4 * x2 * x1 * r12 * r22;
	double A = r12 * y2 - y12 * y2 + y23 + x12 * y2 - y1 * y22 + x12 * y1 + y1 * x22 - 2 * x2 * x1 * y1 - r22 * y2 + y1 * r22 - r12 * y1 + y13 - 2 * x2 * x1 * y2 + x22 * y2;
	
	double denom = -2*x2*x1 - 2*y1*y2 + y12 + y22 + x12 + x22;
	double denom1 = 1.0/denom;
	double C = x12 + y12 - y22 - x22 -r12 + r22;

	double deter_rot, A_rot, tmp;
	if (abs(x1-x2) < 1e-15)
	{
		deter_rot = -0.15e2 * y24 * y12 + 6 * y2 * pow(y1, 5) + 0.20e2 * y23 * y13 + 6 * pow(y2, 5) * y1 - x24 * y12 + 4 * x23 * y12 * x1 + 8 * r12 * x2 * y2 * y1 * x1 - x24 * y22 - y12 * r14 - y12 * r24 + 2 * y24 * r12 + 2 * y24 * r22 - y22 * r14 - y22 * r24 - 2 * x12 * y24 - x14 * y22 - 2 * y24 * x22 - y12 * x14 + 8 * x12 * y23 * y1 + 2 * x12 * y22 * r22 + 2 * x12 * y22 * r12 + 4 * x1 * y24 * x2 + 2 * y2 * y1 * x14 + 2 * r22 * x22 * y22 + 8 * y23 * y1 * x22 - 8 * y23 * y1 * r12 + 0.12e2 * y22 * y12 * r12 - 8 * y2 * y13 * r22 - 8 * y2 * y13 * r12 - 8 * y23 * y1 * r22 + 2 * y2 * y1 * r14 + 2 * y2 * y1 * r24 + 0.12e2 * y22 * y12 * r22 + 8 * y2 * y1 * x1 * r22 * x2 + 2 * y12 * r12 * r22 + 2 * y22 * r12 * r22 + 2 * y14 * r22 + 2 * r12 * y14 - 2 * y14 * x12 + 2 * r12 * x22 * y12 + 2 * r12 * x22 * y22 - 6 * x12 * x22 * y12 + 4 * x13 * x2 * y12 + 4 * x13 * x2 * y22 - 6 * x12 * x22 * y22 + 4 * x23 * x1 * y22 + 2 * x24 * y2 * y1 + 4 * y14 * x2 * x1 + 2 * y12 * x22 * r22 + 8 * y13 * x22 * y2 - 0.12e2 * y12 * x22 * y22 - 0.12e2 * y12 * x12 * y22 + 8 * y13 * x12 * y2 + 2 * y12 * x12 * r22 + 2 * y12 * x12 * r12 - 2 * y14 * x22 - 0.15e2 * y14 * y22 - pow(y1, 6) - pow(y2, 6) - 4 * r12 * x2 * y12 * x1 - 4 * r12 * x2 * x1 * y22 - 4 * r12 * x22 * y2 * y1 - 8 * x13 * x2 * y2 * y1 + 0.12e2 * x12 * x22 * y2 * y1 - 8 * x23 * y2 * y1 * x1 + 0.24e2 * y12 * x2 * x1 * y22 - 0.16e2 * y13 * x2 * y2 * x1 - 4 * y12 * x2 * x1 * r22 - 4 * x1 * y22 * r22 * x2 - 0.16e2 * x1 * y23 * y1 * x2 - 4 * y2 * y1 * x12 * r22 - 4 * y2 * y1 * x12 * r12 - 4 * r22 * x22 * y2 * y1 - 4 * y2 * y1 * r12 * r22;
		A_rot = r12 * x2 - x12 * x2 + x23 + y12 * x2 - x1 * x22 + y12 * x1 + x1 * y22 - 2 * y2 * y1 * x1 - r22 * x2 + x1 * r22 - r12 * x1 + x13 - 2 * y2 * y1 * x2 + y22 * x2;
		tmp = 1.0/(y1-y2)*0.5;

		if (deter_rot > /*TOLERLOW*/TOLER)
		{//two intersections
			double sqrtdeter = sqrt(deter_rot);
			double y = (-(A_rot+sqrtdeter)*x1*denom1 + (A_rot+sqrtdeter)*x2*denom1 + C)*tmp;
			double x = (A_rot+sqrtdeter)*denom1*0.5;
			pts[0] = Point2D(x,y);

			y = (-(A_rot-sqrtdeter)*x1*denom1 + (A_rot-sqrtdeter)*x2*denom1 + C)*tmp;
			x = (A_rot-sqrtdeter)*denom1*0.5;
			pts[1] = Point2D(x,y);
			num = 2;
		}else if (abs(deter) < /*TOLERLOW*/TOLER)
		{
			double y = (-(A_rot)*x1*denom1 + (A_rot)*x2*denom1 + C)*tmp;
			double x = A_rot*denom1*0.5;
			pts[0] = Point2D(x,y);
			num = 1;
		}else{
			num = 0;
			return false;
		}

	}else{
		tmp = 1.0/(x1-x2)*0.5;

		if (deter > /*TOLERLOW*/TOLER)
		{//two intersections
			double sqrtdeter = sqrt(deter);
			double x = (-(A+sqrtdeter)*y1*denom1 + (A+sqrtdeter)*y2*denom1 + C)*tmp;
			double y = (A+sqrtdeter)*denom1*0.5;
			pts[0] = Point2D(x,y);

			x = (-(A-sqrtdeter)*y1*denom1 + (A-sqrtdeter)*y2*denom1 + C)*tmp;
			y = (A-sqrtdeter)*denom1*0.5;
			pts[1] = Point2D(x,y);
			num = 2;
		}else if (abs(deter) < /*TOLERLOW*/TOLER)
		{		
			double x = (-(A)*y1*denom1 + (A)*y2*denom1 + C)*tmp;
			double y = A*denom1*0.5;
			pts[0] = Point2D(x,y);
			num = 1;
		}else{
			num = 0;
			return false;
		}
	}
	return true;
}

//to determine whether (xc, yc) lies on the line (xa,ya) and (xb, yb) 
bool IsPtOnLine(double XA,double YA,double XB,double YB,double XC,double YC)
{
	double l = sqrt((XB-XA)*(XB-XA)+(YB-YA)*(YB-YA));
	double ll = (XB-XA)*(XB-XA)+(YB-YA)*(YB-YA);
	bool res = false;
	if (l!=0){
		double rr = ((YA-YC)*(YA-YB)-(XA-XC)*(XB-XA))/ll;
		double ss = ((YA-YC)*(XB-XA)-(XA-XC)*(YB-YA))/ll;
		if ((rr>0) && (rr<1)){
			double dist = ss*l;
			if(abs(dist) <= TOLER){
				res = true;
			}
		}
	}
	return res;
}

bool IsPtOnSeg(Point P, Point P1, Point P2)
{
	if (IsSamePoint(P1,P2) && IsSamePoint(P,P1)) return true;
	else if (IsSamePoint(P1,P2)) return false;

	if (fabs((P.x()-P1.x())*(P2.y()-P1.y())-(P2.x()-P1.x())*(P.y()-P1.y()))<TOLER)///if (fabs((P.x()-P1.x())/(P2.x()-P1.x())-(P.y()-P1.y())/(P2.y()-P1.y()))<TOLER)
		return true;
	else
		return false;
}


double disPt2Line(double x1,double y1,double x2,double y2,double x0,double y0, double& x, double& y)
{	
	double denom = (x2-x1)*(x2-x1) +(y2-y1)*(y2-y1);
	
	double res = 0;
	double u = 0;
	if (denom != 0){
		u = ((x0-x1)*(x2-x1) + (y0-y1)*(y2-y1))/denom;
	}
	x = x1 + u * (x2 - x1);
	y = y1 + u * (y2 - y1);

	res = sqrt((x0-x)*(x0-x) + (y0-y)*(y0-y));
	return res;
}

//(x1,y1)->(x2,y2)
//return positive value if (x0,y0) lies on the normal direction of the line
double signedDisPt2Line(double x1,double y1,double x2,double y2,double x0,double y0, double& x, double& y)
{	
	double denom = (x2-x1)*(x2-x1) +(y2-y1)*(y2-y1);

	double res = 0;
	double u = 0;
	if (denom != 0){
		u = ((x0-x1)*(x2-x1) + (y0-y1)*(y2-y1))/denom;
	}
	x = x1 + u * (x2 - x1);
	y = y1 + u * (y2 - y1);

	res = sqrt((x0-x)*(x0-x) + (y0-y)*(y0-y));

	double kx = -(y2 - y1);
	double ky = (x2 - x1);
	Vector2D vecK(kx, ky);
	vecK.SetUnit();

	Vector2D vecD(x0-x, y0-y);
	vecD.SetUnit();
	double tmp = vecD*vecK;
	int sign;
	if (tmp > 0)
	{		
		sign = 1;
	}else{
		sign = -1;
	}

	return res*sign;
}

//arc: center(cx, cy)  radius r  start point (x1, y1)   end point(x2, y2)
bool Arc2Contains (double cx, double cy, double r, double x1, double y1, double x2, double y2, Point2D pt)
{
	// Assert: |P-C| = R where P is the input point, C is the circle center,
	// and R is the circle radius.  For P to be on the arc from A to B, it
	// must be on the side of the plane containing A with normal N = Perp(B-A)
	// where Perp(u,v) = (v,-u).
	 
	Vector2D kPmE0(pt.x-x1, pt.y-y1);
	Vector2D kE1mE0(x2-x1, y2-y1);
	double fDotPerp = kPmE0.x*kE1mE0.y - kPmE0.y*kE1mE0.x;
	return fDotPerp >= 0.0;
}

//circle center(cx0, cy0) radius r0
//arc: center(cx1, cy1)  radius r1  start point (x1, y1)   end point(x2, y2)
bool IntrArc2Circle2(double cx0, double cy0, double r0, double cx1, double cy1, double r1, double x1, double y1, double x2, double y2, Point2D* pts, int& num)
{

	Point2D _pts[2];
	int ires;
	bool inter = intersec2Circles(cx0, cy0, r0, cx1, cy1, r1, _pts, ires);

	if (!inter)
	{
		return false;
	}	

	// test if circle-circle intersection points are on the arc
	int id = 0;
	for (int i = 0; i < ires; i++)
	{
		if (Arc2Contains(cx1, cy1, r1, x1, y1, x2, y2, _pts[i]))
		{
			pts[id] = _pts[i];
			id++;
		}
	}
	num = id;
	if (num >0)
	{
		return true;
	}
	return false;	
}

bool LinesCross(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double& x, double& y)
{
	double k1,k2;
	if(x0 == x1)
	{
		if(x2 == x3)
		{
			return false;
		}
		else
		{
			x = x0;
			y = (y3*(x0-x2) + y2*(x3-x0))/(x3-x2);
			return true;
		}
	}
	else if(x2 == x3)
	{
		x = x2;
		y = (y1*(x2-x0) + y0*(x1-x2))/(x1-x0);
		return true;
	}
	else
	{
		k1=(y0-y1)/(x0-x1);
		k2=(y2-y3)/(x2-x3);
		if(k1 == k2)
			return false;
		else
		{
			x=(k1*x0-k2*x2+y2-y0)/(k1-k2);
			y=y0+(x-x0)*k1;
			return true;
		}
	}
	
}

void AngleSplit(double x0, double y0, double x1, double y1, double x2, double y2, double& x, double& y)//angle 012
{//
	double dTmp = (y0-y1)*std::sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)) 
		-(y2-y1)*std::sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
	if(0 == dTmp)
	{
		x = x1;
		y = y1 + 5;
	}
	else
	{
		dTmp = ( (x2-x1)*std::sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1))
			-(x0-x1)*std::sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)) )/ dTmp;
		x = x1 + 5;
		y = y1 + 5*dTmp;
	}
}

double TriEdgeCircleR(const Point point1, const Point point2, const Point point3, const Point point4, Point& pointCenter)
{
	double dX1, dY1, dX2, dY2;
	
	//point3 - point1 - point2 - point4
	double dXTmp1, dYTmp1;
	double dXTmp2, dYTmp2;
	AngleSplit( point3.x(), point3.y(), point1.x(), point1.y(), 
		point2.x(), point2.y(), dXTmp1, dYTmp1);
	AngleSplit(point1.x(), point1.y(), point2.x(), point2.y(),
		point4.x(), point4.y(), dXTmp2, dYTmp2);

	double dX, dY;
	if(LinesCross(dX1, dY1, dX2, dY2, dXTmp1, dYTmp1, dXTmp2, dYTmp2, dX, dY))
	{
		pointCenter = Point(dX, dY);
		return disPt2Line(dX1, dY1, dX2, dY2, dX, dY, dXTmp1, dYTmp1);
	}
	else
	{
		return -1;
	}
	
}