//#include "../stdafx.h"
#include "Polygon.h"
#include "UtilFunc.h"

CPolygon::CPolygon()
{
    count=0;
    acmes=NULL;
}

CPolygon::~CPolygon()
{
    if (acmes!=NULL)
    {
        delete[] acmes;
        acmes=NULL;
    }

}

double CPolygon::TriangleArea(Point p1, Point p2, Point p3)
{
    double area = (p1.x() * p2.y() + p2.x() * p3.y() + p3.x() * p1.y())
        - (p1.x() * p3.y() + p2.x() * p1.y() + p3.x() * p2.y());

    return area / 2;
}

// 计算多边形面积(给定的顶点序列必须能组成多边形)
double CPolygon::PolyArea()
{
    Point p0(0.0, 0.0);

    double area = 0.0;
    for(int i=0; i<this->count; ++i)
        area += TriangleArea(p0, this->acmes[i], this->acmes[i+1]);

    return area > 0.0 ? area : -area;
}

// 判断由（x,y）构成的点 pt 在有向线段 p1p2 的哪一侧)
// 返回值  1 左侧、0 在直线上、-1 右侧
int CPolygon::SideOfLine(Point p1, Point p2, double x, double y)
{
    enum POS {RR =-1, ON = 0, LL = 1}; // 右侧、线上、左侧

    // 叉积的两个子项
    double c1 = (p2.x() - x) * (y - p1.y());
    double c2 = (p2.y() - y) * (x - p1.x());

    return c1 > c2 ? RR : (c1 < c2 ? LL : ON);
}

/*
 * 判断给定的顶点是否能组成多边形(凹多边形判断尚未完成，暂以不是多边形处理)
 * 参数：poly --- 多边形结构指针
 * 返回:
 *     -1 不是多边形
 *      0  给定点共线
 *      1  逆时针凸多边形
 *      2  顺时针凸多边形
 *      3  逆时针凹多边形
 *      4  顺时针凹多边形
 */
int CPolygon::IsPoly()
{
    enum TYPE {
        NOT = -1,               // 不是多边形
        LINE = 0,                // 共线
        ANTI_PROT = 1,      // 逆时针凸多边形
        DEASIL_PROT = 2,   // 顺时针凸多边形
        ANTI_CONC = 3,      // 逆时针凹多边形
        DEASIL_CONC = 4    // 顺时针凹多边形
        };

    int i, j=0, cnt=0, type = NOT;
    int left = 0, right = 0;
    int *dir = new int[this->count];

    if(this->count < 3)
        return NOT;
    else if(this->count == 3)
    {
        int pos = SideOfLine(acmes[0],acmes[1],acmes[2].x(),acmes[2].y());
        return pos == 0 ? LINE : (pos > 0 ? ANTI_PROT : DEASIL_PROT);
    }

    // 每次循环判断一条边
    for(i=0; i<this->count; ++i)
    {
        // 判断所有点在有向线段 P(i)P(i+1) 的哪侧，保存在 dir[] 中
        // 线段的端点不用判断，但为简化程序，不区分这两个点
        for(j=0; j<this->count; ++j)
            dir[j] = SideOfLine(this->acmes[i], this->acmes[i+1], this->acmes[j].x(), this->acmes[j].y());

        // 计算线段两侧点数
        for(left=0,right=0,j=0; j<this->count; ++j)
        {
            if(dir[j] > 0)
                ++left;
            else if(dir[j] < 0)
                ++right;
        }

        // 在一条直线上
        if(left == 0 && right == 0)
        {
            type = LINE;
            break;
        }
        // 分布在线段两侧，不是凸多形
        else if(left > 0 && right > 0)
        {
            // 假定是逆时针凹多边形
            type = ANTI_CONC;
            break;
        }
    }

    // 是凸多边形
    if(type == NOT)
    {
        if(left > 0 && right == 0) // 都在线段左边，逆时针
            type = ANTI_PROT;
        else if(left == 0 && right > 0) // 都在线段右边，顺时针
            type = DEASIL_PROT;
    }
    else if(type == ANTI_CONC)
    {
        // 判断是否凹多边形，尚未实现
        // 暂时作不是多边形处理
        type = NOT;
    }

    delete [] dir;

    return (int)type;
}


// Globals which should be set before calling this function:
//
// int    polySides  =  how many corners the polygon has
// float  polyX[]    =  horizontal coordinates of corners
// float  polyY[]    =  vertical coordinates of corners
// float  x, y       =  point to be tested
//
// (Globals are used in this example for purposes of speed.
// Change as desired.)
//
// The function will return TRUE if the point x,y is inside the
// polygon, or FALSE if it is not. If the point x,y is exactly on
// the edge of the polygon, then the function may return TRUE or
// FALSE.
//
// Note that division by zero is avoided because the division is
// protected by the "if" clause which surrounds it.

bool CPolygon::PtInPoly(double x, double y)
{
  int   i, j=0;
  bool  oddNODES=false ;

    /*CString str;
    str.Format("%d", count);
    AfxMessageBox(str);*/

  assert(count>=0);

  for (i=0; i<this->count; i++) {
    j++; if (j==this->count) j=0;
    if (this->acmes[i].y()<y && this->acmes[j].y()>=y
    ||  this->acmes[j].y()<y && this->acmes[i].y()>=y) {
      if (this->acmes[i].x()+(y-this->acmes[i].y())/(this->acmes[j].y()-this->acmes[i].y())*(this->acmes[j].x()-this->acmes[i].x())<x) {
        oddNODES=!oddNODES; }}}

  return oddNODES;

}

void CPolygon::Init(vector<Point> & ptsBoundary)
{
    this->count = (int)ptsBoundary.size();
    this->acmes = new Point[count+1];

    int i=0;
    // 顶点坐标
    for(vector<Point>::iterator pit = ptsBoundary.begin(); pit != ptsBoundary.end(); pit++ )
    {
        this->acmes[i] = *pit;
        i++;
    }

    // 最后一个顶点与第一个顶点重合
    this->acmes[i] = this->acmes[0];
}

double CPolygon::MaxInscribedCircle(Point& pt)
{
    if (count < 3) 	return -1;

    std::vector<Plane> planes; // # of planes should equal to adjVorVNum;
    std::vector<Point3D> interPts; //# of intersection pts should be less than C(adjVorVNum,3)
    double tmp = M_SQRT2_2;

    for(int i=0; i<count; i++)
    {
        Point vv1 = acmes[i];
        Point vv2 = acmes[(i+1)%count];

        Vector3D vec(vv2.x()-vv1.x(), vv2.y()-vv1.y(),0);
        Vector3D vecN(-vec.y, vec.x, 0);
        vecN.SetUnit();
        vecN= Vector3D(tmp*vecN.x, tmp*vecN.y, tmp);
        //vec.SetUnit();
        //Vector3D nor = vec^vecN;
        Vector3D nor(-vecN.x, -vecN.y, vecN.z);
        //RotatePointAroundVector(vecN, vec, nor);
        Vector3D pt(vv1.x(), vv1.y(), 0);
        double d =  nor*pt;
        Plane pl;
        pl.d = d;
        pl.normal = nor;
        planes.push_back(pl);
    }
    assert((int)planes.size()==count);

    for (int i =0; i<count-2; i++)
    {
        for (int j = i+1; j<count-1; j++)
        {
            for (int k =j+1; k<count; k++)
            {
                Point3D interpt;
                if (ThreePlanesInter(planes[i], planes[j], planes[k], interpt))
                {
                    interPts.push_back(interpt);
                }
            }
        }
    }

    std::vector<Point3D> interPts2; //# should be less than that of interPts
    for (int i =0; i<(int)interPts.size(); i++)
    {
        Vector3D interpt(interPts[i].x, interPts[i].y, interPts[i].z);
        bool isvalid = true;
        for (int j =0; j<count; j++)
        {
            double res = interpt*planes[j].normal - planes[j].d;
            if (res > TOLER)
            {//the point lies on the same side as the normal
                isvalid = false;
                break;
            }
        }
        if (isvalid)
        {
            interPts2.push_back(interPts[i]);
        }
    }

    double maxz = 0;
    Point3D maxzPt;
    for (int i =0; i<(int)interPts2.size(); i++)
    {
        if (abs(interPts2[i].z) > maxz)
        {
            maxzPt = interPts2[i];
            maxz = interPts2[i].z;
        }
    }

    pt = Point(maxzPt.x, maxzPt.y);

    double mind = 10e10;
    for(int i=0; i<count; i++)
    {
        Point vv1 = acmes[i];
        Point vv2 = acmes[(i+1)%count];

        double tmpx, tmpy;
        double dist = disPt2Line(vv1.x(), vv1.y(), vv2.x(), vv2.y(), pt.x(), pt.y(), tmpx, tmpy);

        if (dist < mind)
        {
            mind = dist;
        }
    }

    return mind;
}

std::pair<Point, Point> CPolygon::CommonEdge(const CPolygon& polyRight)
{
	int iLeftSize = this->count;
	int iRightSize = polyRight.count;
	//Point pointTmp1, pointTmp2;
	std::pair<Point, Point> pairLine;
	int bSetFirst = 0;
	for(int i = 0; i<iLeftSize; ++i)
	{
		//pointTmp1 = this->acmes[i];

		for(int j = 0; j<iRightSize; ++j)
		{
			//pointTmp2 = polyRight.acmes[j];
			if(this->acmes[i] == polyRight.acmes[j])
			{
				if(0 == bSetFirst)
				{
					pairLine.first = this->acmes[i];
					++bSetFirst;
				}
				else if(1 == bSetFirst)
				{
					pairLine.second = this->acmes[i];
					++bSetFirst;
				}
				else
				{
					std::cerr<< "More than 2 Points in Polygon is the same."<< endl;
				}
			}

		}//2for
	}//1for

	if(2 == bSetFirst)
	{
		return pairLine;
	}
	else
	{
		std::cerr<< "Less than 2 Points in Polygon is the same."<< endl;
	}
}

bool CPolygon::isCircleInPoly(const Point pointCircleCenter, double dCircleRadii)
{
	double dDistTmp = this->DistPt2Poly(pointCircleCenter);
	if(dDistTmp < dCircleRadii)
		return false;
	else
		return true;
	//if(PtInPoly(pointCircleCenter.x(), pointCircleCenter.y()))
	//{
	//	for(int i = 0; i<this->count; ++i)
	//	{
	//		assert(this->acmes[i%count] != this->acmes[(i+1)%count]);
	//		dDistTmp = std::abs((pointCircleCenter.x() - this->acmes[i%count].x())*
	//							(this->acmes[i%count].y()-this->acmes[(i+1)%count].y()) - 
	//							(pointCircleCenter.y() - this->acmes[i%count].y())*
	//							(this->acmes[i%count].x()-this->acmes[(i+1)%count].x()));
	//		dDistTmp /= std::sqrt((this->acmes[i%count].x()-this->acmes[(i+1)%count].x())*
	//							(this->acmes[i%count].x()-this->acmes[(i+1)%count].x()) +
	//							(this->acmes[i%count].y()-this->acmes[(i+1)%count].y())*
	//							(this->acmes[i%count].y()-this->acmes[(i+1)%count].y()));
	//		if(dDistTmp < dCircleRadii)
	//		{
	//			return false;
	//		}
	//	}
	//1for
	//	return true;
	//}//1if
	//else
	//{
	//	return false;
	//}
}

double CPolygon::DistPt2Poly(const Point pointTmp)
{
	double dDistTmp = 0.0;
	double dDistMin = 0.0;
	for(int i = 0; i<this->count; ++i)
	{
		assert(this->acmes[i%count] != this->acmes[(i+1)%count]);

		dDistTmp = std::abs((pointTmp.x() - this->acmes[i%count].x())*
							(this->acmes[i%count].y()-this->acmes[(i+1)%count].y()) - 
							(pointTmp.y() - this->acmes[i%count].y())*
							(this->acmes[i%count].x()-this->acmes[(i+1)%count].x()));
		dDistTmp /= std::sqrt((this->acmes[i%count].x()-this->acmes[(i+1)%count].x())*
							(this->acmes[i%count].x()-this->acmes[(i+1)%count].x()) +
							(this->acmes[i%count].y()-this->acmes[(i+1)%count].y())*
							(this->acmes[i%count].y()-this->acmes[(i+1)%count].y()));

		if((0 == i) || 
			(0 != i)&&(dDistMin>dDistTmp))
		{
			dDistMin = dDistTmp;
		}

	}//1for

	return dDistMin;
}

std::pair<Point, Point> CPolygon::Neighbour(const Point& pointMid)
{
	int iIndexTmp = -1;
	for(int i = 0; i<this->count; ++i)
	{
		if(pointMid == this->acmes[i])
		{
			iIndexTmp = i;
			break;
		}
	}

	std::pair<Point, Point> pairTmp;
	if(-1 == iIndexTmp)
	{
		std::cerr<< "No such a point"<< std::endl;
		return pairTmp;
	}

	Point pointTmp1, pointTmp2;
	pointTmp1 = this->acmes[(iIndexTmp-1+count)%count];
	pointTmp2 = this->acmes[(iIndexTmp+1)%count];
	pairTmp.first = pointTmp1;
	pairTmp.second = pointTmp2;
}

double CPolygon::TriEdgeCircleR(const Point point1, const Point point2, Point& pointCenter)
{
	std::pair<Point, Point> pairTmp1, pairTmp2;
	pairTmp1 = Neighbour(point1);
	pairTmp2 = Neighbour(point2);
	
	Point pointTmp1 = pairTmp1.first;
	if(pairTmp1.first == point2)
	{
		pointTmp1 = pairTmp1.second;
	}
	Point pointTmp2 = pairTmp2.second;
	if(pairTmp2.second == point1)
	{
		pointTmp2 = pairTmp2.first;
	}
	/*//0 situation
	//damn
	if(x5 == x4)
	{
		assert(y5 != y4 && x6 != x3);
		dYTmp = y3;
		dXTmp = x4 + (y3-y6)*(y3-y4)/(x6-x3);
	}
	else if(x6 == x3)
	{
		assert(y6 != y3 && x5 != x4);
		dYTmp = y4;
		dXTmp = x3 + (y4-y3)*(y4-y5)
	}
	else
	{
		dYTmp = ( (pointTmp1.x()-point2.x())*((pointTmp2.x()-point1.x())*point2.x() - (point1.y()-pointTmp2.y())*point2.y())-
			(pointTmp2.x()-point1.x())*((pointTmp1.x()-point2.x())*point1.x() - (point2.y()-pointTmp1.y())*point1.y()) )/
			( (pointTmp2.x()-point1.x())*(point2.y()-pointTmp1.y())-(pointTmp1.x()-point2.x())*(point1.y()-pointTmp2.y()) );
		dXTmp = ((point2.y()-pointTmp1.y())*dYTmp - (point2.y()-pointTmp1.y())*point1.y() + (pointTmp1.x()-point2.x())*point1.x())/
			(pointTmp1.x()-point2.x());
	}*/
	
	double dXTmp, dXTmp1, dYTmp1;
	double dYTmp, dXTmp2, dYTmp2;
	AngleSplit(pointTmp1.x(), pointTmp1.y(), point1.x(), point1.y(), 
		point2.x(), point2.x(), dXTmp1, dYTmp1);
	AngleSplit(point1.x(), point1.y(), point2.x(), point2.y(), 
		pointTmp2.x(), pointTmp2.y(), dXTmp2, dYTmp2);

	if(LinesCross(dXTmp1, dYTmp1, point1.x(), point1.y(), 
		dXTmp2, dYTmp2, point2.x(), point2.y(), dXTmp, dYTmp))
	{
		pointCenter = Point(dXTmp, dYTmp);
		return disPt2Line(point1.x(), point1.y(), point2.x(), point2.y(),
		pointCenter.x(), pointCenter.y(), dXTmp, dYTmp);
	}
	else
	{
		cerr<< "Not cross?!"<< endl;
		return 0.0;
	}
	
	
}