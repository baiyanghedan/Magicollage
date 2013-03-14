//jerrry
//For circle between 3 edge(4 points)

#include "circle_packing\CCVTVertex.h"

bool CommonEdge(const CCVTVertex vertex1, const CCVTVertex vertex2, 
				Point& point1, Point& point2)
{
	int iSize1 = vertex1.adjVorVertex.GetSize();
	int iSize2 = vertex2.adjVorVertex.GetSize();
	VorVertex * pVertexTmp1 = NULL;
	VorVertex * pVertexTmp2 = NULL;
	bool bFirstSet = false;
	for(int i = 0; i<iSize1; ++i)
	{
		pVertexTmp1 = (VorVertex *)vertex1.adjVorVertex.GetAt(i);
		//pVertexTmp1->id;
		for(int j = 0; j<iSize2; ++j)
		{
			pVertexTmp2 = (VorVertex *)vertex2.adjVorVertex.GetAt(j);
			if(pVertexTmp1->pos == pVertexTmp2->pos)
			{
				if(!bFirstSet)
				{
					point1 = pVertexTmp1->pos;
					bFirstSet = true;
				}
				else
				{
					point2 = pVertexTmp1->pos;
					return true;
				}
			}
		}
	}
	return false;
}

void CCVTVertex::Neighbour(const Point pointMid, Point& pointLeft, Point& pointRite) const
{
	VorVertex * pVertexTmp = NULL;
	int iSize = adjVorVertex.GetSize();
	for(int i = 0; i<iSize; ++i)
	{
		pVertexTmp = (VorVertex*)adjVorVertex.GetAt(i);
		if(pVertexTmp->pos == pointMid)
		{
			pVertexTmp = (VorVertex*)adjVorVertex.GetAt((i-1+iSize)%iSize);
			pointLeft =  pVertexTmp->pos;
			pVertexTmp = (VorVertex*)adjVorVertex.GetAt((i+1)%iSize);
			pointRite = pVertexTmp->pos;
			break;
		}
	}
}

double Radious(const int iTap1, const int iTap2, const double dMicR1, const double dMicR2, const double dRate,
		Point &point1, Point &point2, const CCVTVertex vertex1, const CCVTVertex vertex2)
{//during the computation: point1 & dRadious1 => smaller circle
	//point2 & dRadious2 => larger circle
	//bXchange => swap or not
	//RESULT: iTap1 - point1  iTap2 - point2
	bool bXchange1 = false;
	double dRadi = 0.0;
	double dDist = std::sqrt((point1.x-point2.x)*(point1.x-point2.x) + (point1.y-point2.y)*(point1.y-point2.y));
	double dRadi1 = dMicR1;//min of the given radious
	if(dMicR1 > dMicR2)
	{
		bXchange1 = true;
		dRadi1 = dMicR2;
	}
	
	Point pointCommon1, pointCommon2;
	if(CommonEdge(vertex1, vertex2, pointCommon1, pointCommon2))
	{
		Point pointTmp1, pointTmp2, pointTmp3, pointTmp4;
		Point pointCenter1, pointCenter2;
		//vertex1 cell
		vertex1.Neighbour(pointCommon1, pointTmp1, pointTmp2);
		vertex1.Neighbour(pointCommon2, pointTmp3, pointTmp4);
		if(point2 == pointTmp1)
		{
			pointTmp1 = pointTmp2;
		}
		else if(point2 != pointTmp2)
		{//error neighbour 
			std::cerr<< "Erro neighbour"<< std::endl;
			return -1;
		}

		if(pointTmp3 == point1)
		{
			pointTmp2 = pointTmp4;
		}
		else if(pointTmp4 == point1)
		{
			pointTmp2 = pointTmp3;
		}
		else
		{//error neighbour 
			std::cerr<< "Erro neighbour"<< std::endl;
			return -1;
		}

		double dRTmp1 = TriEdgeCircleR(pointCommon1, pointCommon2,
			pointTmp1, pointTmp2, pointCenter1);
		//vertex2 cell
		vertex2.Neighbour(pointCommon1, pointTmp1, pointTmp2);
		vertex2.Neighbour(pointCommon2, pointTmp3, pointTmp4);
		if(point2 == pointTmp1)
		{
			pointTmp1 = pointTmp2;
		}
		else if(point2 != pointTmp2)
		{//error neighbour 
			std::cerr<< "Erro neighbour"<< std::endl;
			return -1;
		}

		if(pointTmp3 == point1)
		{
			pointTmp2 = pointTmp4;
		}
		else if(pointTmp4 == point1)
		{
			pointTmp2 = pointTmp3;
		}
		else
		{//error neighbour 
			std::cerr<< "Erro neighbour"<< std::endl;
			return -1;
		}

		double dRTmp2 = TriEdgeCircleR(pointCommon1, pointCommon2,
			pointTmp1, pointTmp2, pointCenter2);

		double dRadi2 = dRTmp1;//min of the new radi
		bool bXchange2 = false;
		if(dRTmp1 > dRTmp2)
		{
			bXchange2 = true;
			dRadi2 = dRTmp2;
		}
		//cv::Point2d pointTmp1, pointTmp2;
		bool bXchange = false;
		if(dRadi1 > dRadi2)
		{//3LineCircle is small
			dRadi = dRadi2;
			dDist = std::sqrt((pointCenter1.x()-pointCenter2.x())*(pointCenter1.x()-pointCenter2.x()) +
				(pointCenter1.y()-pointCenter2.y())*(pointCenter1.y()-pointCenter2.y()));
			if(bXchange2)
			{
				bXchange = true;
				point1 = pointCenter2;
				point2 = pointCenter1;
				dRadi1 += dRadi2;
				dRadi2 = dRadi1 - dRadi2;
				dRadi1 = dRadi1 - dRadi2;
			}
			else
			{
				point1 = pointCenter1;
				point2 = pointCenter2;
				//pointTmp1 = pointCenter1;
				//pointTmp2 = pointCenter2;
			}
		}
		else
		{//MIC is small
			dRadi = dRadi1;
			if(bXchange1)
			{
				bXchange = true;
				Point poinTmp = point1;
				point1 = point2;
				point2 = poinTmp;
				//
				dRadi1 = dMicR2;
				dRadi2 = dMicR1;
			}
			else
			{
				dRadi1 = dMicR1;
				dRadi2 = dMicR2;
			}
		}
		//Radi & circle center computation
		//dRadi1 ~ point1 | small M
		//dRadi2 ~ point2 | big   N

		//cross point(dXTmp, dYTmp)
		double dXTmp, dYTmp;
		//vectorMN : point1 -> point2
		double vectorMNx = (point2.x-point1.x)/dDist;
		double vectorMNy = (point2.y-point1.y)/dDist;
		if( LinesCross(pointCommon1.x(), pointCommon1.y(), pointCommon2.x(), pointCommon2.y(),
			point1.x(), point1.y(), point2.x(), point2.y(), dXTmp, dYTmp) )
		{
			if((dDist < dRadi*dRate) && (dDist + std::abs(dRadi2 - dRadi1)) < dRadi*dRate)
			{//not enought space. dRadi = (dDist + dRadious2 - dRadious1)/dRate;
				//dDist = dDist + dRadious2 - dRadiou1
				dRadi = (dDist + std::abs(dRadi2 - dRadi1))/dRate;
			}
			//point1 = crossPoint - 0.5*dRadi*dRate*MN
			//point2 = crossPoint + 0.5*dRadi*dRate*MN
			point1 = Point(dXTmp - vectorMNx*dRadi*dRate*0.5, dYTmp - vectorMNy*dRadi*dRate*0.5);
			point2 = Point(dXTmp + vectorMNx*dRadi*dRate*0.5, dYTmp + vectorMNy*dRadi*dRate*0.5);
		}
		else
		{
			std::cerr<< "No cross point"<< std::endl;
		}
		
		//Result 
		if(bXchange)
		{
			Point poinTmp = point1;
			//pointTmp1 = point2;
			//pointTmp2 = point1;
			point1 = point2;
			point2 = poinTmp;
		}

		//ClearMapRadi(iTap1, iTap2);

		return dRadi;
	}
	else
	{
		return -1;
	}
}