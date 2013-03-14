#include "BoundedPair.h"

bool BoundedPair::Bound(const int iTap1, const int iTap2, double dRate)
{
	if(iTap1 == iTap2)
	{
		return false;
	}
	int iTmp1 = iTap1;
	int iTmp2 = iTap2;
	if(iTap2 < iTap1)
	{
		iTmp1 = iTap2;
		iTmp2 = iTap1;
	}

	if(this->m_mapPairs.count(iTmp1))
	{
		return false;
	}
	else
	{
		this->m_mapPairs[iTmp1] = iTmp2;
		this->m_iSize++;

#ifdef _DEBUG
		if(this->m_mapTapRate.count(iTmp1))
			cout<<"Already got one?"<< endl;
		if(this->m_iSize != this->m_mapPairs.size())
			cout<< "Not equal size"<< endl;
#endif
		
		this->m_mapRadiInfo[iTmp1] = 0;
		if(dRate < 1)
		{
			this->m_mapTapRate[iTmp1] = (double)1.0/dRate;
			//this->m_mapTap2Rate[iTmp2] = (double)1.0/dRate;
		}
		else
		{
			this->m_mapTapRate[iTmp1] = dRate;
			//this->m_mapTap2Rate[iTmp1] = dRate;
		}
	}
	return true;
}

bool BoundedPair::isBounded(const int iTap1, const int iTap2)
{
	bool bResult = false;
	if(this->m_mapPairs.count(iTap1))
	{
		if(this->m_mapPairs[iTap1] == iTap2)
			bResult = true;
	}
	else if(this->m_mapPairs.count(iTap2))
	{
		if(this->m_mapPairs[iTap2] == iTap1)
			bResult = true;
	}

	return bResult;
}

double BoundedPair::RateOf(const int iTap1, const int iTap2)
{
	if(!isBounded(iTap1, iTap2))
	{
		return 0.0;
	}
	else
	{
		if(iTap1 > iTap2)
			return this->m_mapTapRate[iTap2];
		else if(iTap2 > iTap1)
			return this->m_mapTapRate[iTap1];
		else
			return 0.0;
	}
}

bool BoundedPair::BoundedWith(const int iTap, int& iTapTmp)
{//ASSUME: iTapTmp == 0; then iTap is the second
	iTapTmp = -1;
	if(this->m_mapPairs.count(iTap))
	{
		iTapTmp = this->m_mapPairs[iTap];
		return true;
	}
	else
	{
		map<int, int>::iterator iter = m_mapPairs.begin();
		for(; iter != this->m_mapPairs.end(); iter++)
		{
			if(iTap == iter->second)
			{
				iTapTmp = -1*iter->first;
				return true;
			}
		}
		return false;
	}
}

std::pair<int, double> BoundedPair::RadiFresh(const int iTap, const double dRadious)
{
	std::pair<int, double> pairTapRadii(-1, 1.0);

	if(this->m_mapPairs.count(iTap))
	{
		assert(this->m_mapRadiInfo.count(iTap) > 0);
		if(this->m_mapRadiInfo[iTap] <= 1e-12 && this->m_mapRadiInfo[iTap] >= -1e-12)
		{//==0
			this->m_mapRadiInfo[iTap] = dRadious;
			pairTapRadii.first = -1;
			pairTapRadii.second = -1.0;
		}
		else
		{
			pairTapRadii.first = this->m_mapPairs[iTap];
			pairTapRadii.second = this->m_mapRadiInfo[iTap];
		}
	}
	else
	{
		map<int, int>::iterator iter = m_mapPairs.begin();
		for(; iter != this->m_mapPairs.end(); iter++)
		{
			if(iTap == iter->second)
			{
				assert(this->m_mapRadiInfo.count(iter->first) > 0);
				if(this->m_mapRadiInfo[iter->first] <= 1e-12 &&
					this->m_mapRadiInfo[iter->first] >= -1e-12)
				{// == 0
					this->m_mapRadiInfo[iter->first] = -1 * dRadious;
					pairTapRadii.first = -1;
					pairTapRadii.second = -1.0;
				}
				else
				{
					pairTapRadii.first = iter->first;
					pairTapRadii.second = this->m_mapRadiInfo[iter->first];
				}
			}
		}
	}

	return pairTapRadii;
}

double BoundedPair::Radious(const int iTap1, const int iTap2, 
							const double dMicR1, const double dMicR2,
							cv::Point2d &point1, cv::Point2d &point2)//Mic center
{//
	double dDist = std::sqrt((point1.x-point2.x)*(point1.x-point2.x) + (point1.y-point2.y)*(point1.y-point2.y));
	if(dDist <= 1e-12 && dDist >= -1e-12)
	{
		return 0.0;
	}
	double dRTmp = std::abs(dMicR1);
	bool bXchange = false;
	double dRate = RateOf(iTap1, iTap2);

	int iTapTmp1 = iTap1, iTapTmp2 = iTap2;
	cv::Point2d vectorMN;
	cv::Point2d point0;
	if(std::abs(dMicR1) > std::abs(dMicR2))//set origin
	{
		dRTmp = std::abs(dMicR2);
		bXchange = true;
		iTapTmp1 = iTap2;
		iTapTmp2 = iTap1;
		vectorMN = cv::Point2d((point1.x-point2.x)/dDist, (point1.y-point2.y)/dDist);
		point0 = cv::Point2d((point2.x - vectorMN.x*dMicR2), (point2.y - vectorMN.y*dMicR2));
	}
	else
	{
		vectorMN = cv::Point2d((point2.x-point1.x)/dDist, (point2.y-point1.y)/dDist);
		point0 = cv::Point2d((point1.x - vectorMN.x*dMicR1), (point1.y - vectorMN.y*dMicR1));		
	}

	//check doable: dRate>2 ?
	
	assert(dRate >= 2.0);
	if((dDist < dRTmp*dRate) 
		&& ((dDist + std::abs(dMicR1-dMicR2)) > dRTmp*dRate))
	{
		double dRadiOld = dRTmp;
		dRTmp = (dDist - std::abs(dMicR1-dMicR2))/dRate;
		//circle heart refresh
		point1 = point0 + (2*dRadiOld-dRTmp)*vectorMN;
		point2 = point1 + dRate*dRTmp*vectorMN;		
	}
	else if((dDist > dRTmp*dRate) 
		&& ((dDist-std::abs(dMicR1-dMicR2)) > dRTmp*dRate))
	{//?
		double dRadiOld = dRTmp;
		dRTmp = (dDist - (dMicR1+dMicR2-dRTmp))/(dRate-1.0);
		//circle heart refresh
		point1 = point0 + (2*dRadiOld-dRTmp)*vectorMN;
		point2 = point1 + dRate*dRTmp*vectorMN;
	}
	else
	{
		//circle heart refresh
		point1 = point0 + dRTmp*vectorMN;
		point2 = point0 + (1+dRate)*dRTmp*vectorMN;
	}

	if(bXchange)
	{
		cv::Point2d pointTmp = point1;
		point1 = point2;
		point2 = pointTmp;
	}
#ifdef _DEBUG
//	if(dRTmp > dRadious1 || dRTmp > dRadious2)
//		cout<< "Error : recomputaion"<< endl;
#endif

	ClearMapRadi(iTap1, iTap2);

	return dRTmp;
}

double BoundedPair::Radious(const int iTap1, const int iTap2, const double dMicR1, const double dMicR2,
		cv::Point2d &point1, cv::Point2d &point2, CPolygon* poly1, CPolygon *poly2)
{//during the computation: point1 & dRadious1 => smaller circle
	//point2 & dRadious2 => larger circle
	//bXchange => swap or not
	bool bXchange1 = false;
	double dRadi = 0.0;
	double dDist = std::sqrt((point1.x-point2.x)*(point1.x-point2.x) + (point1.y-point2.y)*(point1.y-point2.y));
	double dRadi1 = dMicR1;//min of the given radious
	if(dMicR1 > dMicR2)
	{
		bXchange1 = true;
		dRadi1 = dMicR2;
	}
	
	double dRate = RateOf(iTap1, iTap2);
	Point point1Tmp, point2Tmp;
	std::pair<Point, Point> pairTmp = poly1->CommonEdge(*poly2);
	Point pointCenter1, pointCenter2;
	double dRTmp1 = poly1->TriEdgeCircleR(pairTmp.first, pairTmp.second, pointCenter1);
	double dRTmp2 = poly2->TriEdgeCircleR(pairTmp.first, pairTmp.second, pointCenter2);
	double dRadi2 = dRTmp1;//min of the new radi
	bool bXchange2 = false;
	if(dRTmp1 > dRTmp2)
	{
		bXchange2 = true;
		dRadi2 = dRTmp2;
	}
	cv::Point2d pointTmp1, pointTmp2;
	bool bXchange = false;
	if(dRadi1 > dRadi2)
	{//3LineCircle is small
		dRadi = dRadi2;
		dDist = std::sqrt((pointCenter1.x()-pointCenter2.x())*(pointCenter1.x()-pointCenter2.x()) +
			(pointCenter1.y()-pointCenter2.y())*(pointCenter1.y()-pointCenter2.y()));
		if(bXchange2)
		{
			bXchange = true;
			//pointTmp1 = cv::Point2d(pointCenter2.x(), pointCenter2.y());
			//pointTmp2 = cv::Point2d(pointCenter1.x(), pointCenter1.y());
			point1 = cv::Point2d(pointCenter2.x(), pointCenter2.y());
			point2 = cv::Point2d(pointCenter1.x(), pointCenter1.y());
			dRadi1 += dRadi2;
			dRadi2 = dRadi1 - dRadi2;
			dRadi1 = dRadi1 - dRadi2;
		}
		else
		{
			point1 = cv::Point2d(pointCenter1.x(), pointCenter1.y());
			point2 = cv::Point2d(pointCenter2.x(), pointCenter2.y());
			pointTmp1 = cv::Point2d(pointCenter1.x(), pointCenter1.y());
			pointTmp2 = cv::Point2d(pointCenter2.x(), pointCenter2.y());
		}
	}
	else
	{//MIC is small
		dRadi = dRadi1;
		if(bXchange1)
		{
			bXchange = true;
			cv::Point2d poinTmp = point1;
			//pointTmp1 = point2;
			//pointTmp2 = point1;
			point1 = point2;
			point2 = poinTmp;
			//
			dRadi1 = dMicR2;
			dRadi2 = dMicR1;
		}
		else
		{
			//pointTmp1 = point1;
			//pointTmp2 = point2;
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
	cv::Point2d vectorMN = cv::Point2d((point2.x-point1.x)/dDist, (point2.y-point1.y)/dDist);
	if( LinesCross(pairTmp.first.x(), pairTmp.first.y(), pairTmp.second.x(), pairTmp.second.y(),
		pointTmp1.x, pointTmp1.y, pointTmp2.x, pointTmp2.y, dXTmp, dYTmp) )
	{
		if((dDist < dRadi*dRate) && (dDist + std::abs(dRadi2 - dRadi1)) < dRadi*dRate)
		{//not enought space. dRadi = (dDist + dRadious2 - dRadious1)/dRate;
				//dDist = dDist + dRadious2 - dRadiou1
				dRadi = (dDist + std::abs(dRadi2 - dRadi1))/dRate;
		}
		//point1 = crossPoint - 0.5*dRadi*dRate*MN
		//point2 = crossPoint + 0.5*dRadi*dRate*MN
		point1.x = dXTmp - vectorMN.x*dRadi*dRate*0.5;
		point1.y = dYTmp - vectorMN.y*dRadi*dRate*0.5;
		point2.x = dXTmp + vectorMN.x*dRadi*dRate*0.5;
		point2.y = dYTmp + vectorMN.y*dRadi*dRate*0.5;		
	}
	else
	{
		std::cerr<< "No cross point"<< std::endl;
	}
	//if(dRate == 2)
	//{
	//	if((dDist - dRadious1 - dRadious2) > 0 )
	//	{// got a gap between the 2 circles
	//		std::pair<Point, Point> pairTmp = poly1->CommonEdge(*poly2);

	//	}
	//}
	//else
	//{

	//}
	if(bXchange)
	{
		cv::Point2d poinTmp = point1;
		//pointTmp1 = point2;
		//pointTmp2 = point1;
		point1 = point2;
		point2 = poinTmp;
	}

	ClearMapRadi(iTap1, iTap2);

	return dRadi;
}

bool BoundedPair::ClearMapRadi(const int iTap1, const int iTap2)
{
	if(this->m_mapRadiInfo.count(iTap1))
	{
		this->m_mapRadiInfo[iTap1] = 0.0;
		return true;
	}
	else if(this->m_mapRadiInfo.count(iTap2))
	{
		this->m_mapRadiInfo[iTap2] = 0.0;
		return true;
	}
	else
	{//not bounded
		return false;
	}
}

bool BoundedPair::Draw(cv::Mat& matCanvas, vector<cv::Point2d>& aPoints)
{//aPoints' point order should be unchangeable;
	if(matCanvas.empty() || 0 == aPoints.size())
	{
		return false;
	}
	else
	{
		map<int, int>::iterator iter = this->m_mapPairs.begin();
		for(; iter != this->m_mapPairs.end();)
		{
			cv::line(matCanvas, aPoints[iter->first], aPoints[iter->second], cv::Scalar::all(255), 2);
			++iter;
		}
		
		return true;
	}
}