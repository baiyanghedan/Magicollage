//#include "../stdafx.h"
#include <algorithm>
#include "CPCC2D.h"
#include "UtilFunc.h"
#include "ContCircle2.h"
#include "CCVTVertex.h"

#include <tnt.h>
//#include <tnt_linalg_addon.h>
//#include <tnt_stopwatch.h>

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

//#include <Windows.h>

#define CPK 1

CPCC2D::CPCC2D(void)
{
    m_DiagramType = PDIAGRAM;
    m_nIter = 0;
    pLastPD = NULL;
    vMinCPRate = NULL;
    vMaxCPRate = NULL;
    poly = NULL;
	
	////////////////////////////////////////////
	//jerrry
	//////////////////////////////////////////////
	pPairInfo = NULL;
    _CrtDumpMemoryLeaks();
}

CPCC2D::CPCC2D(std::vector<WeightedPoint>& pts, BoundType boundtype)
{
    m_DiagramType = PDIAGRAM;
    m_nIter = 0;
    pLastPD = NULL;
    vMinCPRate = NULL;
    vMaxCPRate = NULL;
    poly = NULL;
	////////////////////////////////////////////
	//jerrry
	//////////////////////////////////////////////
	pPairInfo = NULL;

    m_boundType = boundtype;
    for(int i=0; i<(int)pts.size(); ++i)
        AddPoint(pts[i]);	
}
///////////////////////////////////////////////////
//jerrry
//////////////////////////////////////////////
CPCC2D::CPCC2D(std::vector<WeightedPoint>& pts, BoundType boundtype, BoundedPair& pairInfo)
{
    m_DiagramType = PDIAGRAM;
    m_nIter = 0;
    pLastPD = NULL;
    vMinCPRate = NULL;
    vMaxCPRate = NULL;
    poly = NULL;
	pPairInfo = &pairInfo;

	m_boundType = boundtype;
    for(int i=0; i<(int)pts.size(); ++i)
        AddPoint(pts[i]);	
}

CPCC2D::~CPCC2D(void)
{
    if (pLastPD != NULL)  delete pLastPD;
    if (poly != NULL)  delete poly;
	///////////////////////////////////////////////////////////////
	//jerrry
	///////////////////////////////////////////////////////////////
	if (pPairInfo != NULL) delete pPairInfo;

    FreeMemory();	
}

Vertex_handle CPCC2D::AddPoint(WeightedPoint& wpt)
{
    CLloydBase2D::AddPoint(wpt);
    Vertex_handle vh;
    m_nIter = 0;
    if (m_DiagramType == PDIAGRAM)
    {		
        vh = m_pdt.insert(wpt);
        AddCVTVertex(vh);
        Update();
    }else if(m_DiagramType == MDIAGRAM)
	{
        AddCVTVertex(wpt);
    }
    return vh;
}

Vertex_handle CPCC2D::AddRestrictedPoint(WeightedPoint& wpt)
{
    CLloydBase2D::AddPoint(wpt);
    Vertex_handle vh;
    m_nIter = 0;
    if (m_DiagramType == PDIAGRAM)
    {		
        vh = m_pdt.insert(wpt);
        AddCVTVertex(vh, true);
        
    }else if(m_DiagramType == MDIAGRAM){
        AddCVTVertex(wpt);
    }
    return vh;
}

void CPCC2D::ShuffleSeeds()
{
    if(!m_bCirclePacking) return;

    std::vector<double> weights;
    for(int i=0; i<(int)m_points.size(); ++i){
        weights.push_back(m_points[i].weight());
    }
    random_shuffle( weights.begin(), weights.end() );
    for(int i=0; i<(int)m_points.size(); ++i){
        m_points[i] = WeightedPoint(m_points[i].point(), weights[i]);
    }
    FreeMemory();
    Vertex_handle vh;
    for(int i=0; i<(int)m_points.size(); ++i){
        vh = m_pdt.insert(m_points[i]);
        AddCVTVertex(vh);
    }
    m_nIter = 0;
}

bool comp_greater(double a, double b)
{
    return a>b;
};

void CPCC2D::GenerateSeedsInitDS1(const int num)
{	//init according to the "DS1" refer to Hifi08 	[for CIRCLE bound only]
    if((num <1)||(!m_bSetRadii) ||(!m_bCirclePacking) ||(m_boundType != CIRCLEBOUND))
    {
        return;
    }	
    assert(m_cpSize == num);
    if(m_cpSize < 1) return;

    m_pdt.clear();
    m_points.clear();
    FreeMemory();

    std::sort(m_radii.begin(), m_radii.end(), comp_greater);
    
    Vector2D* initpts = new Vector2D[num];

    double cx=m_radii[1];	double cy=0;
    double cr = m_radii[0]+m_radii[1];
    initpts[0] = Vector2D(0, 0);
    initpts[1] = Vector2D(m_radii[0]+m_radii[1], 0);

    std::vector<Vector2D> possibleS;

    //******
    /*CString strFileR;
    strFileR.Format("c:\\debug_r.txt");
    std::ofstream	fpr(strFileR,std::ios::out);*/

    for (int i=2; i < num; i++)
    {
        possibleS.clear();

        double resx[2]={0};
        double resy[2]={0};
        
        double r = m_radii[i];
        double mindelta = 10e10;

        /*fpr << i+1 <<"\n";*/
        //step1 
        for(int j=0; j<i-1; j++){
            for(int l=j+1; l<i; l++){
                bool res = tangentto2Circles(initpts[j].x, initpts[j].y, m_radii[j], initpts[l].x, initpts[l].y, m_radii[l], r, resx, resy);
                if (!res)
                {
                    continue;
                }
                //make sure that no overlapping	
                for (int m = 0; m<2; m++)
                {
                    bool overlapped = false;
                    for (int k = 0; k<i; k++)
                    {
                        if((k!=j)&&(k!=l)&& (is2CirclesOverlapped(resx[m],resy[m],r,initpts[k].x, initpts[k].y, m_radii[k])))
                        {
                            overlapped = true;
                            break;
                        }
                    }
                    if (!overlapped)
                    {
                        possibleS.push_back(Vector2D(resx[m], resy[m]));
                    }
                }
            }
        }//finished filling possibleS
        assert(possibleS.size()>0);

        //step2
        double minrc = 10e10;
        Vector2D mincenter, minpos;

        int currNum = i+1;
        Vector2D* currPts = new Vector2D[currNum];
        for (int j=0; j<i; j++)
        {
            currPts[j] = initpts[j];
        }
        int tmpK =0;
        for(int j = 0; j<(int)possibleS.size(); j++)
        {
            currPts[i] = possibleS[j];
            Circle2 minCir;
            FindMinContainingCircle(i+1, currPts, minCir);
            /*fpr<< std::scientific << minCir.Radius<<"\n";*/
            if (minCir.Radius + TOLER < minrc)
            {
                minrc = minCir.Radius;
                mincenter = minCir.Center;
                minpos = currPts[i];
                tmpK = j;
            }
        }
        /*fpr<<minrc << "  #" <<tmpK << "\n\n";*/
        delete[] currPts;

        initpts[i] = minpos;
        cr = minrc;
        cx = mincenter.x;
        cy = mincenter.y;

        ////print to a file to testify the result 
        //CString strFile;
        //strFile.Format("c:\\debug_%d.txt", currNum);
        //std::ofstream	fp(strFile,std::ios::out);	
        //assert(fp.good());
        //double tmpk = radius_r/cr;
        //fp << 2 << ' ' << currNum << ' ' << m_boundType << std::endl; 
        //for(int j=0; j<currNum; j++)
        //{
        //	fp << std::scientific << (double)(tmpk*(initpts[j].x-cx)) << ' ' << (double)(tmpk*(initpts[j].y-cy))<<' ' << (tmpk*tmpk*m_radii[j]*m_radii[j]) <<"\n";
        //}
        //fp.close();
    }

    /*fpr.close();*/
    
    //reconstruct the bounding circle
    double k = radius_r/cr;
    
    m_k = k/**0.5*/;

    /*std::ofstream mout("c:\\test.txt");
    mout<<"2  "<<num<<"  1\n";*/
    for (int i=0; i<num; i++)
    {
        //mout<<k*(initpts[i].x-cx) << "    "<<k*(initpts[i].y-cy)<<"   "<<m_k*m_radii[i]*m_k*m_radii[i]<<"\n";
        AddPoint(WeightedPoint(Point(k*(initpts[i].x-cx), k*(initpts[i].y-cy)), m_radii[i]*m_radii[i]));		
    }
    //mout.close();

    delete[] initpts;	
}

bool CPCC2D::FindMinContainingCircle(const int num, const Vector2D* initpts,  Circle2& rkMinimal)
{
    //Using CGAL to find the closest bounding circle
    std::vector<Sphere> S;
    for (int i=0; i<num; i++)
    {
        S.push_back(Sphere(Point_2(initpts[i].x, initpts[i].y), m_radii[i]));
    }
    Min_sphere ms(S.begin(),S.end());
    
    double cgalc[2] = {0};
    int d = 0;
    for (Min_sphere::Cartesian_const_iterator c_it = ms.center_cartesian_begin();
        c_it != ms.center_cartesian_end(); 	++c_it){
            cgalc[d] = *c_it;
            d++;
    }	
    rkMinimal.Radius = ms.radius();
    rkMinimal.Center = Vector2D(cgalc[0], cgalc[1]);
    return (ms.is_valid());
}

void CPCC2D::GenerateSeedsInitDS2(const int num)
{	//init according to the "DS2" refer to Hifi08  [for CIRCLE bound only]
    if((num <1)||(!m_bSetRadii) ||(!m_bCirclePacking) ||(m_boundType != CIRCLEBOUND))
    {
        return;
    }

    assert(m_cpSize == num);
    if(m_cpSize < 1) return;

    m_pdt.clear();
    m_points.clear();
    FreeMemory();

    std::sort(m_radii.begin(), m_radii.end(), comp_greater);

    double cr = 0;
    for (int i =0; i<num; i++)
    {
        cr += m_radii[i]*m_radii[i];
    }
    cr = sqrt(cr);
    double cx=0;	double cy=0;

    Vector2D* initpts = new Vector2D[num];

    initpts[0] = Vector2D(0, 0);
    initpts[1] = Vector2D(m_radii[0]+m_radii[1], 0);
    
    for (int i=2; i < num; i++)
    {
        double resx[2]={0};
        double resy[2]={0};
        double delta=0, x, y, currx, curry;
        double r = m_radii[i];
        bool bfound = false;
        bool bsucc = false;
        double mindelta = 10e10;
        for(int j=0; j<i-1; j++){
            if(bfound) break;
            for(int l=j+1; l<i; l++){
                bool res = tangentto2Circles(initpts[j].x, initpts[j].y, m_radii[j], initpts[l].x, initpts[l].y, m_radii[l], r, resx, resy);
                if (!res)
                {
                    continue;
                }
                //make sure that no overlapping
                bool succtangent = false;
                for (int m = 0; m<2; m++)
                {
                    bool overlapped = false;
                    for (int k = 0; k<i; k++)
                    {
                        if ((k!=j)&&(k!=l)&& (is2CirclesOverlapped(resx[m],resy[m],r,initpts[k].x, initpts[k].y, m_radii[k])))
                        {
                            overlapped = true;						
                            break;
                        }
                    }
                    if (overlapped)
                    {
                        continue;
                    }else{
                        x = resx[m];
                        y = resy[m];
                        succtangent = true;
                        bsucc = true;
                        break;
                    }
                }
                if (!succtangent)
                {
                    continue;
                }
                if (isCircleInCircle(x,y,r,cx,cy,cr,delta))
                {
                    initpts[i] = Vector2D(x,y);
                    bfound = true;
                    break;
                }else{
                    if (delta < mindelta)
                    {
                        mindelta = delta;
                        currx = x;
                        curry = y;
                    }
                }
            }
        }

        if (!bsucc)
        {
            //CString str;
            //str.Format(_T("Error occurs in DS2, didn't find the appropriate position for circle %d"), i);
            //AfxMessageBox(str);
            //printf("%s\n", str);
            printf("Error occurs in DS2, didn't find the appropriate position for circle %d", i);
        }
        if(!bfound){
            initpts[i] = Vector2D(currx,curry);
            cr += mindelta;
        }
    }

    Circle2 minCir;
    FindMinContainingCircle(num, initpts, minCir);
    double newr = minCir.Radius;
    double newcx = minCir.Center.x;
    double newcy = minCir.Center.y;

    ////find the closest bounding circle only using all centers (From WM4)
    //Circle2 minc;
    //MinCircle2* mincircle = new MinCircle2(num, initpts, minc);
    //double newr = -10e10;
    //for (int i=0; i<num; i++)
    //{
    //	double dist = SqrtDistance(initpts[i], minc.Center) + m_radii[i];
    //	if (dist > newr)
    //	{
    //		newr = dist;
    //	}
    //}
    //assert(newr > minc.Radius);	

    //double newcx = minc.Center.x;
    //double newcy = minc.Center.y;	
    //delete mincircle;
     
    //reconstruct the bounding circle
    double k = radius_r/newr;
    m_k = k/**0.01*/;
    for (int i=0; i<num; i++)
    {
        AddPoint(WeightedPoint(Point(k*(initpts[i].x-newcx), k*(initpts[i].y-newcy)), m_radii[i]*m_radii[i]));		
    }

    delete[] initpts;
}

void CPCC2D::CVT_MBH(const int num, double h, std::string strFile, int& iternum)
{//Monotonic basin hopping
    if (num <1) return;

    bool bwritefile = true;
    if (strFile == ""){
        bwritefile = false;
    }
    
    h = geth();

    bool bimproved = false;
    int tmpnum = 0;
    const int MaxNoImp = 2000;
    int failNum = 0;
    iternum = 0;
    Vector2D* ptcoords = new Vector2D[m_cpSize];	

    for (unsigned int j=0; j<m_cpSize; j++)
    {
        ptcoords[j] = Vector2D(m_points[j].point().x(), m_points[j].point().y());
    }

    double minE = EvaluateErrorFunc();

    std::ofstream	fp(strFile, std::ios::out);	
    if (bwritefile)
    {
        fp << 2 << ' ' << (int)m_points.size()<< ' ' << m_boundType  << "  MBH  "<< m_initMethod <<"  total:"<<num<<std::endl; 
        fp <<"h = "<< h <<"   start energy: "<<minE<< std::endl;
    }		

    TNT::Stopwatch mwatchCP;
    mwatchCP.start();	
    
    for (int i=0; i<num; i++)
    {
        DoLloydAll(tmpnum);
        double currentE = EvaluateErrorFunc();
        fp <<tmpnum<<"  ";
        iternum += tmpnum;

        bimproved = false;
        
        if (currentE < minE)
        {
            minE = currentE;
            for (unsigned int j=0; j<m_cpSize; j++)
            {
                ptcoords[j] = Vector2D(m_points[j].point().x(), m_points[j].point().y());
            }
            bimproved = true;
            failNum = 0;
        }else{
            for (unsigned int j=0; j<m_cpSize; j++)
            {
                m_points[j] = WeightedPoint(Point(ptcoords[j].x, ptcoords[j].y), m_points[j].weight());
            }
            failNum ++;
        }
        mwatchCP.stop();
        if (bwritefile)
        {
            fp << std::setprecision(16) << std::scientific << currentE;
            if (bimproved)
            {
                fp << "    Improved";
            }else{
                fp << "    not improved";
            }
            fp<<"    time: "<< mwatchCP.read();
        }
        fp << std::endl; 
        mwatchCP.resume();
        if (failNum > MaxNoImp)
        {
            /*h *= 2;
            fp <<"h increases!    h = "<< h << std::endl;*/
            fp << "  stop with MaxNoImp = "<< MaxNoImp <<std::endl;
            break;
        }

        if (i<num-1){
            Perturbation(h);
            Update();
            //fp << "         perturbed";
        }
    }

    mwatchCP.stop();

    Update();	

    if (bwritefile){		
        fp <<"Number of tests is:  "<< num <<"\n";
        fp <<"h = "<< h <<"\n";
        fp <<"The running time is:  "<< mwatchCP.read() <<"\n";
        fp <<"The min Energy is:  "<< std::setprecision(16) << std::scientific << minE <<"\n";		
        fp.close();
    }

    delete[] ptcoords;
}

void CPCC2D::CP_A10or15Global(const int num, std::string strFile, bool useBFGS, int& iternum, InitMethod type)
{
    if (num <1) return;
    
    std::ofstream	fp(strFile, std::ios::out);	

    bool testTime = false;
    if (!testTime)
    {
        fp << 2 << ' ' << (int)m_points.size()<< ' ' << m_boundType  ;
        if (type==A10)
        {
            fp<<"  A10 times"<<std::endl; 
        }else if (type==A15)
        {
            fp<<"  A15 times"<<std::endl; 
        }
    }

    double delta = 10e10;
    int tmpnum, totalnum = 0;
    iternum = 0;	

    if (!testTime)
    {
        fp <<"The start k is:  "<<std::setprecision(16) << std::scientific << m_k <<"\n";
    }

    TNT::Stopwatch mwatchCP, mwatchA;
    mwatchCP.start();   mwatchCP.stop();
    mwatchA.start(); 	mwatchA.stop();
    while(delta > 1e-8 && totalnum < num)
    {
        double oldAlpha = m_k;
        
        mwatchCP.resume();

        mwatchA.resume();
        if (type==A10)
        {
            GenerateSeedsInitA10New(m_cpSize);
        }else if (type==A15)
        {
            GenerateSeedsInitA15(m_cpSize);
        }
        mwatchA.stop();

        CirclePackingAll(tmpnum, true);
        
        mwatchCP.stop();
        iternum += tmpnum;

        if (!testTime)
        {
            fp <<tmpnum;
            fp <<"  The max k:  "<< std::setprecision(16) <<std::scientific << m_k;
            fp <<"  The min radius:  "<< std::setprecision(16) << std::scientific << 1.0/m_k ;
            fp<<" A10/15 t: "<<mwatchA.read();
            fp <<"  Time:"<<mwatchCP.read()<<"\n";
        }

        delta = abs(m_k-oldAlpha);
        totalnum++;
    }	

    fp <<"A10/15 Time used:  "<<mwatchA.read()<<"\n";
    fp <<"The running time is:  "<< mwatchCP.read() <<"\n";
    fp <<"Number of tests is:  "<< totalnum <<"\n";	
    fp.close();	
}

void CPCC2D::CPEqualCircleTests(const int num, const std::string&  strFile, int& iternum)
{
    if (num <1) return;

    bool bwritefile = true;
    if (strFile == "")
    {
        bwritefile = false;
    }
    
    std::ofstream	fp(strFile.c_str(), std::ios::out);	

    if (bwritefile)
    {
        fp << 2 << ' ' << (int)m_points.size()<< ' ' << m_boundType  << "  random tests\n";
    }

    double maxalpha = -10e10;
    Vector2D* ptcoords = new Vector2D[m_cpSize];
    bool found = false;
    iternum = 0;

    TNT::Stopwatch mwatchCP;
    mwatchCP.start();

    for (int i=0; i<num; i++)
    {
        GenerateRandomSeeds(m_cpSize);
        int tmpnum = 0;
        DoLloydAll(tmpnum);
        
        CirclePackingAll(tmpnum, false);
        iternum += tmpnum;

        if (m_k > maxalpha)
        {
            maxalpha = m_k;
            for (unsigned int j=0; j<m_cpSize; j++)
            {
                ptcoords[j] = Vector2D(m_points[j].point().x(), m_points[j].point().y());
            }
            found = true;
        }
        if (bwritefile)
        {
            fp << std::setprecision(16) <<std::scientific << m_k <<"\n";
        }
    }

    mwatchCP.stop();

    if (bwritefile)
    {
        fp <<"The running time is:  "<< mwatchCP.read() <<"\n";
        fp <<"Number of tests is:  "<< num <<"\n";
        fp <<"The max k is:  "<< std::setprecision(16) <<std::scientific << maxalpha <<"\n";
        if (m_boundType == CIRCLEBOUND)
        {
            fp <<"The min radius is:  "<< std::setprecision(16) << std::scientific << 1.0/m_k <<"\n";
        }
        fp.close();
    }

    if (!found){
        delete[] ptcoords;
        return;
    }
    m_pdt.clear();
    m_points.clear();
    FreeMemory();

    m_k = maxalpha;
    for (unsigned int i=0; i<m_cpSize; i++)
    {
        AddPoint(WeightedPoint(Point(ptcoords[i].x, ptcoords[i].y), m_radii[i]*m_radii[i]));		
    }
    delete[] ptcoords;
}

double CPCC2D::geth(){
    //double h_ = 1.0/sqrt((double)m_cpSize)*0.5;  //h2
    //double h_ = 1.0/sqrt((double)m_cpSize)*0.8;  //h8
    //double h_ = 1.0/sqrt((double)m_cpSize)*0.9;  //h9
    //double h_ = pow(1.0/sqrt(PI*m_cpSize)*0.5,(0.5/(double)m_cpSize)) *0.5*((exp(1.0)/(double)m_cpSize));  //h
    //double h_ = 1.0/sqrt((double)m_cpSize)*0.25;   //h4
    //h_ *= 0.5;
    //if (m_boundType == RECTBOUND){
    //	//h_ *= 1.2;
    //	h_ *= (bound_x_max-bound_x_min);
    //}else{
    //	h_ *= 2;
    //}
    //return h_;
    return 0.5;
}

//compute the dissimilarity of two local minimizers
void CPCC2D::dissimilarity_measure(std::vector<double> & delta)
{
    
    Vector2D bcenter = Vector2D(0,0);
    for (unsigned int j=0; j<m_cpSize; j++){
        bcenter += Vector2D(m_points[j].point().x(), m_points[j].point().y());
    }
    bcenter /= m_cpSize;

    for (unsigned int j=0; j<m_cpSize; j++){
        double dis = (bcenter.x - m_points[j].point().x())*(bcenter.x - m_points[j].point().x()) + (bcenter.y - m_points[j].point().y())*(bcenter.y - m_points[j].point().y());
        delta.push_back(sqrt(dis));
    }
    std::sort(delta.begin(), delta.end());
}

double get_dm(std::vector<double> & delta1, std::vector<double> & delta2){
    double res =0;
    assert(delta1.size()==delta2.size());
    for (unsigned int i=0; i<delta1.size(); i++){
        res += fabs(delta1[i]-delta2[i]);
    }
    return res;
}

int get_min_k_id(std::vector<double> & ks){
    double mink = 10e10;
    int minid = -1;
    for (unsigned int i=0; i<ks.size(); i++){
        if (ks[i] < mink){
            mink = ks[i];
            minid = i;
        }
    }
    return minid;
}

int get_max_k_id(std::vector<double> & ks, double & maxk){
    maxk = -10e10;
    int maxid = -1;
    for (unsigned int i=0; i<ks.size(); i++){
        if (ks[i] > maxk){
            maxk = ks[i];
            maxid = i;
        }
    }
    return maxid;
}

int CPCC2D::CP_MBH(const int num, double h, bool lloydfirst, bool bperturb, std::string strFile, bool useBFGS, int& iternum)
{//Monotonic basin hopping
    if (num <1) return 0;

    bool bwritefile = true;
    if (strFile == ""){
        bwritefile = false;
    }
    
    //h = geth();
    
    m_k = 1e-6;
    double ferr;
    doCirclePackingIteration(ferr, CPK);
    double max_k = m_k;//-10e10;

    std::ofstream	fp(strFile, std::ios::out);	

    if (bwritefile)
    {
        fp << 2 << ' ' << (int)m_points.size()<< ' ' << m_boundType  << "  MBH  "<< m_initMethod <<"  total:"<<num<<std::endl; 
        fp <<"h = "<< h << "   Lloyd first: "<< lloydfirst << "   use BFGS: "<< useBFGS <<"   start k: "<<m_k<< std::endl;
    }

    bool bimproved = false;
    const int MaxNoImp = 200;
    int tmpnum = 0;
    int ImpNum = 0;
    int failNum = 0;
    iternum = 0;
    Vector2D* ptcoords = new Vector2D[m_cpSize];

    for (unsigned int j=0; j<m_cpSize; j++)
    {
        ptcoords[j] = Vector2D(m_points[j].point().x(), m_points[j].point().y());
    }

    TNT::Stopwatch mwatchCP;
    mwatchCP.start();
    
    if (lloydfirst){
        DoLloydAll(tmpnum);
        fp<<"Lloyd time: "<<mwatchCP.read() << std::endl;
    }
    for (int i=0; i<num; i++)
    {
        CirclePackingAll(tmpnum, false);
        
        fp <<tmpnum<<"  ";
        iternum += tmpnum;

        bimproved = false;
        
        if (m_k > max_k)
        {
            max_k = m_k;
            for (unsigned int j=0; j<m_cpSize; j++)
            {
                ptcoords[j] = Vector2D(m_points[j].point().x(), m_points[j].point().y());
            }
            bimproved = true;
            failNum = 0;
            ImpNum ++;
        }else{
            for (unsigned int j=0; j<m_cpSize; j++)
            {
                m_points[j] = WeightedPoint(Point(ptcoords[j].x, ptcoords[j].y), m_points[j].weight());
            }
            failNum ++;
        }
        mwatchCP.stop();
        if (bwritefile)
        {
            fp << std::setprecision(14) << std::scientific << max_k<<"   "<<m_k;
            if (bimproved)
            {
                fp << "    Improved";
            }else{
                fp << "    not improved";
            }
            fp<<"    time: "<< mwatchCP.read();
        }
        fp << std::endl; 
        mwatchCP.resume();
        if (failNum > MaxNoImp)
        {
            /*h *= 2;
            fp <<"h increases!    h = "<< h << std::endl;*/
            fp << "  stop with MaxNoImp = "<< MaxNoImp <<std::endl;
            break;
        }

        if (i<num-1){
            if (bperturb){
                Perturbation(h);
            }else{
                ShuffleSeeds();
            }
            //m_k *= 0.6;
            Update();
            //fp << "         perturbed";
        }
    }

    mwatchCP.stop();

    m_k=max_k;
    Update();
    
    iternum += tmpnum;

    if (bwritefile)
    {
        fp <<"last #iter:  "<< tmpnum <<"\n";
        fp <<"Number of tests is:  "<< num <<"\n";
        fp <<"h = "<< h <<"\n";
        fp <<"The running time is:  "<< mwatchCP.read() <<"\n";
        fp <<"The max k is:   "<< std::setprecision(16) << std::scientific << max_k <<"\n";
        fp <<"The last k is:  "<< std::setprecision(16) << std::scientific << m_k <<"\n";
        if (m_boundType == CIRCLEBOUND)
        {
            fp <<"The min radius is:  "<< std::setprecision(16) << std::scientific << 1.0/m_k <<"\n";
        }
        fp.close();
    }

    delete[] ptcoords;
    return ImpNum;
}

void CPCC2D::CP_SA(const int num, double h, bool lloydfirst, std::string strFile, bool useBFGS, int& iternum)
{//SimulatedAnnealing
    if (num <1) return;

    bool bwritefile = true;
    if (strFile == ""){
        bwritefile = false;
    }
    
    std::ofstream	fp(strFile, std::ios::out);	
    double T = 0.5; //the smaller T is, the smaller probability SA accepts an smaller k

    h = geth();

    if (bwritefile)
    {
        fp << 2 << ' ' << (int)m_points.size()<< ' ' << m_boundType  <<"  SA"<< std::endl; 
        fp <<"h = "<< h << "    T0 = "<<T<< "   Lloyd first: "<< lloydfirst << "   use BFGS: "<< useBFGS << std::endl;
    }
    
    double lastalpha = -10e10;
    double bestalpha = -10e10;
    
    Vector2D* ptcoords = new Vector2D[m_cpSize];
    Vector2D* bestptcoords = new Vector2D[m_cpSize];

    bool bsa,bimproved = false;
    iternum = 0;
    int tmpnum = 0;
    int MaxNoImp = 30;
    int failNum = 0;

    TNT::Stopwatch mwatchCP;
    mwatchCP.start();
    
    if (lloydfirst){
        DoLloydAll(tmpnum);
        fp<<"Lloyd time: "<<mwatchCP.read() << std::endl;
    }

    for (int i=0; i<num; i++)
    {		
        CirclePackingAll(tmpnum, false);
        
        fp <<tmpnum<<"  ";		
        iternum += tmpnum;

        bimproved = false; bsa = false;
        double delta  = exp((m_k-lastalpha)/T);
        T *= 0.9;

        if (m_k > bestalpha){
            bestalpha = m_k;
            for (unsigned int j=0; j<m_cpSize; j++)
            {
                bestptcoords[j] = Vector2D(m_points[j].point().x(), m_points[j].point().y());
            }
        }
        if ((m_k > lastalpha) ||(delta> (rand()*1.0/RAND_MAX)))
        {
            if (m_k < lastalpha){
                bsa = true;	
            }
            lastalpha = m_k;
            for (unsigned int j=0; j<m_cpSize; j++)
            {
                ptcoords[j] = Vector2D(m_points[j].point().x(), m_points[j].point().y());
            }
            bimproved = true;
            failNum = 0;
        }else{
            for (unsigned int j=0; j<m_cpSize; j++)
            {
                m_points[j] = WeightedPoint(Point(ptcoords[j].x, ptcoords[j].y), m_points[j].weight());
            }
            failNum ++;
        }
        if (bwritefile)
        {
            fp << std::setprecision(16) <<std::scientific << m_k;
            if (bimproved)
            {
                if (bsa){
                    fp << "    SA  delta:" <<delta;
                }else
                    fp << "     improved";
            }else{
                fp << "     not improved";
            }
            fp<<"    time: "<<mwatchCP.read();
        }

        if (i<num-1){			
            Perturbation(h);
            m_k *= 0.6;
            Update();
            //fp << "         perturbed";
        }
        fp << std::endl; 

        /*if (failNum > MaxNoImp)
        {
            h *= 2;
            fp <<"h increases!    h = "<< h << std::endl;
        }*/
    }	
    mwatchCP.stop();

    if (m_k < bestalpha){
        for (unsigned int j=0; j<m_cpSize; j++)
        {
            m_points[j] = WeightedPoint(Point(bestptcoords[j].x, bestptcoords[j].y), m_points[j].weight());
        }
    }

    Update();
    CirclePackingAll(tmpnum, false);
    iternum += tmpnum;

    if (bwritefile)
    {
        fp <<"Number of tests is:  "<< num <<"\n";
        fp <<"h = "<< h <<"\n";
        fp <<"The running time is:  "<< mwatchCP.read() <<"\n";
        fp <<"The max k is:  "<< std::setprecision(16) <<std::scientific << m_k <<"\n";
        if (m_boundType == CIRCLEBOUND)
        {
            fp <<"The min radius is:  "<< std::setprecision(16) << std::scientific << 1.0/m_k <<"\n";
        }
        fp.close();
    }

    delete[] ptcoords;
    delete[] bestptcoords;
}


void CPCC2D::CPRandomSeedsTimes(const int num, std::string strFile, int& iternum)
{
    if (num <1) return;
    std::cout.precision(16);
    std::ofstream	fp(strFile, std::ios::out);	
    
    fp << 2 << ' ' << (int)m_points.size()<< ' ' << m_boundType  << std::endl; 

    double maxalpha = -10e10;
    Vector2D* ptcoords = new Vector2D[m_cpSize];
    bool found = false;
    iternum = 0;

    TNT::Stopwatch mwatchCP;
    mwatchCP.start();

    for (int i=0; i<num; i++)
    {
        GenerateRandomSeeds(m_cpSize);
        int tmpnum = 0;
        CirclePackingAll(tmpnum, true);
        iternum += tmpnum;

        if (m_k > maxalpha)
        {
            maxalpha = m_k;
            for (unsigned int j=0; j<m_cpSize; j++)
            {
                ptcoords[j] = Vector2D(m_points[j].point().x(), m_points[j].point().y());
            }
            found = true;
        }
        fp << std::scientific << m_k <<"\n";
    }
    
    mwatchCP.stop();

    fp <<"The running time is:  "<< mwatchCP.read() <<"\n";
    fp <<"Number of tests is:  "<< num <<"\n";
    fp <<"The max k is:  "<< std::scientific << maxalpha <<"\n";
    if (m_boundType == CIRCLEBOUND)
    {
        fp <<"The min radius is:  "<< std::setprecision(16) << std::scientific << 1.0/m_k <<"\n";
    }
    fp.close();

    if (!found){
        delete[] ptcoords;
        return;
    }
    m_pdt.clear();
    m_points.clear();
    FreeMemory();

    m_k = maxalpha;
    for (unsigned int i=0; i<m_cpSize; i++)
    {
        AddPoint(WeightedPoint(Point(ptcoords[i].x, ptcoords[i].y), m_radii[i]*m_radii[i]));		
    }
    delete[] ptcoords;
}

void CPCC2D::GenerateBoundarySeeds(const int num)
{
    /*double theta = PI2/num;
    for (int i = 0; i < num; i++)
    {
        double x = center_x + radius_r*cos(theta*i);
        double y = center_y + radius_r*sin(theta*i);
        AddRestrictedPoint(WeightedPoint(Point(x, y), 0.0001));
        m_radii.push_back(0.01);
    }
    m_cpSize = (int)m_radii.size();*/
    
    int tmpnum = (int)ptsBoundary.size();	
    int count = -1;
    for (int i = 0; i < tmpnum; i++)
    {
        double x = ptsBoundary[i].x();
        double y = ptsBoundary[i].y();
        double x1 = ptsBoundary[(i+1)%tmpnum].x();
        double y1 = ptsBoundary[(i+1)%tmpnum].y();
        double dist = (x-x1)*(x-x1) + (y-y1)*(y-y1);
        if (dist < 0.005/*0.02*/)
        {
            count ++;
            if (count%3 != 0)
            {
                continue;
            }			
        }else if (dist > 2.6)//add points on the two longest edges of the boundary
        {
            double deltax = abs(x - x1)/7;
            double deltay = abs(y - y1)/7;
            double tmpx, tmpy;
            for (int j=1; j<7; j++)
            {
                if (x> x1)
                {
                    tmpx = x - deltax*j;
                }else{
                    tmpx = x + deltax*j;
                }
                if (y> y1)
                {
                    tmpy = y - deltay*j;
                }else{
                    tmpy = y + deltay*j;
                }
                AddRestrictedPoint(WeightedPoint(Point(tmpx, tmpy), 0.0001));
                m_radii.push_back(0.01);
            }
            continue;
        }

        AddRestrictedPoint(WeightedPoint(Point(x, y), 0.0001));
        m_radii.push_back(0.01);
    }
    m_cpSize = (int)m_radii.size();
}

//--------------------------------------------------------
//modified by magiccatq
void CPCC2D::GenerateRandomSeeds(const int num)
{
    if (num <1)  return;
    
    m_pdt.clear();
    m_points.clear();
    FreeMemory();
    if (!m_bSetRadii)
    {
        m_radii.clear();
    }

    double coord[2];	
    if(this->m_aOuterInitPos.size() == num)
    {
        for (int i = 0; i < num; i++)
        {//printf("use outer\n");	
            coord[0] = this->m_aOuterInitPos[i].x();
            coord[1] = this->m_aOuterInitPos[i].y();
            while(!isInBound(coord[0],coord[1]))
            {//printf("use rand\n");
                for (int j = 0; j < 2; j++) 
                {
                    coord[j] = rand()*2.0/RAND_MAX - 1.0;
                } 
            }

            if (m_bCirclePacking)
            {
                if (m_bSetRadii)
                {
                    AddPoint(WeightedPoint(Point(coord[0], coord[1]), m_radii[i]));
                }else{
                    AddPoint(WeightedPoint(Point(coord[0], coord[1]), 0.0001));
                    m_radii.push_back(0.01);
                }
                m_k = INIT_K;
            } 
            else
            {
                AddPoint(WeightedPoint(Point(coord[0], coord[1]), 0.0));
            }		
        }
    }
    else
    {
        for (int i = 0; i < num; i++)
        {		
            do{
                for (int j = 0; j < 2; j++) 
                {
                    coord[j] = rand()*2.0/RAND_MAX - 1.0;
                } 
            }while(!isInBound(coord[0],coord[1]));

            if (m_bCirclePacking)
            {
                if (m_bSetRadii)
                {
                    AddPoint(WeightedPoint(Point(coord[0], coord[1]), m_radii[i]*m_radii[i]));
                }else{
                    AddPoint(WeightedPoint(Point(coord[0], coord[1]), 0.0001));
                    m_radii.push_back(0.01);
                }
                m_k = INIT_K;
            } 
            else
            {
                AddPoint(WeightedPoint(Point(coord[0], coord[1]), 0.0));
            }		
        }
    }//else
    m_bSetRadii = true;
    m_cpSize = (int)m_radii.size();	
}
//--------------------------------------------------------
/*
void CPCC2D::GenerateRandomSeeds(const int num)
{
    if (num <1)  return;
    
    m_pdt.clear();
    m_points.clear();
    FreeMemory();
    if (!m_bSetRadii)
    {
        m_radii.clear();
    }

    double coord[2];	
    
    for (int i = 0; i < num; i++)
    {		
        do{
            for (int j = 0; j < 2; j++) 
            {
                coord[j] = rand()*2.0/RAND_MAX - 1.0;
            } 
        }while(!isInBound(coord[0],coord[1]));

        if (m_bCirclePacking)
        {
            if (m_bSetRadii)
            {
                AddPoint(WeightedPoint(Point(coord[0], coord[1]), m_radii[i]*m_radii[i]));
            }else{
                AddPoint(WeightedPoint(Point(coord[0], coord[1]), 0.0001));
                m_radii.push_back(0.01);
            }
            m_k = INIT_K;
        } 
        else
        {
            AddPoint(WeightedPoint(Point(coord[0], coord[1]), 0.0));
        }		
    }
    m_bSetRadii = true;
    m_cpSize = (int)m_radii.size();	
}
*/

void CPCC2D::GenerateSeedsinGrid(const int rows, const int columns)
{
    if (m_boundType != RECTBOUND)
    {
        return;
    }
    const int num = rows * columns;
    if (num <=0)
    {
        return;
    }	
    const double height = bound_y_max - bound_y_min;
    const double width = bound_x_max - bound_x_min;
    
    m_pdt.clear();
    m_points.clear();
    FreeMemory();
    m_radii.clear();	

    double cellwidth = width/columns;
    double cellheight = height/rows;
    double coord[2];
    for (int i=0; i<rows; i++)
    {
        double tmpy = bound_y_max - i*cellheight - cellheight*0.5;
        for (int j=0; j<columns; j++)
        {
            coord[0] = bound_x_min + j*cellwidth + cellwidth*0.5;
            coord[1] = tmpy;

            if (m_bCirclePacking)
            {
                AddPoint(WeightedPoint(Point(coord[0], coord[1]), 0.0001));
                m_radii.push_back(0.01);
            }else
            {
                AddPoint(WeightedPoint(Point(coord[0], coord[1]), 0.0));
            }
        }
    }

    if (m_bCirclePacking)
    {
        m_bSetRadii = true;
        m_k =  INIT_K;
        m_cpSize = (int)m_radii.size();
    }
}

void CPCC2D::Update()
{//vertices pos
    int nv = (int)m_pdt.number_of_vertices();
    if( nv < 2)	return;	
    
    m_pdt.clear();

    double k2 = m_k*m_k;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        if (m_bCirclePacking)
        {//vit: CVT region info
            double w = k2*m_points[(*vit)->id].weight();
            (*vit)->vh  = m_pdt.insert(WeightedPoint(m_points[(*vit)->id].point(), w));
            (*vit)->weight = w;
        } 
        else
        {
            (*vit)->vh  = m_pdt.insert(m_points[(*vit)->id]);
            (*vit)->weight = m_points[(*vit)->id].weight();
        }		
        (*vit)->footPrints[m_nIter] = (*vit)->pos;
        (*vit)->pos = m_points[(*vit)->id].point();
        (*vit)->type = VALID;
        (*vit)->isInsideCell = true;
        (*vit)->ClearNeighbor();
    }//point weight info refresh
    
	//memory release
    FreeCVTFaces();
    FreeVorVertices();

    if (m_pdt.number_of_hidden_vertices() > 0)//?
    {
        for(Hidden_verticex_iterator vit = m_pdt.hidden_vertices_begin(); vit != m_pdt.hidden_vertices_end(); vit++)
        {
            CCVTVertex * vtx = FindMyVertex(vit);
            vtx->type = HIDDEN;
        }
    }

    if(nv == 2)
    {		
        if ((m_boundType == RECTBOUND) ||(m_boundType == POLYBOUND))
        {
            InitTwoVertex_RectPoly();
        }
		else if(m_boundType == CIRCLEBOUND)
        {
            InitTwoVertex_Circle();
        }
        ComputeCentroids();
        return;
    }
    
	//leave them alone
    CreateCVTFaces();		// create faces and 
    if (m_boundType == RECTBOUND)
    {		
        CreateVorVertices_Rect();	// create corresponding voronoi vertex
    }else if(m_boundType == CIRCLEBOUND)
    {		
        CreateVorVertices_Circle();	// create corresponding voronoi vertex
    }else if (m_boundType == POLYBOUND)
    {
        CreateVorVertices_Poly();
    }

    ComputeCentroids();
}

void CPCC2D::AddCVTVertex(Vertex_handle vh, bool restricted)
{
    CCVTVertex * v = new CCVTVertex;
    v->id = (int)vertices.size();
    v->vh = vh;
    v->pos = vh->point().point();
    v->weight = vh->point().weight();
    v->diagramType = PDIAGRAM;
    v->isRestricted = restricted;
    vertices.push_back(v);
}

void CPCC2D::AddCVTVertex(WeightedPoint pt)
{
    CCVTVertex * v = new CCVTVertex;
    v->id = (int)vertices.size();
    v->vh = NULL;
    v->pos = pt.point();
    v->weight = pt.weight();
    v->diagramType = MDIAGRAM;
    vertices.push_back(v);
}

//not used actually
void CPCC2D::CreateCVTVertices()
{
    // create CVT vertex
    int index = 0;
    for(Finite_vertices_iterator vit = m_pdt.finite_vertices_begin(); vit != m_pdt.finite_vertices_end(); vit++)
    {
        CCVTVertex * v = new CCVTVertex;
        v->id = index;
        v->vh = vit;
        Vertex_handle vh = vit;
        v->pos = vh->point().point();
        v->weight = vh->point().weight();
        vertices.push_back(v);
        index ++;
    }
}

bool CPCC2D::isInBound(double x, double y)
{
    bool inside;
    if (m_boundType == RECTBOUND)
    {
        inside = isInBox(x,y);
    }else if(m_boundType == CIRCLEBOUND)
    {
        inside = isInCircle(x,y);
    }else if(m_boundType == POLYBOUND)
    {
        inside = PtInPoly(x,y);
    }
    return inside;
}

void CPCC2D::CreateCVTFaces()
{
    if( (int)m_pdt.number_of_faces() <= 0 ) return;
	//face info refresh
    int index = 0;
	////jerrry
	//std::vector<std::vector<cv::Point2d>> aTriangles;
	//std::vector<cv::Point2d> aTriangle(3);
    for(Face_iterator fit = m_pdt.faces_begin(); fit != m_pdt.faces_end(); fit++)
    {
        CCVTFace * f = new CCVTFace;
        f->id = index ++;
        f->fh = fit;
        faces.push_back(f);

        Face fd = * fit;
		////jerrry
		//aTriangle.clear();
        for( int i=0; i<3; i++)
        {
            Vertex_handle vh = fd.vertex(i);
            for(std::list<CCVTVertex*>::iterator vit=vertices.begin(); vit != vertices.end(); vit++)
            {
                CCVTVertex * v = *vit; 
                if(vh == v->vh){
                    assert(v->type != HIDDEN);
                    f->v[i] = v;
					////jerrry
					//aTriangle.push_back(cv::Point2d(v->pos.x(), v->pos.y()));
                }
            }
        }
		////jerrry
		//aTriangles.push_back(aTriangle);

        // create voronoi vertex which are the circum center of triangle in the bounding box
        Point centerpt = m_pdt.weighted_circumcenter(fit);
        if( isInBound(centerpt.x(), centerpt.y()) ) {
            if (! FindVorVertex(centerpt) ){
                VorVertex * vv = new VorVertex;
                vv->type = CIRCUMCENTER;
                vv->pos = centerpt;
                vv->face = f;
                vorvertices.push_back(vv);
            }
        }
    }
	/*cv::Mat matTmp = cv::Mat::zeros((cv::Size(1050, 750)), CV_8UC3);
	cv::imshow("Blank", matTmp);
	cv::waitKey();
	cv::drawContours(matTmp, aTriangles, -1, cv::Scalar(0, 255, 0), 2);
	cv::imshow("CVTFaces", matTmp);
	cv::waitKey();*/
}

void CPCC2D::CreateVorVertices_Poly()
{//	
    if (poly == NULL) return;
    
    int index = 0;

    // add corner vertex of bounding box to the voronoi vertex list
    for(int i=0; i<poly->count; i++)
	{
        AddVoronoiVertex(CORNER, poly->acmes[i]); 
    }

    int nocellNum = 0;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = (CCVTVertex*)*vit;

        v->isAdjBorder = false;
        v->isTangent2Border = false;
        v->isActive = true;

        if (v->type == HIDDEN)//?
        {
            v->ClearNeighbor();
            continue;
        }

        VorVertex			*vorv;
        Edge_circulator		ec = m_pdt.incident_edges(v->vh), done = ec;
        Vertex_circulator	vc = m_pdt.incident_vertices(v->vh), donevc = vc;
        
        Segment				s, edge;
        Ray					r;
        Line				l;
        bool				bIntersect = false;	// whether an one ring edge of a vertex intersect with boundary

        CGAL::Object		result;
        Point				ppos;

        // add one ring vertex to current vertex's adjacent vertex list.
        do {
            if (m_pdt.is_infinite(vc) ) 
			{
                vc ++;
                continue;
            }

            CCVTVertex * vtx = FindMyVertex(vc);
            assert(vtx != NULL);
            v->addAdjVertex(vtx);

            vc ++;
        } while ( vc != donevc );

        // create one ring voronoi vertex for each vertex	
        do {
            if( m_pdt.is_infinite(ec) ){
                ec ++;
                continue;
            }

            edge = m_pdt.segment(ec);	// edge
            index = -1;

            CGAL::Object o = m_pdt.dual(ec);

            if (CGAL::assign(s,o)) 
			{
                if (s.is_degenerate())
                {
                    int tmp = 2;
                }
                //judge whether the two weighted pts lie on the same side of the segment
                bool b1 = IsTriangleCCW(s.target(), edge.target(), s.source());
                bool b2 = IsTriangleCCW(s.target(), edge.source(), s.source());
                if ((b1&&b2) || !(b1||b2))
                {//b1 and b2 have the same signal, i.e., on the same side
                    CCVTVertex * v1 = FindMyVertex(edge.target());
                    CCVTVertex * v2 = FindMyVertex(edge.source());
                    assert(v1 != NULL);
                    assert(v2 != NULL);
                    if(v1->weight != v2->weight)
                    {					
                        if (v1->weight < v2->weight)
                        {
                            v1->isInsideCell = false;
                        }else
						{
                            v2->isInsideCell = false;
                        }
                    }
                }

                if(isInBound(s.source().x(),s.source().y()) && isInBound(s.target().x(),s.target().y()))
                {
                    vorv = FindVorVertex(s.source());
                    assert(vorv != NULL);
                    v->addVorVertex(vorv);
                    vorv = FindVorVertex(s.target());
                    assert(vorv != NULL);
                    v->addVorVertex(vorv);
                }
                else if(isInBound(s.target().x(),s.target().y()))
				{
                    for (int i=0; i<poly->count; i++)
                    {
                        Segment bdedge(poly->acmes[i], poly->acmes[i+1]);
                        result = CGAL::intersection(s,bdedge);

                        if (CGAL::assign(ppos, result))
                        {
                            // add a new voronoi vertex	
                            vorv = FindVorVertex(ppos);
                            if( !vorv )	{
                                vorv = new VorVertex;
                                vorv->type = INTERSECT;
                                vorv->pos  = ppos;
                                vorv->v[0] = FindMyVertex(edge.target());
                                vorv->v[1] = FindMyVertex(edge.source());
                                vorv->bindex = i;
                                vorvertices.push_back(vorv);
                            }
                            assert(vorv != NULL);
                            if (vorv->type == INTERSECT)
                            {
                                v->isAdjBorder = true;
                            }
                            v->addVorVertex(vorv);
                            vorv = FindVorVertex(s.target());
                            assert(vorv != NULL);
                            v->addVorVertex(vorv);
                            bIntersect = true;
                        }
                    }
                }
                else if(isInBound(s.source().x(),s.source().y()))
                {
                    for (int i=0; i<poly->count;i++)
                    {
                        Segment bdedge(poly->acmes[i], poly->acmes[i+1]);
                        result = CGAL::intersection(s,bdedge);

                        if (CGAL::assign(ppos, result))
                        {
                            // add a new voronoi vertex	
                            vorv = FindVorVertex(ppos);
                            if( vorv != NULL ){
                                assert(vorv != NULL);
                                double x = vorv->pos.x();
                            }
                            else if( vorv == NULL )
                            {
                                vorv = new VorVertex;
                                vorv->type = INTERSECT;
                                vorv->pos  = ppos;
                                vorv->v[0] = FindMyVertex(edge.target());
                                vorv->v[1] = FindMyVertex(edge.source());
                                vorv->bindex= i;
                                vorvertices.push_back(vorv);
                                assert(vorv != NULL);
                                double x = vorv->pos.x();
                            }
                            assert(vorv != NULL);
                            if (vorv->type == INTERSECT)
                            {
                                v->isAdjBorder = true;
                            }
                            v->addVorVertex(vorv);
                            vorv = FindVorVertex(s.source());
                            assert(vorv != NULL);
                            v->addVorVertex(vorv);
                            bIntersect = true;
                        }
                    }
                }
                else 
				{	// the start point the ray lie outside the bounding box. so has two intersection voronoi points
                    Point2D linesp(s.source().x(),s.source().y()), pa, pb, pts;
                    Vector2D dir(s.target().x() - s.source().x(), s.target().y()-s.source().y());

                    for(int i=0; i<poly->count; i++)
                    {
                        pa = Point2D(poly->acmes[i].x(),poly->acmes[i].y());
                        pb = Point2D(poly->acmes[i+1].x(),poly->acmes[i+1].y());

                        int n = LinesegLineInter(linesp, dir, pa, pb, pts);
                        if ( n < 1 ) continue;
                        else	{
                            Point pos(pts.x, pts.y);
                            double t1 = (pts.x - s.source().x())/(s.target().x() - s.source().x());
                            double t2 = (pts.y - s.source().y())/(s.target().y() - s.source().y());
                            if( t1 < 0 || t1 > 1 || t2 < 0 || t2 > 1) continue;
                            vorv = FindVorVertex(pos);
                            if( !vorv )	{
                                vorv = new VorVertex;
                                vorv->type = INTERSECT;
                                vorv->pos  = pos;
                                vorv->v[0] = FindMyVertex(edge.target());
                                vorv->v[1] = FindMyVertex(edge.source());
                                vorv->bindex= i;								
                                vorvertices.push_back(vorv);
                            }
                            assert(vorv != NULL);
                            if (vorv->type == INTERSECT)
                            {
                                v->isAdjBorder = true;
                            }
                            v->addVorVertex(vorv);
                            bIntersect = true;
                        }
                    }
                }
            }
            else if (CGAL::assign(r,o)) 
			{
                if (r.is_degenerate())
                {
                    ec++;
                    continue;
                }else
				{
                    //judge whether the two weighted pts lie on the same side of the segment
                    Point pt_on_r = Point(r.source().x()+r.direction().dx(),r.source().y()+r.direction().dy());
                    bool b1 = IsTriangleCCW(pt_on_r, edge.target(), r.source());
                    bool b2 = IsTriangleCCW(pt_on_r, edge.source(), r.source());
                    if ((b1&&b2) || !(b1||b2))
                    {//b1 and b2 have the same signal, i.e., on the same side
                        CCVTVertex * v1 = FindMyVertex(edge.target());
                        CCVTVertex * v2 = FindMyVertex(edge.source());
                        assert(v1 != NULL);
                        assert(v2 != NULL);
                        if(v1->weight != v2->weight)
                        {					
                            if (v1->weight < v2->weight)
                            {
                                v1->isInsideCell = false;
                            }else{
                                v2->isInsideCell = false;
                            }
                        }
                    }

                    if(isInBound(r.source().x(),r.source().y()))
                    {
                        Point2D linesp(r.source().x(),r.source().y()), pa, pb, pts;
                        Vector2D dir(r.direction().dx(),r.direction().dy());

                        for (int i=0; i<poly->count; i++)
                        {
                            pa = Point2D(poly->acmes[i].x(), poly->acmes[i].y());
                            pb = Point2D(poly->acmes[i+1].x(), poly->acmes[i+1].y());

                            int n = LinesegLineInter(linesp, dir, pa, pb, pts);
                            if ( n < 1 ) continue;
                            else{
                                if(dir * Vector2D(pts - linesp) <= 0)	continue;

                                Point pos(pts.x, pts.y);

                            /*Segment bdedge(poly->acmes[i],poly->acmes[i+1]);
                            result = CGAL::intersection(r,bdedge);

                            if (CGAL::assign(ppos,result))
                            {*/
                                // add a new voronoi vertex	
                                vorv = FindVorVertex(pos);
                                if( !vorv )	{
                                    vorv = new VorVertex;
                                    vorv->type = INTERSECT;
                                    vorv->pos  = pos;
                                    vorv->v[0] = FindMyVertex(edge.target());
                                    vorv->v[1] = FindMyVertex(edge.source());
                                    vorv->bindex= i;
                                    vorvertices.push_back(vorv);
                                }
                                assert(vorv != NULL);
                                if (vorv->type == INTERSECT)
                                {
                                    v->isAdjBorder = true;
                                }
                                v->addVorVertex(vorv);
                                vorv = FindVorVertex(r.source());
                                assert(vorv != NULL);
                                v->addVorVertex(vorv);
                                bIntersect = true;
                            }
                        }
                    }
                    else {	// the start point the ray lie outside the bounding box. so has two intersection voronoi points
                        Point2D linesp(r.source().x(),r.source().y()), pa, pb, pts;
                        Vector2D dir(r.direction().dx(),r.direction().dy());

                        for(int i=0; i<poly->count; i++)
                        {
                            pa = Point2D(poly->acmes[i].x(),poly->acmes[i].y());
                            pb = Point2D(poly->acmes[i+1].x(),poly->acmes[i+1].y());

                            int n = LinesegLineInter(linesp, dir, pa, pb, pts);
                            if ( n < 1 ) continue;
                            else{
                                if(dir * Vector2D(pts - linesp) <= 0)	continue;

                                Point pos(pts.x, pts.y);
                                vorv = FindVorVertex(pos);
                                if( !vorv )	{
                                    vorv = new VorVertex;
                                    vorv->type = INTERSECT;
                                    vorv->pos  = pos;
                                    vorv->v[0] = FindMyVertex(edge.target());
                                    vorv->v[1] = FindMyVertex(edge.source());
                                    vorv->bindex= i;
                                    vorvertices.push_back(vorv);
                                }
                                assert(vorv != NULL);
                                if (vorv->type == INTERSECT)
                                {
                                    v->isAdjBorder = true;
                                }
                                v->addVorVertex(vorv);
                            }
                        }

                        bIntersect = true;
                    }
                }
            }
            else if (CGAL::assign(l,o))
            {
                Point2D pa, pb, onpt, pts;
                Vector2D dir(l.direction().dx(), l.direction().dy());
                if (l.is_degenerate())
                {
                    ec ++;
                    continue;
                }else if (l.is_vertical())
                {		
                    onpt.x = -l.c()/l.a();
                    onpt.y = 0;
                }else if (l.is_horizontal())
                {
                    onpt.x = 0;
                    onpt.y = -l.c()/l.b();
                }else{
                    onpt.x = l.point().x();
                    onpt.y = l.point().y();	
                }

                for(int i=0; i<poly->count; i++)
                {
                    pa = Point2D(poly->acmes[i].x(),poly->acmes[i].y());
                    pb = Point2D(poly->acmes[i+1].x(),poly->acmes[i+1].y());

                    int n = LinesegLineInter(onpt, dir, pa, pb, pts);
                    if ( n < 1 ) continue;
                    else{		
                        Point pos(pts.x, pts.y);
                        vorv = FindVorVertex(pos);
                        if( !vorv )	{
                            vorv = new VorVertex;
                            vorv->type = INTERSECT;
                            vorv->pos  = pos;
                            vorv->v[0] = FindMyVertex(edge.target());
                            vorv->v[1] = FindMyVertex(edge.source());
                            vorv->bindex= i;								
                            vorvertices.push_back(vorv);
                        }
                        assert(vorv != NULL);
                        if (vorv->type == INTERSECT)
                        {
                            v->isAdjBorder = true;
                        }
                        v->addVorVertex(vorv);
                    }
                }
                bIntersect = true;
            }
            ec ++;
        } while (ec != done);

        v->SortVorVertexCCW();

        if (v->GetVorVertexNumber() < 1)
        {
            v->type = NOCELL;
            nocellNum++;
        }

        // if the voronoi diagram around the vertex is intersect with boundary, 
        // we should add more corner voronoi vertex to its adjacent array
        if( bIntersect ) {
            for( int i=0; i<v->GetVorVertexNumber(); i++) {
                VorVertex *v1, *v2, *vnext;
                v1 = v->GetVorVertex(i);
                v2 = v->GetVorVertex((i+1)%v->GetVorVertexNumber());

                if( v1->type != CIRCUMCENTER && v2->type != CIRCUMCENTER ){	
                    // need add corner voronoi vertex
                    if (v->isInsideCell)
                    {
                        if (v->adjVorVertex.GetSize() == 2 && IsTriangleCCW(v->pos, v1->pos, v2->pos)) {
                            // decide which side to insert new voronoi vertex. dmyan 2005-05-10
                            VorVertex * temp;
                            temp = v1;
                            v->adjVorVertex[i] = v2;
                            v1 = v2;
                            v->adjVorVertex[(i+1)%v->GetVorVertexNumber()] = temp;
                            v2 = temp;
                            while( !IsPointOnSameBoundary(v1->pos, v2->pos) ) {
                                vnext = FindNextVorVertexOnBoundaryCCW(v1);
                                v->adjVorVertex.InsertAt(i+1, vnext);
                                i ++;
                                v1 = vnext;
                            }
                            break;		// dmyan 2005 05 14
                        }
                        else if (! IsSameDualEdge(v1, v2) ) {
                            while( !IsPointOnSameBoundary(v1->pos, v2->pos) ) {
                                vnext = FindNextVorVertexOnBoundaryCCW(v1);
                                v->adjVorVertex.InsertAt(i+1, vnext);
                                i ++;
                                v1 = vnext;
                            }
                            //						break;		// dmyan 2005 05 14
                        }
                    }else{
                        if ((int)v->adjVorVertex.GetSize() == 2 && !IsTriangleCCW(v->pos, v1->pos, v2->pos)) {
                            // decide which side to insert new voronoi vertex. dmyan 2005-05-10
                            VorVertex * temp;
                            temp = v1;
                            v->adjVorVertex[i] = v2;
                            v1 = v2;
                            v->adjVorVertex[(i+1)%v->GetVorVertexNumber()] = temp;
                            v2 = temp;
                            while( !IsPointOnSameBoundary(v1->pos, v2->pos) ) {
                                vnext = FindNextVorVertexOnBoundaryCCW(v1);
                                v->adjVorVertex.InsertAt(i+1, vnext);
                                i ++;
                                v1 = vnext;
                            }
                            break;		// dmyan 2005 05 14
                        }
                        else if (! IsSameDualEdge(v1, v2) ) {
                            while( !IsPointOnSameBoundary(v1->pos, v2->pos) ) {
                                vnext = FindNextVorVertexOnBoundaryCCW(v1);
                                v->adjVorVertex.InsertAt(i+1, vnext);
                                i ++;
                                v1 = vnext;
                            }
                            //						break;		// dmyan 2005 05 14
                        }
                    }					
                }
            }
        } // current vertex's voronoi diagram
    }	

    if (nocellNum == (int)m_pdt.number_of_vertices())
    {
        for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
        {		
            CCVTVertex	*v = (CCVTVertex*)*vit;
            if (v->isInsideCell)
            {
                v->type = VALID;
                for(std::list<VorVertex*>::iterator vvit = vorvertices.begin(); vvit!=vorvertices.end(); vvit++)
                {
                    VorVertex * vv = (VorVertex *)*vvit;
                    if (vv->type == CORNER)
                    {
                        v->addVorVertex(vv);
                    }
                }
                v->SortVorVertexCCW();
            }
        }
    }	
}


void CPCC2D::CreateVorVertices_Rect()
{	
    int index = 0;

    // add corner vertex of bounding box to the voronoi vertex list
    
    AddVoronoiVertex(CORNER, Point(bound_x_min, bound_y_min));
    AddVoronoiVertex(CORNER, Point(bound_x_max, bound_y_min));
    AddVoronoiVertex(CORNER, Point(bound_x_max, bound_y_max));
    AddVoronoiVertex(CORNER, Point(bound_x_min, bound_y_max));	
    
    int nocellNum = 0;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = (CCVTVertex*)*vit;

        v->isAdjBorder = false;
        v->isTangent2Border = false;
        v->isActive = true;

        if (v->type == HIDDEN)
        {
            v->ClearNeighbor();
            continue;
        }

        VorVertex			*vorv;
        Edge_circulator		ec = m_pdt.incident_edges(v->vh), done = ec;
        Vertex_circulator	vc = m_pdt.incident_vertices(v->vh), donevc = vc;
        double				t,dx,dy;
        Segment				s, edge;
        Ray					r;
        Line				l;
        bool				bIntersect = false;	// whether an one ring edge of a vertex intersect with boundary

        // add one ring vertex to current vertex's adjacent vertex list.
        do {
            if (m_pdt.is_infinite(vc) ) {
                vc ++;
                continue;
            }
            
            CCVTVertex * vtx = FindMyVertex(vc);
            assert(vtx != NULL);
            v->addAdjVertex(vtx);

            vc ++;
        } while ( vc != donevc );

        // create one ring voronoi vertex for each vertex	
        do {
            if( m_pdt.is_infinite(ec) ){
                ec ++;
                continue;
            }
            
            edge = m_pdt.segment(ec);	// edge
            index= -1;

            CGAL::Object	o = m_pdt.dual(ec);

            if (CGAL::assign(s,o)) {
                if (s.is_degenerate())
                {
                    int tmp =2;
                }
                //judge whether the two weighted pts lie on the same side of the segment
                bool b1 = IsTriangleCCW(s.target(), edge.target(), s.source());
                bool b2 = IsTriangleCCW(s.target(), edge.source(), s.source());
                if ((b1&&b2) || !(b1||b2))
                {//b1 and b2 have the same signal, i.e., on the same side
                    CCVTVertex * v1 = FindMyVertex(edge.target());
                    CCVTVertex * v2 = FindMyVertex(edge.source());
                    assert(v1 != NULL);
                    assert(v2 != NULL);
                    if(v1->weight != v2->weight)
                    {					
                        if (v1->weight < v2->weight)
                        {
                            v1->isInsideCell = false;
                        }else{
                            v2->isInsideCell = false;
                        }
                    }
                }

                if(isInBox(s.source().x(),s.source().y()) && isInBox(s.target().x(),s.target().y()))
                {
                    vorv = FindVorVertex(s.source());
                    assert(vorv != NULL);
                    v->addVorVertex(vorv);
                    vorv = FindVorVertex(s.target());
                    assert(vorv != NULL);
                    v->addVorVertex(vorv);
                }
                else if(isInBox(s.target().x(),s.target().y())) {
                    dx = s.source().x()-s.target().x();
                    dy = s.source().y()-s.target().y();
                    t = intersect_rect_bound(s.target().x(),s.target().y(),dx,dy,index);

                    // add a new voronoi vertex
                    Point pos = Point(s.target().x()+t*dx,s.target().y()+t*dy);
                    vorv = FindVorVertex(pos);
                    if( !vorv )	{
                        vorv = new VorVertex;
                        vorv->type = INTERSECT;
                        vorv->pos  = pos;
                        vorv->v[0] = FindMyVertex(edge.target());
                        vorv->v[1] = FindMyVertex(edge.source());
                        vorv->bindex= index;
                        vorvertices.push_back(vorv);
                    }					
                    assert(vorv != NULL);
                    if (vorv->type == INTERSECT)
                    {
                        v->isAdjBorder = true;
                    }
                    v->addVorVertex(vorv);
                    vorv = FindVorVertex(s.target());
                    assert(vorv != NULL);
                    v->addVorVertex(vorv);
                    bIntersect = true;
                }
                else if(isInBox(s.source().x(),s.source().y()))
                {
                    dx = s.target().x()-s.source().x();
                    dy = s.target().y()-s.source().y();
                    t = intersect_rect_bound(s.source().x(),s.source().y(),dx,dy, index);
                    // add a new voronoi vertex
                    Point pos = Point(s.source().x()+t*dx,s.source().y()+t*dy);
                    vorv = FindVorVertex(pos);
                    if( vorv != NULL ){
                        assert(vorv != NULL);
                        double x = vorv->pos.x();
                    }
                    else if( vorv == NULL )
                    {
                        vorv = new VorVertex;
                        vorv->type = INTERSECT;
                        vorv->pos  = Point(s.source().x()+t*dx,s.source().y()+t*dy);
                        vorv->v[0] = FindMyVertex(edge.target());
                        vorv->v[1] = FindMyVertex(edge.source());
                        vorv->bindex= index;
                        vorvertices.push_back(vorv);
                        assert(vorv != NULL);
                        double x = vorv->pos.x();
                    }
                    assert(vorv != NULL);
                    if (vorv->type == INTERSECT)
                    {
                        v->isAdjBorder = true;
                    }
                    v->addVorVertex(vorv);
                    vorv = FindVorVertex(s.source());
                    assert(vorv != NULL);
                    v->addVorVertex(vorv);
                    bIntersect = true;
                }

                else {	// the start point the ray lie outside the bounding box. so has two intersection voronoi points
                    Point2D linesp(s.source().x(),s.source().y()), pa, pb, pts;
                    Vector2D dir(s.target().x() - s.source().x(), s.target().y()-s.source().y());

                    for(int i=0; i<4; i++)	{
                        switch (i) {
                            case 0:
                                pa = Point2D(bound_x_min, bound_y_min);
                                pb = Point2D(bound_x_max, bound_y_min);
                                break;
                            case 1:
                                pa = Point2D(bound_x_max, bound_y_min);
                                pb = Point2D(bound_x_max, bound_y_max);
                                break;
                            case 2:
                                pa = Point2D(bound_x_max, bound_y_max);
                                pb = Point2D(bound_x_min, bound_y_max);
                                break;
                            case 3:
                                pa = Point2D(bound_x_min, bound_y_max);
                                pb = Point2D(bound_x_min, bound_y_min);
                                break;
                        }

                        int n = LinesegLineInter(linesp, dir, pa, pb, pts);
                        if ( n < 1 ) continue;
                        else	{
                            Point pos(pts.x, pts.y);
                            double t1 = (pts.x - s.source().x())/(s.target().x() - s.source().x());
                            double t2 = (pts.y - s.source().y())/(s.target().y() - s.source().y());
                            if( t1 < 0 || t1 > 1 || t2 < 0 || t2 > 1) continue;
                            vorv = FindVorVertex(pos);
                            if( !vorv )	{
                                vorv = new VorVertex;
                                vorv->type = INTERSECT;
                                vorv->pos  = pos;
                                vorv->v[0] = FindMyVertex(edge.target());
                                vorv->v[1] = FindMyVertex(edge.source());
                                vorv->bindex= i;								
                                vorvertices.push_back(vorv);
                            }
                            assert(vorv != NULL);
                            if (vorv->type == INTERSECT)
                            {
                                v->isAdjBorder = true;
                            }
                            v->addVorVertex(vorv);
                            bIntersect = true;
                        }
                    }
                }
            }	

            else if (CGAL::assign(r,o)) {
                if (r.is_degenerate())
                {
                    ec++;
                    continue;
                }else{
                    //judge whether the two weighted pts lie on the same side of the segment
                    Point pt_on_r = Point(r.source().x()+r.direction().dx(),r.source().y()+r.direction().dy());
                    bool b1 = IsTriangleCCW(pt_on_r, edge.target(), r.source());
                    bool b2 = IsTriangleCCW(pt_on_r, edge.source(), r.source());
                    if ((b1&&b2) || !(b1||b2))
                    {//b1 and b2 have the same signal, i.e., on the same side
                        CCVTVertex * v1 = FindMyVertex(edge.target());
                        CCVTVertex * v2 = FindMyVertex(edge.source());
                        assert(v1 != NULL);
                        assert(v2 != NULL);
                        if(v1->weight != v2->weight)
                        {					
                            if (v1->weight < v2->weight)
                            {
                                v1->isInsideCell = false;
                            }else{
                                v2->isInsideCell = false;
                            }
                        }
                    }

                    if(isInBox(r.source().x(),r.source().y()))
                    {
                        t = intersect_rect_bound(r.source().x(),r.source().y(),r.direction().dx(),r.direction().dy(),index);					
                        // add a new voronoi vertex
                        Point pos = Point(r.source().x()+t*r.direction().dx(),r.source().y()+t*r.direction().dy());
                        vorv = FindVorVertex(pos);
                        if( !vorv )	{
                            vorv = new VorVertex;
                            vorv->type = INTERSECT;
                            vorv->pos  = pos;
                            vorv->v[0] = FindMyVertex(edge.target());
                            vorv->v[1] = FindMyVertex(edge.source());
                            vorv->bindex= index;
                            vorvertices.push_back(vorv);
                        }
                        assert(vorv != NULL);
                        if (vorv->type == INTERSECT)
                        {
                            v->isAdjBorder = true;
                        }
                        v->addVorVertex(vorv);
                        vorv = FindVorVertex(r.source());
                        assert(vorv != NULL);
                        v->addVorVertex(vorv);
                        bIntersect = true;
                    }
                    else {	// the start point the ray lie outside the bounding box. so has two intersection voronoi points
                        Point2D linesp(r.source().x(),r.source().y()), pa, pb, pts;
                        Vector2D dir(r.direction().dx(),r.direction().dy());

                        for(int i=0; i<4; i++)	{
                            switch (i) {
                                case 0:
                                    pa = Point2D(bound_x_min, bound_y_min);
                                    pb = Point2D(bound_x_max, bound_y_min);
                                    break;
                                case 1:
                                    pa = Point2D(bound_x_max, bound_y_min);
                                    pb = Point2D(bound_x_max, bound_y_max);
                                    break;
                                case 2:
                                    pa = Point2D(bound_x_max, bound_y_max);
                                    pb = Point2D(bound_x_min, bound_y_max);
                                    break;
                                case 3:
                                    pa = Point2D(bound_x_min, bound_y_max);
                                    pb = Point2D(bound_x_min, bound_y_min);
                                    break;
                            }

                            int n = LinesegLineInter(linesp, dir, pa, pb, pts);
                            if ( n < 1 ) continue;
                            else{
                                if(dir * Vector2D(pts - linesp) <= 0)	continue;

                                Point pos(pts.x, pts.y);
                                vorv = FindVorVertex(pos);
                                if( !vorv )	{
                                    vorv = new VorVertex;
                                    vorv->type = INTERSECT;
                                    vorv->pos  = pos;
                                    vorv->v[0] = FindMyVertex(edge.target());
                                    vorv->v[1] = FindMyVertex(edge.source());
                                    vorv->bindex= i;								
                                    vorvertices.push_back(vorv);
                                }
                                assert(vorv != NULL);
                                if (vorv->type == INTERSECT)
                                {
                                    v->isAdjBorder = true;
                                }
                                v->addVorVertex(vorv);
                            }
                        }
                        bIntersect = true;
                    }
                }
            }
            else if (CGAL::assign(l,o))
            {
                Point2D pa, pb, onpt, pts;
                Vector2D dir(l.direction().dx(), l.direction().dy());
                if (l.is_degenerate())
                {
                    ec ++;
                    continue;
                }else if (l.is_vertical())
                {		
                    onpt.x = -l.c()/l.a();
                    onpt.y = 0;
                }else if (l.is_horizontal())
                {
                    onpt.x = 0;
                    onpt.y = -l.c()/l.b();
                }else{
                    onpt.x = l.point().x();
                    onpt.y = l.point().y();	
                }

                for(int i=0; i<4; i++)	{
                    switch (i) {
                        case 0:
                            pa = Point2D(bound_x_min, bound_y_min);
                            pb = Point2D(bound_x_max, bound_y_min);
                            break;
                        case 1:
                            pa = Point2D(bound_x_max, bound_y_min);
                            pb = Point2D(bound_x_max, bound_y_max);
                            break;
                        case 2:
                            pa = Point2D(bound_x_max, bound_y_max);
                            pb = Point2D(bound_x_min, bound_y_max);
                            break;
                        case 3:
                            pa = Point2D(bound_x_min, bound_y_max);
                            pb = Point2D(bound_x_min, bound_y_min);
                            break;
                    }
                    int n = LinesegLineInter(onpt, dir, pa, pb, pts);
                    if ( n < 1 ) continue;
                    else{		
                        Point pos(pts.x, pts.y);
                        vorv = FindVorVertex(pos);
                        if( !vorv )	{
                            vorv = new VorVertex;
                            vorv->type = INTERSECT;
                            vorv->pos  = pos;
                            vorv->v[0] = FindMyVertex(edge.target());
                            vorv->v[1] = FindMyVertex(edge.source());
                            vorv->bindex= i;								
                            vorvertices.push_back(vorv);
                        }
                        assert(vorv != NULL);
                        if (vorv->type == INTERSECT)
                        {
                            v->isAdjBorder = true;
                        }
                        v->addVorVertex(vorv);
                    }
                }
                bIntersect = true;
            }
            ec ++;
        } while (ec != done);

        v->SortVorVertexCCW();
                 
        if (v->GetVorVertexNumber() < 1)
        {
            v->type = NOCELL;
            nocellNum++;
        }

        // if the voronoi diagram around the vertex is intersect with boundary, 
        // we should add more corner voronoi vertex to its adjacent array
        if( bIntersect ) {
            for( int i=0; i<v->GetVorVertexNumber(); i++) {
                VorVertex *v1, *v2, *vnext;
                v1 = v->GetVorVertex(i);
                v2 = v->GetVorVertex((i+1)%v->GetVorVertexNumber());

                if( v1->type != CIRCUMCENTER && v2->type != CIRCUMCENTER ){	
                    // need add corner voronoi vertex
                    if (v->isInsideCell)
                    {
                        if (v->adjVorVertex.GetSize() == 2 && IsTriangleCCW(v->pos, v1->pos, v2->pos)) {
                            // decide which side to insert new voronoi vertex. dmyan 2005-05-10
                            VorVertex * temp;
                            temp = v1;
                            v->adjVorVertex[i] = v2;
                            v1 = v2;
                            v->adjVorVertex[(i+1)%v->GetVorVertexNumber()] = temp;
                            v2 = temp;
                            while( !IsPointOnSameBoundary(v1->pos, v2->pos) ) {
                                vnext = FindNextVorVertexOnBoundaryCCW(v1);
                                v->adjVorVertex.InsertAt(i+1, vnext);
                                i ++;
                                v1 = vnext;
                            }
                            break;		// dmyan 2005 05 14
                        }
                        else if (! IsSameDualEdge(v1, v2) ) {
                            while( !IsPointOnSameBoundary(v1->pos, v2->pos) ) {
                                vnext = FindNextVorVertexOnBoundaryCCW(v1);
                                v->adjVorVertex.InsertAt(i+1, vnext);
                                i ++;
                                v1 = vnext;
                            }
    //						break;		// dmyan 2005 05 14
                        }
                    }else{
                        if ((int)v->adjVorVertex.GetSize() == 2 && !IsTriangleCCW(v->pos, v1->pos, v2->pos)) {
                            // decide which side to insert new voronoi vertex. dmyan 2005-05-10
                            VorVertex * temp;
                            temp = v1;
                            v->adjVorVertex[i] = v2;
                            v1 = v2;
                            v->adjVorVertex[(i+1)%v->GetVorVertexNumber()] = temp;
                            v2 = temp;
                            while( !IsPointOnSameBoundary(v1->pos, v2->pos) ) {
                                vnext = FindNextVorVertexOnBoundaryCCW(v1);
                                v->adjVorVertex.InsertAt(i+1, vnext);
                                i ++;
                                v1 = vnext;
                            }
                            break;		// dmyan 2005 05 14
                        }
                        else if (! IsSameDualEdge(v1, v2) ) {
                            while( !IsPointOnSameBoundary(v1->pos, v2->pos) ) {
                                vnext = FindNextVorVertexOnBoundaryCCW(v1);
                                v->adjVorVertex.InsertAt(i+1, vnext);
                                i ++;
                                v1 = vnext;
                            }
                            //						break;		// dmyan 2005 05 14
                        }
                    }					
                }
            }
        } // current vertex's voronoi diagram
    }	

    if (nocellNum == (int)m_pdt.number_of_vertices())
    {
        for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
        {		
            CCVTVertex	*v = (CCVTVertex*)*vit;
            if (v->isInsideCell)
            {
                v->type = VALID;
                for(std::list<VorVertex*>::iterator vvit = vorvertices.begin(); vvit!=vorvertices.end(); vvit++)
                {
                    VorVertex * vv = (VorVertex *)*vvit;
                    if (vv->type == CORNER)
                    {
                        v->addVorVertex(vv);
                    }
                }
                v->SortVorVertexCCW();
            }
        }
    }	
}

void CPCC2D::ComputeCentroids()
{
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = (CCVTVertex*)*vit;
        
        if (v->type != VALID){
            continue;
        }
        
        const int adjVorVNum = (int)v->adjVorVertex.GetSize();
        if (adjVorVNum <= 1)
        {
            continue;
        }

        Kernel::FT	w(0.0);
        double x = 0, y =0;
        double len = 0;
        Point2D p = Point2D( v->pos.x(),  v->pos.y());

        if ((m_boundType == RECTBOUND) ||(m_boundType == POLYBOUND)||((m_boundType == CIRCLEBOUND) && (!v->isAdjBorder)))
        {
            for(int i=0; i<adjVorVNum; i++)
            {
                VorVertex * vv1, * vv2;
                vv1 = (VorVertex *)v->adjVorVertex.GetAt(i);
                vv2 = (VorVertex *)v->adjVorVertex.GetAt( (i+1)%adjVorVNum);

                Kernel::FT a = Triangle(vv1->pos, vv2->pos, v->vh->point()).area();

                //		//******test for weighted CVT******
                //		/*double wx = 0, wy = 0;
                //		wcvtCenter(vv1->pos, vv2->pos, v->pos, wx, wy);
                //		Point c = Point(a*wx, a*wy);*/
                //		//*********************************

                Point c = Point(a*(vv1->pos.x() + vv2->pos.x() + v->pos.x())/3, a*(vv1->pos.y() + vv2->pos.y() + v->pos.y())/3);
                x += c.x(); y += c.y();
                w += a;

                Point2D p1 = Point2D( vv1->pos.x(),  vv1->pos.y());
                len += Distance(p1,p);
            }
            v->area = w;
        }
        else if (m_boundType == CIRCLEBOUND)
        {
            assert (v->isAdjBorder);
            
            if (adjVorVNum < 3)
            {
                assert(adjVorVNum == 2);
                
                VorVertex * vv1, * vv2;
                vv1 = (VorVertex *)v->adjVorVertex.GetAt(0);
                vv2 = (VorVertex *)v->adjVorVertex.GetAt(1);

                double cx, cy;
                double dist = disPt2Line(vv1->pos.x(), vv1->pos.y(), vv2->pos.x(), vv2->pos.y(), center_x, center_y, cx, cy);
                double tmpd = (cx-center_x)*(cx-center_x) + (cy-center_y)*(cy-center_y);

                Point2D p1 = Point2D( vv1->pos.x(),  vv1->pos.y());
                Point2D p2 = Point2D( vv2->pos.x(),  vv2->pos.y());				
                len = Distance(p1,p) + Distance(p2,p);

                if (abs(tmpd) > TOLER)
                {
                    Point cent;
                    CircularSegmentCentroidCCW(vv1->pos, vv2->pos, w, cent);
                    x = cent.x()*w; y = cent.y()*w;
                    v->area = w;
                }else{
                    double nord = 0.5*radius_r/ sqrt((v->pos.x()-center_x)*(v->pos.x()-center_x)+(v->pos.y()-center_y)*(v->pos.y()-center_y));
                    x = center_x + (v->pos.x()-center_x)*nord;
                    y = center_y + (v->pos.y()-center_y)*nord;
                    w = 1;
                    v->area = 0.5*radius_r*radius_r*PI;
                }
            }else
            {
                double newx = 0, newy = 0;
                for (int i =0; i<adjVorVNum; i++)
                {
                    VorVertex * vv;
                    vv = (VorVertex *)v->adjVorVertex.GetAt(i);
                    newx += vv->pos.x();
                    newy += vv->pos.y();

                    Point2D p1 = Point2D( vv->pos.x(),  vv->pos.y());					
                    len += Distance(p1,p);
                }
                newx /= adjVorVNum;
                newy /= adjVorVNum;
                Point newc = Point(newx, newy);

                for(int i=0; i<adjVorVNum; i++)
                {
                    VorVertex * vv1, * vv2;
                    vv1 = (VorVertex *)v->adjVorVertex.GetAt(i);
                    vv2 = (VorVertex *)v->adjVorVertex.GetAt((i+1)%adjVorVNum);

                    Kernel::FT a = abs(Triangle(vv1->pos, vv2->pos, newc).area());

                    Point c = Point(a*(vv1->pos.x() + vv2->pos.x() + newx)/3, a*(vv1->pos.y() + vv2->pos.y() + newy)/3);
                    x += c.x(); y += c.y();
                    w += a;

                    if (vv1->type == INTERSECT && vv2->type == INTERSECT)
                    {
                        if ((vv1->v[0] == vv2->v[0]) &&(vv1->v[1] == vv2->v[1])){
                            continue;
                        }
                        
                        Point cent;
                        double area;
                        CircularSegmentCentroidCCW(vv1->pos, vv2->pos, area, cent);
                        x += cent.x()*area; y += cent.y()*area;
                        w += area;
                    }
                }
                v->area = w;
            }
        }
        x /= w; y /= w;
        v->centroid = Point(x,y);
        v->avglen = len/adjVorVNum;
    }
}

//ignore the sequence of pt1 and pt2, always keep the smaller circular segment
void CPCC2D::CircularSegmentCentroid(Point pt1, Point pt2, double& area, Point& centroid)
{
    double cx, cy;
    double cx_mean = (pt1.x() + pt2.x())*0.5;
    double cy_mean = (pt1.y() + pt2.y())*0.5;

    double dist = disPt2Line(pt1.x(), pt1.y(), pt2.x(), pt2.y(), center_x, center_y, cx, cy);
    double h = radius_r - dist;
    if (h==0)
    {
        area = 0;
        centroid = Point(cx_mean, cy_mean);
        return;
    }

    double c = 2*sqrt(h*(2*radius_r-h));
    double theta = 2*asin(c/radius_r*0.5);
    assert(theta >0);

    Kernel::FT a_sector = 0.5*radius_r*radius_r*theta;
    Kernel::FT a_tri = abs(Triangle(pt1, pt2, Point(center_x, center_y)).area());
    
    area = a_sector - a_tri;	

    double cy_sector_tmp = 4*radius_r*sin(theta*0.5)/(3*theta);
    assert(cy_sector_tmp>0);

    double tmpd = sqrt((cx_mean-center_x)*(cx_mean-center_x)+(cy_mean-center_y)*(cy_mean-center_y));

    double nord = cy_sector_tmp/tmpd;
    double cx_sector = center_x + (cx_mean-center_x)*nord;
    double cy_sector = center_y + (cy_mean-center_y)*nord;

    double cx_tri = (pt1.x() + pt2.x() + center_x)/3;
    double cy_tri = (pt1.y() + pt2.y() + center_y)/3;
    
    //centroid  = Point((a_sector*cx_sector-a_tri*cx_tri), (a_sector*cy_sector-a_tri*cy_tri));
    //the real centroid should be:
    centroid  = Point((a_sector*cx_sector-a_tri*cx_tri)/area, (a_sector*cy_sector-a_tri*cy_tri)/area);
}

//suppose that pt1 and pt2 are on the circle in counterclockwise sequence
void CPCC2D::CircularSegmentCentroidCCW(Point pt1, Point pt2, double& area, Point& centroid)
{
    double cx, cy;
    double cx_mean = (pt1.x() + pt2.x())*0.5;
    double cy_mean = (pt1.y() + pt2.y())*0.5;

    double dist = disPt2Line(pt1.x(), pt1.y(), pt2.x(), pt2.y(), center_x, center_y, cx, cy);
    double h = radius_r - dist;
    if (h==0)
    {
        area = 0;
        centroid = Point(cx_mean, cy_mean);
        return;
    }

    double c = 2*sqrt(h*(2*radius_r-h));
    double theta = 2*asin(c/radius_r*0.5);
    assert(theta >0);

    Kernel::FT a_sector = 0.5*radius_r*radius_r*theta;
    Kernel::FT a_tri = abs(Triangle(pt1, pt2, Point(center_x, center_y)).area());

    area = a_sector - a_tri;	

    double cy_sector_tmp = 4*radius_r*sin(theta*0.5)/(3*theta);
    assert(cy_sector_tmp>0);

    double tmpd = sqrt((cx_mean-center_x)*(cx_mean-center_x)+(cy_mean-center_y)*(cy_mean-center_y));

    double nord = cy_sector_tmp/tmpd;
    double cx_sector = center_x + (cx_mean-center_x)*nord;
    double cy_sector = center_y + (cy_mean-center_y)*nord;

    double cx_tri = (pt1.x() + pt2.x() + center_x)/3;
    double cy_tri = (pt1.y() + pt2.y() + center_y)/3;

    //centroid  = Point((a_sector*cx_sector-a_tri*cx_tri), (a_sector*cy_sector-a_tri*cy_tri));
    //the real centroid should be:
    centroid  = Point((a_sector*cx_sector-a_tri*cx_tri)/area, (a_sector*cy_sector-a_tri*cy_tri)/area);

    //determine whether we should take the bigger circular segment or the smaller one
    float x1 = (float)(pt1.x() - center_x);
    float y1 = (float)(pt1.y() - center_y);
    float x2 = (float)(pt2.x() - center_x);
    float y2 = (float)(pt2.y() - center_y);

    float da1 = atan2 (y1, x1);
    float da2 = atan2 (y2, x2);

    float d_angle = da2 - da1; 
    if (d_angle < 0)
    {
        d_angle = static_cast<float>(PI2 + d_angle);
    }

    if(d_angle > PI){
        double ca = PI*radius_r*radius_r;
        double x = center_x*ca - centroid.x()*area; 
        double y = center_y*ca - centroid.y()*area;
        area = ca - area;
        centroid  = Point(x/area, y/area);
    }
}

void CPCC2D::CreateVorVertices_Circle()
{	
    int index = 0;

    int nocellNum = 0;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {		
        CCVTVertex	*v = (CCVTVertex*)*vit;

        v->isAdjBorder = false;
        v->isTangent2Border = false;
        v->isActive = true;

        if (v->type == HIDDEN)
        {
            v->ClearNeighbor();
            continue;
        }

        VorVertex			*vorv;
        Edge_circulator		ec = m_pdt.incident_edges(v->vh), done = ec;
        Vertex_circulator	vc = m_pdt.incident_vertices(v->vh), donevc = vc;
        double				t,dx,dy;
        Segment				s, edge;
        Ray					r;
        Line				l;
        bool				bIntersect = false;	// whether an one ring edge of a vertex intersect with boundary

        // add one ring vertex to current vertex's adjacent vertex list.
        do {
            if (m_pdt.is_infinite(vc) ) {
                vc ++;
                continue;
            }

            CCVTVertex * vtx = FindMyVertex(vc);
            assert(vtx != NULL);
            v->addAdjVertex(vtx);

            vc ++;
        } while ( vc != donevc );

        // create one ring voronoi vertex for each vertex	
        do {
            if( m_pdt.is_infinite(ec) ){
                ec ++;
                continue;
            }

            edge = m_pdt.segment(ec);	// edge
            index= -1;

            CGAL::Object o = m_pdt.dual(ec);

            if (CGAL::assign(s,o)) {
                if (s.is_degenerate())
                {
                    int tmp =2;
                }
                //judge whether the two weighted pts lie on the same side of the segment
                bool b1 = IsTriangleCCW(s.target(), edge.target(), s.source());
                bool b2 = IsTriangleCCW(s.target(), edge.source(), s.source());
                if ((b1&&b2) || !(b1||b2))
                {//b1 and b2 have the same signal, i.e., on the same side
                    CCVTVertex * v1 = FindMyVertex(edge.target());
                    CCVTVertex * v2 = FindMyVertex(edge.source());
                    assert(v1 != NULL);
                    assert(v2 != NULL);
                    if(v1->weight != v2->weight)
                    {					
                        if (v1->weight < v2->weight)
                        {
                            v1->isInsideCell = false;
                        }else{
                            v2->isInsideCell = false;
                        }
                    }
                }

                if(isInCircle(s.source().x(),s.source().y()) && isInCircle(s.target().x(),s.target().y()))
                {
                    vorv = FindVorVertex(s.source());
                    assert(vorv != NULL);
                    v->addVorVertex(vorv);
                    vorv = FindVorVertex(s.target());
                    assert(vorv != NULL);
                    v->addVorVertex(vorv);
                }
                else if(isInCircle(s.target().x(),s.target().y())) {
                    dx =s.source().x()-s.target().x();
                    dy = s.source().y()-s.target().y();
                    
                    bool binter = intersect_circle_bound(s.target().x(),s.target().y(),dx,dy, t);
                    assert(binter);
                    assert(t>0);

                    // add a new voronoi vertex
                    Point	pos = Point(s.target().x()+t*dx,s.target().y()+t*dy);
                    vorv = FindVorVertex(pos);
                    if( !vorv )	{
                        vorv = new VorVertex;
                        vorv->type = INTERSECT;
                        vorv->pos  = pos;
                        vorv->v[0] = FindMyVertex(edge.target());
                        vorv->v[1] = FindMyVertex(edge.source());
                        vorv->bindex= -1;
                        vorvertices.push_back(vorv);
                    }					
                    assert(vorv != NULL);
                    if (vorv->type == INTERSECT)
                    {
                        v->isAdjBorder = true;
                    }
                    v->addVorVertex(vorv);
                    vorv = FindVorVertex(s.target());
                    assert(vorv != NULL);
                    v->addVorVertex(vorv);
                    bIntersect = true;
                }
                else if(isInCircle(s.source().x(),s.source().y()))
                {
                    dx = s.target().x()-s.source().x();
                    dy = s.target().y()-s.source().y();
                    
                    bool binter = intersect_circle_bound(s.source().x(),s.source().y(),dx,dy,t);
                    assert(binter);
                    assert(t>0);

                    // add a new voronoi vertex
                    Point pos = Point(s.source().x()+t*dx,s.source().y()+t*dy);
                    vorv = FindVorVertex(pos);
                    if( vorv != NULL ){
                        assert(vorv != NULL);
                        double x = vorv->pos.x();
                    }
                    else if( vorv == NULL )
                    {
                        vorv = new VorVertex;
                        vorv->type = INTERSECT;
                        vorv->pos  = Point(s.source().x()+t*dx,s.source().y()+t*dy);
                        vorv->v[0] = FindMyVertex(edge.target());
                        vorv->v[1] = FindMyVertex(edge.source());
                        vorv->bindex= -1;
                        vorvertices.push_back(vorv);
                    }
                    assert(vorv != NULL);
                    if (vorv->type == INTERSECT)
                    {
                        v->isAdjBorder = true;
                    }
                    v->addVorVertex(vorv);
                    vorv = FindVorVertex(s.source());
                    assert(vorv != NULL);
                    v->addVorVertex(vorv);
                    bIntersect = true;
                }
                else {	// the start point the ray lie outside the bounding circle. so has two intersection voronoi points
                    //OR the segment has no intersection with the boundary!
                    dx = s.target().x()-s.source().x();
                    dy = s.target().y()-s.source().y();
                    double t1, t2;
                    Point pt[2];
                    bool binter = intersect_circle_bound(s.source().x(),s.source().y(),dx,dy, t1,t2);
                    if(binter && (t2>0) &&(t1<1))
                    {
                        pt[0] = Point(s.source().x()+t1*dx, s.source().y()+t1*dy);
                        pt[1] = Point(s.source().x()+t2*dx, s.source().y()+t2*dy);

                        for(int i=0; i<2; i++)	{
                            vorv = FindVorVertex(pt[i]);
                            if( !vorv )	{
                                vorv = new VorVertex;
                                vorv->type = INTERSECT;
                                vorv->pos  = pt[i];
                                vorv->v[0] = FindMyVertex(edge.target());
                                vorv->v[1] = FindMyVertex(edge.source());
                                vorv->bindex= -1;
                                vorvertices.push_back(vorv);
                            }							
                            assert(vorv != NULL);
                            if (vorv->type == INTERSECT)
                            {
                                v->isAdjBorder = true;
                            }
                            v->addVorVertex(vorv);
                        }
                    }
                }
            }	

            else if (CGAL::assign(r,o)) {
                if (r.is_degenerate())
                {
                    ec ++;
                    continue;
                //}else if (r.is_vertical())
                //{//how to deal with?	
                //	double dx = r.direction().dx();
                //	double dy = r.direction().dy();
                //	int tmp =2;
                //}else if (r.is_horizontal())
                //{//how to deal with?
                //	double dx = r.direction().dx();
                //	double dy = r.direction().dy();
                //	double tmpx1 = r.point(0).x();
                //	double tmpy1 = r.point(0).y();
                //	double tmpx2 = r.point(1).x();
                //	double tmpy2 = r.point(1).y();
                //	int tmp =2;
                }else{
                    //judge whether the two weighted pts lie on the same side of the segment
                    Point pt_on_r = Point(r.source().x()+r.direction().dx(),r.source().y()+r.direction().dy());
                    bool b1 = IsTriangleCCW(pt_on_r, edge.target(), r.source());
                    bool b2 = IsTriangleCCW(pt_on_r, edge.source(), r.source());
                    if ((b1&&b2) || !(b1||b2))
                    {//b1 and b2 have the same signal, i.e., on the same side
                        CCVTVertex * v1 = FindMyVertex(edge.target());
                        CCVTVertex * v2 = FindMyVertex(edge.source());
                        assert(v1 != NULL);
                        assert(v2 != NULL);
                        if(v1->weight != v2->weight)
                        {					
                            if (v1->weight < v2->weight)
                            {
                                v1->isInsideCell = false;
                            }else{
                                v2->isInsideCell = false;
                            }
                        }
                    }

                    if(isInCircle(r.source().x(),r.source().y()))
                    {
                        bool binter = intersect_circle_bound(r.source().x(),r.source().y(),r.direction().dx(),r.direction().dy(), t);					
                        assert(binter);
                        assert(t>0);

                        // add a new voronoi vertex
                        Point pos = Point(r.source().x()+t*r.direction().dx(),r.source().y()+t*r.direction().dy());
                        vorv = FindVorVertex(pos);
                        if( !vorv )	{
                            vorv = new VorVertex;
                            vorv->type = INTERSECT;
                            vorv->pos  = pos;
                            vorv->v[0] = FindMyVertex(edge.target());
                            vorv->v[1] = FindMyVertex(edge.source());
                            vorv->bindex= -1;
                            vorvertices.push_back(vorv);
                        }
                        assert(vorv != NULL);
                        if (vorv->type == INTERSECT)
                        {
                            v->isAdjBorder = true;
                        }
                        v->addVorVertex(vorv);
                        vorv = FindVorVertex(r.source());
                        assert(vorv != NULL);
                        v->addVorVertex(vorv);
                        bIntersect = true;
                    }
                    else {
                        // the start point the ray lie outside the bounding box. so has two intersection voronoi points					
                        //OR the ray has no intersection with the boundary!
                        double dx = r.direction().dx();
                        double dy = r.direction().dy();
                        double t1, t2;
                        bool binter = intersect_circle_bound(r.source().x(),r.source().y(), dx, dy, t1, t2);
                        Point pt[2];
                        if (t2>0 && binter)
                        {
                            pt[0] = Point(r.source().x()+t1*dx, r.source().y()+t1*dy);
                            pt[1] = Point(r.source().x()+t2*dx, r.source().y()+t2*dy);

                            for(int i=0; i<2; i++)	{
                                vorv = FindVorVertex(pt[i]);
                                if( !vorv )	{
                                    vorv = new VorVertex;
                                    vorv->type = INTERSECT;
                                    vorv->pos  = pt[i];
                                    vorv->v[0] = FindMyVertex(edge.target());
                                    vorv->v[1] = FindMyVertex(edge.source());
                                    vorv->bindex= -1;
                                    vorvertices.push_back(vorv);
                                }							
                                assert(vorv != NULL);
                                if (vorv->type == INTERSECT)
                                {
                                    v->isAdjBorder = true;
                                }
                                v->addVorVertex(vorv);
                            }
                        }
                        bIntersect = true;
                    }
                }
                
            }
            else if (CGAL::assign(l,o))
            {
                double onptx, onpty;
                if (l.is_degenerate())
                {
                    ec ++;
                    continue;
                }else if (l.is_vertical())
                {
                    onptx = -l.c()/l.a();
                    onpty = 0;
                }else if (l.is_horizontal())
                {
                    onptx = 0;
                    onpty = -l.c()/l.b();
                }else{
                    onptx = l.point().x();
                    onpty = l.point().y();	
                }
                double dx = l.direction().dx();
                double dy = l.direction().dy();
                double t1, t2;
                bool binter = intersect_circle_bound(onptx, onpty, dx, dy, t1, t2);
                Point pt[2];
                if (binter)
                {
                    pt[0] = Point(onptx+t1*dx, onpty+t1*dy);
                    pt[1] = Point(onptx+t2*dx, onpty+t2*dy);

                    for(int i=0; i<2; i++)	{
                        vorv = FindVorVertex(pt[i]);
                        if( !vorv )	{
                            vorv = new VorVertex;
                            vorv->type = INTERSECT;
                            vorv->pos  = pt[i];
                            vorv->v[0] = FindMyVertex(edge.target());
                            vorv->v[1] = FindMyVertex(edge.source());
                            vorv->bindex= -1;
                            vorvertices.push_back(vorv);
                        }
                        assert(vorv != NULL);
                        if (vorv->type == INTERSECT)
                        {
                            v->isAdjBorder = true;
                        }
                        v->addVorVertex(vorv);
                    }
                    bIntersect = true;
                }
            }
            ec ++;
        } while (ec != done);			

        if (v->GetVorVertexNumber() ==2)
        {
            VorVertex * vv1, * vv2;
            vv1 = (VorVertex *)v->adjVorVertex.GetAt(0);
            vv2 = (VorVertex *)v->adjVorVertex.GetAt(1);
            bool isCCW = IsTriangleCCW(vv1->pos, v->pos, vv2->pos);
            if (!isCCW)
            {
                v->ClearNeighbor();	
                v->addVorVertex(vv2);
                v->addVorVertex(vv1);
            }
        }else{
            v->SortVorVertexCCW();
        }

        if (v->GetVorVertexNumber() < 1)
        {
            v->type = NOCELL;
            nocellNum++;
        }
    }

    if (nocellNum == (int)m_pdt.number_of_vertices())
    {
        for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
        {		
            CCVTVertex	*v = (CCVTVertex*)*vit;
            if (v->isInsideCell)
            {
                v->type = VALID;
                for(std::list<VorVertex*>::iterator vvit = vorvertices.begin(); vvit!=vorvertices.end(); vvit++)
                {
                    VorVertex * vv = (VorVertex *)*vvit;
                    if (vv->type == CORNER)
                    {
                        v->addVorVertex(vv);
                    }
                }
                v->SortVorVertexCCW();
            }
        }
    }
    //OncircleAssertion(); //only for test, comment later
}

void CPCC2D::OncircleAssertion() //only for test
{
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = (CCVTVertex*)*vit;
        int num =0;
        for(int i=0; i<v->adjVorVertex.GetSize(); i++)
        {		
            VorVertex * vv1;
            vv1 = (VorVertex *)v->adjVorVertex.GetAt(i);
            if (vv1->type == INTERSECT)
            {
                num++;
            }
        }
        //it is possible to be here, as cut the circle into three pieces by two parallel lines
        assert(num <=2);
    }
}

void CPCC2D::InitTwoVertex_Circle()
{
    assert((int)m_pdt.number_of_faces()==0);

    CCVTVertex * vertex[2];
    VorVertex * vorv;
    VorVertex * vorv1[2];

    int index = 0;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit!=vertices.end(); ++vit)
    {
        vertex[index++] = *vit;
    }

    Line l;
    double dx, dy;

    Edge_iterator eit = m_pdt.edges_begin();
    for ( ; eit !=m_pdt.edges_end(); ++eit)
    {//should be only one edge
        CGAL::Object o = m_pdt.dual(eit);
        if(CGAL::assign(l,o)){
            dx = l.direction().dx();
            dy = l.direction().dy();
        }
    }

    // calculate bisect perpendicular
    Point2D	 onpt;

    if (l.is_degenerate())
    {
        vertex[0]->type = NOCELL;
        vertex[1]->type = NOCELL;
        return;
    }else if (l.is_vertical())
    {		
        onpt.x = -l.c()/l.a();
        onpt.y = 0;
    }else if (l.is_horizontal())
    {
        onpt.x = 0;
        onpt.y = -l.c()/l.b();
    }else{
        onpt.x = l.point().x();
        onpt.y = l.point().y();	
    }

    double	t1, t2;
    bool bintersect = intersect_circle_bound(onpt.x, onpt.y, dx, dy,t1,t2);
    if (bintersect){
        Point	pos[2];
        pos[0] = Point(onpt.x+t1*dx, onpt.y+t1*dy);
        pos[1] = Point(onpt.x+t2*dx, onpt.y+t2*dy);
        for(int i=0; i<2; i++)
        {
            vorv = FindVorVertex(pos[i]);
            if( !vorv )	{
                vorv = new VorVertex;
                vorv->type = INTERSECT;
                vorv->pos  = pos[i];
                vorv->v[0] = vertex[0];
                vorv->v[1] = vertex[1];
                vorv->bindex = -1;
                vorvertices.push_back(vorv);
            }
            assert(vorv != NULL);
            vertex[0]->addVorVertex(vorv);
            vertex[1]->addVorVertex(vorv);
            vorv1[i] = vorv;
        }

        bool isCCW[2];
        for (int i =0; i<2; i++)
        {
            isCCW[i] = IsTriangleCCW(vorv1[0]->pos, vertex[i]->pos, vorv1[1]->pos);
        }
        
        if ((isCCW[0]&&isCCW[1]) || !(isCCW[0]||isCCW[1]))
        {//b1 and b2 have the same signal
            if (vertex[0]->weight < vertex[1]->weight)
            {
                vertex[0]->isInsideCell = false;
            }else{
                vertex[1]->isInsideCell = false;
            }
        }

        //adjust the sequence of the two vorv for each vertex
        for (int i=0; i<2; i++)
        {
            if (!isCCW[i])
            {
                vertex[i]->ClearNeighbor();
                vertex[i]->addVorVertex(vorv1[1]);
                vertex[i]->addVorVertex(vorv1[0]);
            }
        }	
        
    }
    else{//the line has no intersection with the circle boundary
        if (vertex[0]->weight < vertex[1]->weight)
        {
            vertex[0]->type = NOCELL;
        }else{
            vertex[1]->type = NOCELL;
        }
    }	

    for(int i=0; i<2; i++)
    {
        CCVTVertex	*v = vertex[i];
        v->addAdjVertex(vertex[(i+1)%2]);
        v->isAdjBorder = true;
    }	
}

void CPCC2D::InitTwoVertex_RectPoly()
{//check the intersection with the border , if so , add new vorvertex
    assert((int)m_pdt.number_of_faces()==0);

    CCVTVertex * vertex[2];
    VorVertex * vorv, * vorv1, * vorv2;
    
    int index = 0;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit!=vertices.end(); ++vit)
    {
        vertex[index++] = *vit;
    }

    Vector2D	dir;
    Line l;

    Edge_iterator eit = m_pdt.edges_begin();
    for ( ; eit !=m_pdt.edges_end(); ++eit)
    {//should be only one edge
        CGAL::Object o = m_pdt.dual(eit);
        if(CGAL::assign(l,o)){
            dir.x = l.direction().dx();
            dir.y = l.direction().dy();
        }
    }

    // add corner vertex of bounding box to the voronoi vertex list
    if(m_boundType == RECTBOUND){
        AddVoronoiVertex(CORNER, Point(bound_x_min, bound_y_min));
        AddVoronoiVertex(CORNER, Point(bound_x_max, bound_y_min));
        AddVoronoiVertex(CORNER, Point(bound_x_max, bound_y_max));
        AddVoronoiVertex(CORNER, Point(bound_x_min, bound_y_max));
    }else if (m_boundType == POLYBOUND)
    {
        for(int i=0; i<poly->count; i++){
            AddVoronoiVertex(CORNER, poly->acmes[i]); 
        }
    }

    // calculate bisect perpendicular
    Point2D		pa, pb, pts, onpt;	
    int num  = 0;

    if (l.is_degenerate())
    {
        vertex[0]->type = NOCELL;
        vertex[1]->type = NOCELL;
        return;
    }else if (l.is_vertical())
    {		
        onpt.x = -l.c()/l.a();
        onpt.y = 0;
    }else if (l.is_horizontal())
    {
        onpt.x = 0;
        onpt.y = -l.c()/l.b();
    }else{
        onpt.x = l.point().x();
        onpt.y = l.point().y();	
        if (abs(onpt.y)>1e6)
        {
            onpt.y = 1.0;
            onpt.x = -(l.b()+l.c())/l.a();
        }
    }

    int boundedgeNum;
    if(m_boundType == RECTBOUND){
        boundedgeNum = 4;
    }else if (m_boundType == POLYBOUND)
    {
        boundedgeNum = poly->count;
    }

    for(int i=0; i<boundedgeNum; i++)
    {
        if(m_boundType == RECTBOUND){
            switch (i) {
            case 0:
                pa = Point2D(bound_x_min, bound_y_min);
                pb = Point2D(bound_x_max, bound_y_min);
                break;
            case 1:
                pa = Point2D(bound_x_max, bound_y_min);
                pb = Point2D(bound_x_max, bound_y_max);
                break;
            case 2:
                pa = Point2D(bound_x_max, bound_y_max);
                pb = Point2D(bound_x_min, bound_y_max);
                break;
            case 3:
                pa = Point2D(bound_x_min, bound_y_max);
                pb = Point2D(bound_x_min, bound_y_min);
                break;
            }
        }else if (m_boundType == POLYBOUND)
        {
            pa = Point2D(poly->acmes[i].x(), poly->acmes[i].y());
            pb = Point2D(poly->acmes[i+1].x(), poly->acmes[i+1].y());
        }		

        int n = LinesegLineInter(onpt, dir, pa, pb, pts);//for what?
        if ( n < 1 ) continue;
        else
		{
            num ++;
            Point pos(pts.x, pts.y);
            vorv = FindVorVertex(pos);
            if( !vorv )	{
                vorv = new VorVertex;
                vorv->type = INTERSECT;
                vorv->pos  = pos;
                vorv->v[0] = vertex[0];
                vorv->v[1] = vertex[1];
                vorv->bindex = i;
                vorvertices.push_back(vorv);
            }
            assert(vorv != NULL);
            vertex[0]->addVorVertex(vorv);
            vertex[1]->addVorVertex(vorv);
            if (num == 1)
            {
                vorv1 = vorv;
            }else if (num == 2)
            {
                vorv2 = vorv;
            }
        }
    }

    bool onenocell = false;
    if (num == 0)
    {//no intersection with any boundary, set the wpt with smaller weight as hidden point
        if (vertex[0]->weight < vertex[1]->weight)
        {
            vertex[0]->type = NOCELL;
        }else{
            vertex[1]->type = NOCELL;
        }
        onenocell = true;
    }else if(num==1)
    {//one intersection with the boundary
     //very special case, e.g., the line goes across one corner of the boundary
        onenocell = true;
    }

    if (!onenocell)
    {
        bool b1 = IsTriangleCCW(vorv1->pos, vertex[0]->pos, vorv2->pos);
        bool b2 = IsTriangleCCW(vorv1->pos, vertex[1]->pos, vorv2->pos);
        if ((b1&&b2) || !(b1||b2))
        {//b1 and b2 have the same signal
            if (vertex[0]->weight != vertex[1]->weight)
            {
                if (vertex[0]->weight < vertex[1]->weight)
                {
                    vertex[0]->isInsideCell = false;
                }else{
                    vertex[1]->isInsideCell = false;
                }
            }
        }
    }
    
    for(int i=0; i<2; i++)
    {
        CCVTVertex	*v = vertex[i];

        v->addAdjVertex(vertex[(i+1)%2]);
        if (v->type != VALID){
            continue;
        }

        if (onenocell)
        {
            std::list<VorVertex*>::iterator vit;
            for(vit = vorvertices.begin(); vit != vorvertices.end(); vit ++)
            {
                VorVertex * vor = (VorVertex *)*vit;
                v->addVorVertex(vor);
            }
            return;
        }

        v->isAdjBorder = true;

        // if the voronoi diagram around the vertex is intersect with boundary, 
        // we should add more corner voronoi vertex to its adjacent array

        for( int i=0; i<v->GetVorVertexNumber(); i++) {
            VorVertex *v1, *v2, *vnext;
            
            v1 = v->GetVorVertex(i);
            v2 = v->GetVorVertex((i+1)%v->GetVorVertexNumber());

            // need add corner voronoi vertex
            if (v->isInsideCell)
            {
                if (v->adjVorVertex.GetSize() == 2 && IsTriangleCCW(v->pos, v1->pos, v2->pos))
                {
                    // decide which side to insert new voronoi vertex. dmyan 2005-05-10
                    VorVertex * temp;
                    temp = v1;
                    v->adjVorVertex[i] = v2;
                    v1 = v2;
                    v->adjVorVertex[(i+1)%v->GetVorVertexNumber()] = temp;
                    v2 = temp;	
                }
                while( !IsPointOnSameBoundary(v1->pos, v2->pos) ) {
                    vnext = FindNextVorVertexOnBoundaryCCW(v1);
                    v->adjVorVertex.InsertAt(i+1, vnext);
                    i ++;
                    v1 = vnext;
                }
            }else{
                if (v->adjVorVertex.GetSize() == 2 && !IsTriangleCCW(v->pos, v1->pos, v2->pos))
                {
                    // decide which side to insert new voronoi vertex. dmyan 2005-05-10
                    VorVertex * temp;
                    temp = v1;
                    v->adjVorVertex[i] = v2;
                    v1 = v2;
                    v->adjVorVertex[(i+1)%v->GetVorVertexNumber()] = temp;
                    v2 = temp;	
                }
                while( !IsPointOnSameBoundary(v1->pos, v2->pos) ) {
                    vnext = FindNextVorVertexOnBoundaryCCW(v1);
                    v->adjVorVertex.InsertAt(i+1, vnext);
                    i ++;
                    v1 = vnext;
                }
            }
            break;		// dmyan 2005 05 14
        }
    }	
}
 
bool CPCC2D::IfCirclesIntersect()
{//return true if two circles intersect
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = (CCVTVertex*)*vit;

        for(int j=0; j<v->adjVertice.GetSize(); j++) 
        {
            CCVTVertex	*vj = (CCVTVertex *)v->adjVertice.GetAt(j);
            if (If2CirclesTouchEachOther(v,vj))
            {
                return true;
            }
        }
    }
    return false;
}

bool CPCC2D::If2CirclesTouchEachOther(CCVTVertex *v1, CCVTVertex *v2)
{//return true if two circles intersect
    double dist = sqrt(CGAL::squared_distance(v1->pos, v2->pos));
    double rr = sqrt(v1->weight) + sqrt(v2->weight);
    double delta = dist - rr;
    if ( delta < TOLER )
    {//two circles intersect
        return true;
    }
    return false;
}

bool CPCC2D::If1CircleTouchesBoundary(CCVTVertex* v)
{//return true if the circle touches the boundary
    if (m_boundType == RECTBOUND)
    {
        Point pa, pb;
        double r = sqrt(v->weight);
        for(int i=0; i<4; i++){
            switch (i) {
            case 0:
                pa = Point(bound_x_min, bound_y_min);
                pb = Point(bound_x_max, bound_y_min);
                break;
            case 1:
                pa = Point(bound_x_max, bound_y_min);
                pb = Point(bound_x_max, bound_y_max);
                break;
            case 2:
                pa = Point(bound_x_max, bound_y_max);
                pb = Point(bound_x_min, bound_y_max);
                break;
            case 3:
                pa = Point(bound_x_min, bound_y_max);
                pb = Point(bound_x_min, bound_y_min);
                break;
            }
            
            double tmpx, tmpy;
            double dist = disPt2Line(pa.x(),pa.y(),pb.x(),pb.y(),v->pos.x(),v->pos.y(), tmpx, tmpy);
            if ((dist - TOLER)< r)
            {//intersects
                return true;
            }
        }//end for i=0-4
    }else if (m_boundType == CIRCLEBOUND)
    {
        Point cent = Point(center_x, center_y);
        double dist = sqrt(CGAL::squared_distance(v->pos, cent));
        double vr = sqrt(v->weight);
        //double delta_internal = dist - (radius_r - vr);
        //double delta_external = dist - (radius_r + vr);		
        //if((abs(delta_internal) < TOLERLOW) ||(abs(delta_external) < TOLERLOW))
        if ((dist+/*TOLERLOW*/TOLER)>(radius_r - vr))
        {
            return true;
        }
    }
    
    return false;
}

bool CPCC2D::IfCirclesTouchBoundary()
{//return true if some circle touches the boundary
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = (CCVTVertex*)*vit;
        if (If1CircleTouchesBoundary(v))
        {
            return true;
        }
    }
    return false;
}

void CPCC2D::FindCentroid(CCVTVertex* v, Point& pt)
{
    double x = 0, y =0;
    Kernel::FT	w(0.0);

    for(int i=0; i<v->adjVorVertex.GetSize(); i++)
    {		
        VorVertex * vv1, * vv2;
        vv1 = (VorVertex *)v->adjVorVertex.GetAt(i);
        vv2 = (VorVertex *)v->adjVorVertex.GetAt( (i+1)%v->adjVorVertex.GetSize());
        
        Kernel::FT a = Triangle(vv1->pos, vv2->pos, v->vh->point()).area();
        Point c = Point(a*(vv1->pos.x() + vv2->pos.x() + v->pos.x())/3, a*(vv1->pos.y() + vv2->pos.y() + v->pos.y())/3);
        x += c.x(); y += c.y();
        w += a;
    }
    x /= w; y /= w;
    pt= Point(x, y);
}

double CPCC2D::FindMaxInscribedCircle_Polygon(CCVTVertex* v, Point& pt)
{
    const int adjVorVNum = (int)v->adjVorVertex.GetSize();
    if (adjVorVNum < 2){
        return -1;
    }
    //check whether the voronoi cell is convex

    std::vector<Point_2> vorcellnodes;
    for(int i=0; i<adjVorVNum; i++)
    {
        VorVertex * vv = (VorVertex *)v->adjVorVertex.GetAt(i);
        vorcellnodes.push_back(Point_2(vv->pos.x(), vv->pos.y()));
    }
    std::vector< Point_2 > convexPts;
    CGAL::ch_graham_andrew( vorcellnodes.begin(), vorcellnodes.end(), std::back_inserter(convexPts));

    CPolygon* v_poly = new CPolygon();
    v_poly->Init(convexPts);
    double res = v_poly->MaxInscribedCircle(pt);

    if (v_poly != NULL)  delete v_poly;
    return res;
}

struct PositionDis
{	
    Vector2D pos;	
    double dis;
};
bool PositionDis_comp_smaller(PositionDis a, PositionDis b)
{
    return a.dis<b.dis;
};
int CPCC2D::FindMinDistanceEdge_Polygon(CCVTVertex* v, Point& pt, Vector2D& minNor, double& mindis)
{
    const int adjVorVNum = (int)v->adjVorVertex.GetSize();
    if (adjVorVNum < 2){
        return -1;
    }
    
    double x0, y0, minx0, miny0;
    mindis = 10e10;
    int minid = -1;
    std::vector<Vector2D> allDirs;
    std::vector<PositionDis> allPos;
    for(int i=0; i<adjVorVNum; i++)
    {
        VorVertex * vv1, * vv2;
        vv1 = (VorVertex *)v->adjVorVertex.GetAt(i);
        vv2 = (VorVertex *)v->adjVorVertex.GetAt( (i+1)%adjVorVNum);

        //double dis = disPt2Line(vv1->pos.x(),vv1->pos.y(),vv2->pos.x(),vv2->pos.y(),v->pos.x(),v->pos.y(),x0,y0);
        double dis = signedDisPt2Line(vv1->pos.x(),vv1->pos.y(),vv2->pos.x(),vv2->pos.y(),v->pos.x(),v->pos.y(),x0,y0);
        if (dis<0){
            continue;
        }
        PositionDis tmppd;
        tmppd.pos = Vector2D(x0 - v->pos.x(), y0 - v->pos.y());
        tmppd.dis = dis;
        allPos.push_back(tmppd);

        if (dis < mindis){
            mindis = dis;
            minid = i;
            minx0 = x0;
            miny0 = y0;
        }
    }

    std::sort(allPos.begin(), allPos.end(), PositionDis_comp_smaller);
    for (int i=0; i<(int)allPos.size(); i++){
        if (fabs(allPos[i].dis-mindis) < TOLERLOW2){
            allDirs.push_back(allPos[i].pos);
        }
    }
    Vector2D leastnormV;
    double leastnorm = 10e10;
    int n = (int)allDirs.size();
    for (int i=0; i<n; i++)
    {
        Vector2D v0 = allDirs[i];
        Vector2D v1 = allDirs[(i+1)%n];
        double alpha = (v1*v1-v0*v1)/((v0-v1)*(v0-v1));
        Vector2D tmpv = alpha*v0 + (1-alpha)*v1;
        double tmpd = tmpv*tmpv;
        if (tmpd<leastnorm){
            leastnorm = tmpd;
            leastnormV = tmpv;
        }
    }
    if ((int)allDirs.size()>1){
        minNor = leastnormV;
    }else{
        minNor = Vector2D(minx0 - v->pos.x(), miny0 - v->pos.y());
    }	
    pt = Point(minx0, miny0);
    return minid;
}

void CPCC2D::FreeMemory()
{
    FreeCVTVertices();
    FreeCVTFaces();
    FreeVorVertices();
}

void CPCC2D::FreeCVTVertices()
{
    for(std::list<CCVTVertex *>::iterator vit=vertices.begin(); vit!=vertices.end(); vit++)
    {
        CCVTVertex * v = (CCVTVertex *)*vit;
        delete v;
    }
    vertices.clear();
}

void CPCC2D::FreeCVTFaces()
{
    for(std::list<CCVTFace*>::iterator fit = faces.begin(); fit!=faces.end(); fit++)
    {
        CCVTFace * f = (CCVTFace *)*fit;
        delete f;
    }
    faces.clear();
}

void CPCC2D::FreeVorVertices()
{
    for(std::list<VorVertex*>::iterator vvit = vorvertices.begin(); vvit!=vorvertices.end(); vvit++)
    {
        VorVertex * vv = (VorVertex *)*vvit;
        delete vv;
    }
    vorvertices.clear();
}

double CPCC2D::EvaluateErrorFuncIncircle()
{
    double	totalerr = 0;

    //for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit ++) 
    //{
    //	CCVTVertex * v = *vit;
    //	//totalerr += v->inRadius;
    //	totalerr += abs(v->inRadius - sqrt(v->weight));
    //}

    //f = 1/k
    totalerr = 1.0/m_k;
    return totalerr;
}

double CPCC2D::EvaluateErrorFunc()
{
    double	totalerr = 0, trierr, vorerr;	

    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit ++) 
    {
        CCVTVertex * _v = *vit;
        int n = _v->GetVorVertexNumber();
        double w = _v->weight;
        double area = 0;
        VorVertex * P1, *P2;
            
        vorerr = 0;
        for(int i=0; i<n; i++){
            P1 = _v->GetVorVertex(i);
            P2 = _v->GetVorVertex((i+1)%n);
            Vector2D V1(P1->x() - _v->x(), P1->y() - _v->y());
            Vector2D V2(P2->x() - _v->x(), P2->y() - _v->y());
            double area_qua = fabs(V1^V2);
            trierr =  area_qua* (V1*V1 + V1*V2 + V2*V2)/ 12;
            area += area_qua;
            vorerr += trierr;
        }

        totalerr += vorerr;
        //totalerr -= 0.5*w*area; //??what is this for???
    }
    return totalerr;
}
    
double CPCC2D::EvaluateErrorDist() 
{
    double	totalerror = 0.0;

    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = *vit;

        if (v->type != VALID){			
            continue;
        }		

        double error = CGAL::squared_distance(v->pos, v->centroid);
        totalerror += error;
    }

    return totalerror;
}

CCVTVertex* CPCC2D::FindMyVertex(Vertex_handle vhandle)
{
    for(std::list<CCVTVertex *>::iterator vit=vertices.begin(); vit!=vertices.end(); vit++)
    {
        CCVTVertex * v = (CCVTVertex *)*vit;
        if( vhandle == v->vh )
            return v;
    }
    return NULL;
}

CCVTVertex* CPCC2D::FindMyVertex(Point pt)
{
    for(std::list<CCVTVertex *>::iterator vit=vertices.begin(); vit!=vertices.end(); vit++)
    {
        CCVTVertex * v = (CCVTVertex *)*vit;
        if( IsSamePoint(pt, v->vh->point()) )
            return v;
    }
    return NULL;
}


CCVTFace*	CPCC2D::FindMyFace(Face_handle fhandle)
{
    for(std::list<CCVTFace*>::iterator fit = faces.begin(); fit!=faces.end(); fit++)
    {
        CCVTFace * f = (CCVTFace *)*fit;
        if( fhandle == f->fh )
            return f;
    }
    return NULL;
}

VorVertex*	CPCC2D::FindVorVertex(Point pt)
{
    std::list<VorVertex*>::iterator vit;
    for(vit = vorvertices.begin(); vit != vorvertices.end(); vit ++)
    {
        VorVertex * v = (VorVertex *)*vit;
        if ( IsSamePoint(v->pos, pt) )
            return v;
    }
    return NULL;
}

void CPCC2D::AddVoronoiVertex(VertexType type, Point pos)
{
    VorVertex * vv = new VorVertex;
    vv->type = type;
    vv->pos = pos;	
    vorvertices.push_back(vv);
}

bool CPCC2D::IsPointOnSameBoundary(Point p1, Point p2)
{
    if (m_boundType == RECTBOUND){
        if( fabs( p1.x() - bound_x_max ) < TOLERLOW && fabs(p2.x() - bound_x_max) < TOLERLOW)
            return true;
        else if( fabs( p1.x() - bound_x_min ) < TOLERLOW && fabs(p2.x() - bound_x_min) < TOLERLOW)
            return true;
        else if( fabs( p1.y() - bound_y_max ) < TOLERLOW && fabs(p2.y() - bound_y_max) < TOLERLOW)
            return true;
        else if( fabs( p1.y() - bound_y_min ) < TOLERLOW && fabs(p2.y() - bound_y_min) < TOLERLOW)
            return true;
    }else if (m_boundType == POLYBOUND){
        for(int i=0; i<poly->count; i++)
        {
            if (IsPtOnSeg(p1, poly->acmes[i], poly->acmes[i+1]) && IsPtOnSeg(p2, poly->acmes[i], poly->acmes[i+1]))
                return true;
        } 
    }

    return false;
}

bool CPCC2D::IsTriangleCCW(Point p1, Point p2, Point p3)
{
    Triangle t(p1, p2, p3);
    if ( t.orientation() > 0 )
        return true;
    return false;
}

VorVertex*  CPCC2D::FindNextVorVertexOnBoundaryCCW(VorVertex * vcur)
{
    VorVertex * vnext = NULL;

    if (m_boundType == RECTBOUND)
    {
        // if current is a corner vertex
        if (IsSamePoint(Point(bound_x_min, bound_y_min), vcur->pos))
            vnext = FindVorVertex(Point(bound_x_max, bound_y_min) );
        if (IsSamePoint(Point(bound_x_max, bound_y_min), vcur->pos))
            vnext = FindVorVertex(Point(bound_x_max, bound_y_max) );
        if (IsSamePoint(Point(bound_x_max, bound_y_max), vcur->pos))
            vnext = FindVorVertex(Point(bound_x_min, bound_y_max) );
        if (IsSamePoint(Point(bound_x_min, bound_y_max), vcur->pos))
            vnext = FindVorVertex(Point(bound_x_min, bound_y_min) );

        if( vnext )	return vnext;

        // if current lies on boundary
        if ( fabs(vcur->pos.y() - bound_y_min) < TOLERLOW )
            vnext = FindVorVertex(Point(bound_x_max, bound_y_min) );
        if ( fabs(vcur->pos.y() - bound_y_max) < TOLERLOW )
            vnext = FindVorVertex(Point(bound_x_min, bound_y_max) );
        if ( fabs(vcur->pos.x() - bound_x_min) < TOLERLOW )
            vnext = FindVorVertex(Point(bound_x_min, bound_y_min) );
        if ( fabs(vcur->pos.x() - bound_x_max) < TOLERLOW )
            vnext = FindVorVertex(Point(bound_x_max, bound_y_max) );		
    }
    else if (m_boundType == POLYBOUND)
    {
        // if current is a corner vertex
        for(int i=0; i<poly->count; i++)
        {
            if (IsSamePoint(poly->acmes[i], vcur->pos)){
                vnext = FindVorVertex(poly->acmes[i+1]);
                break;
            }
        }
        if( vnext )	return vnext;

        // if current lies on boundary
        for(int i=0; i<poly->count; i++)
        {
            int j = i+1;
            double tmp = fabs((poly->acmes[j].y()-poly->acmes[i].y())*vcur->pos.x() - 
                (poly->acmes[j].x()-poly->acmes[i].x())*vcur->pos.y() - 
                poly->acmes[j].y()*poly->acmes[i].x()+poly->acmes[j].x()*poly->acmes[i].y());
            if (tmp < TOLERLOW)
            {
                vnext = FindVorVertex(poly->acmes[j]);
                break;
            }
        }
    }
    assert(vnext != NULL);
    return vnext;
}

bool CPCC2D::IsSameDualEdge(VorVertex* v1, VorVertex *v2)
{
    if( (v1->v[0] == v2->v[0] && v1->v[1] == v2->v[1]) 
        || (v1->v[0] == v2->v[1] && v1->v[1] == v2->v[0] ) )
        return true;
    return false;
}


//not used actually
void CPCC2D::GeneratePDVertex()	// for each non boundary vertex, create its one ring PD vertex.
{	
    bool		isboundary = false;
    int			index = 0;
    m_points.clear();

    if(m_pdt.number_of_vertices() < 1) return;

    // find all non-boundary vertex to be moved
    for( PDtri::Finite_vertices_iterator vit = m_pdt.finite_vertices_begin();
        vit != m_pdt.finite_vertices_end();  ++vit)
    {
        isboundary = false;
        /* */
        //LLoyd  Centre de gravite de la cellule de Voronoi
        PDtri::Edge_circulator ec=m_pdt.incident_edges(vit),
            done=ec;
        do{ 
            Segment s;
            if ( ( m_pdt.is_infinite( ec) ) || ( ! CGAL::assign( s, m_pdt.dual(ec)) ) )
            { 
                isboundary = true; 
                break;  
            }
            ++ec ; 
        }
        while( ec != done);

        if( !isboundary )
        {
            CCVTVertex * v = new CCVTVertex();
            v->id = index ++;
            v->vh = vit;
            vertices.push_back(v);
        }
        else m_points.push_back(vit->point());
    }

    // 
    index = 0;
    for( std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit ++)
    {
        CCVTVertex * v = *vit;
        if (v->type != VALID)
        {
            continue;
        }
        PDtri::Face_circulator fc = m_pdt.incident_faces(v->vh), done=fc;
        do {
            CCVTFace *f = FindMyFace(fc);
            if( !f ) {

                // add a new face
                f = new CCVTFace();
                f->id = index ++;
                f->fh = fc;
                Face fd = *fc;
                for( int i=0; i<3; i++)	{
                    f->v[i] = FindMyVertex(fd.vertex(i));
                }			
                faces.push_back(f);
            }

            // add face's dual voronoi vertex
            VorVertex * vorv;
            Point pt = m_pdt.weighted_circumcenter(fc);
            vorv = FindVorVertex(pt);
            if ( !vorv )
            {
                vorv = new VorVertex;
                vorv->type = CIRCUMCENTER;
                vorv->pos  = pt;				
                vorv->face = f;
                vorvertices.push_back(vorv);
            }
            // add voronoi vertex to current vertex
            v->adjVorVertex.Add(vorv);

            ++fc;
        } while(fc != done);
    }
}

void CPCC2D::wcvtCenter(Point p1, Point p2, Point p3, double& x, double& y)
{//used for weighted CVT
    //density function p=exp(-5x-5y)
    double x1=p1.x();
    double y1=p1.y();
    double x2=p2.x();
    double y2=p2.y();
    double x3=p3.x();
    double y3=p3.y();
    double x12 = x1*x1;
    x = (4 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y3 + 0.5e1 * x2 * pow(x3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - 0.10e2 * x2 * x2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - 0.10e2 * pow(x2, 0.3e1) * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 + 0.10e2 * pow(x2, 0.3e1) * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 - 0.5e1 * x2 * pow(y1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 - 0.5e1 * x2 * pow(x3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 + 0.15e2 * x2 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y3 + 0.10e2 * pow(x2, 0.3e1) * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 + 0.15e2 * x2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * x1 - 0.10e2 * x2 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 + 0.15e2 * x2 * x2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 + 0.5e1 * x2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 - 2 * x2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - 0.5e1 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 + 0.10e2 * pow(x2, 0.3e1) * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 - 0.5e1 * x2 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 - 0.5e1 * x2 * pow(x3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 - 0.5e1 * x2 * x2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 - 4 * x2 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 - 0.5e1 * x2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 - 0.10e2 * pow(x2, 0.3e1) * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 + 0.5e1 * x2 * pow(x1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - 0.5e1 * x2 * pow(x1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 - 0.5e1 * x2 * x2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 - 0.10e2 * x2 * x2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - 0.5e1 * x1 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x3 + 0.10e2 * x1 * x1 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 - 0.5e1 * x2 * x2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 + 0.10e2 * pow(x2, 0.3e1) * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 - 0.5e1 * x2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 - 0.5e1 * x2 * pow(x1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 - 0.5e1 * x2 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 - 4 * x1 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 - x1 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - 0.5e1 * x1 * pow(y3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 + 4 * x1 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x2 + 4 * x1 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x2 + 4 * x2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 - 2 * x1 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 + x1 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 + 0.10e2 * x1 * x1 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 - 0.5e1 * x1 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x2 - 0.5e1 * x2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x3 - 0.5e1 * x2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y3 + 4 * x2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 + 0.10e2 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x3 - x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y3 + 0.10e2 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y3 - 0.3e1 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x3 - 2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 - 2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x3 - 2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * y2 - x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * x3 - 0.3e1 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * y3 + 6 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * x3 + 0.5e1 * x2 * x2 * pow(x1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) - 2 * x2 * x2 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) + x2 * pow(y3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) - x1 * pow(y3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) + pow(x1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - x3 * pow(y1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) - 0.5e1 * pow(x2, 0.3e1) * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) - 2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x3 * x2 * x2 + 4 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x3 * x1 * x1 + 2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x3 * y1 * y1 - 0.5e1 * pow(x2, 0.3e1) * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) - 2 * x2 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) + x2 * pow(y1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) - 2 * pow(x1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 - 0.5e1 * pow(x2, 0.3e1) * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) + 6 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x3 + 2 * pow(x1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 - 2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x2 - pow(x1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 + 0.5e1 * x2 * x2 * pow(y3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) + 0.5e1 * x2 * x2 * pow(x3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) + 2 * x2 * pow(x3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) - 0.5e1 * pow(x2, 0.3e1) * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) - 2 * x1 * pow(x3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) + pow(x3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - pow(x3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 + 2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y3 + 0.5e1 * x2 * x2 * pow(y1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) + 0.5e1 * x2 * pow(y3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - 0.5e1 * x2 * pow(y3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 + 2 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y3 - 2 * x2 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - 4 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x2 + 0.20e2 * x2 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * y2 - 0.15e2 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y3 * y3 - 4 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 + 0.15e2 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 * y3 - 0.5e1 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y2 - 0.15e2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * x1 * x3 + 0.20e2 * x2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 + exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * y1 * y1 * y2 - 0.5e1 * x2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * y2 - 0.20e2 * x2 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 - 0.15e2 * x2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y3 + 0.15e2 * x2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 - 0.10e2 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y2 * y2 - 0.10e2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 * x1 + 0.10e2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 * x3 - 0.10e2 * x2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * x3 - 0.15e2 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * x1 + 2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * y1 * y1 * y3 - 0.5e1 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * x3 + 0.15e2 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * x1 + 0.10e2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y2 * y2 - 0.5e1 * x2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * y2 - 0.15e2 * x1 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x3 * x3 + 0.10e2 * x1 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 * x2 - 0.15e2 * x1 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x3 + 0.10e2 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 * x1 + 0.20e2 * x3 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * x2 - 0.10e2 * x3 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x2 * y3 - 2 * x1 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x3 - 0.10e2 * x1 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x2 * y1 + 0.20e2 * x2 * x2 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 - 0.5e1 * x2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 - 0.5e1 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * x3 * x3 + 4 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 + 0.20e2 * x1 * x1 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * y3 + 4 * x1 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * y2 - 0.5e1 * x1 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * x2 + 0.20e2 * x2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * y2 - 0.5e1 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 * y3 - 2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x3 * y2 * y1 - 0.5e1 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * x3 - 0.10e2 * x1 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x2 * y3 + 0.20e2 * x1 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * x2 - x2 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 + 0.20e2 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x3 * x3 + 4 * x1 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * y2 - 2 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 - 2 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x3 + 8 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * y1 * y3 * x1 - 0.20e2 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x2 * y2 - 4 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * y2 - x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 * y3 - 4 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * y3 - 2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * y3 * y1 * y2 + 4 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * y2 * y1 * x2 - 0.15e2 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * x3 - 4 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * y1 * y3 * x2 + 0.5e1 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * x2 + 0.15e2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * x1 * y2 - 0.15e2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * x1 * y3 - 4 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 * x3 + 4 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 * x2 - 0.10e2 * x2 * x2 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 + x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 * y2 - 2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y3 + 2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y1 - 0.5e1 * x2 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 - 2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y3 + x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 * y3 - x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y3 * y3 + 2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 + 0.15e2 * x2 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 - 0.5e1 * x2 * pow(y1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 + 0.5e1 * x2 * pow(y1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - 2 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 + exp(-0.5e1 * x3 - 0.5e1 * y3) * x3 * pow(y2, 0.3e1) + 0.5e1 * x1 * x1 * pow(y3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) + 4 * x2 * x2 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) + x3 * pow(y1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) - 0.5e1 * pow(x3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 + 0.5e1 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(x2, 0.3e1) - 0.5e1 * pow(x3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * x1 + 0.5e1 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(y2, 0.3e1) - 0.5e1 * pow(x3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 + 2 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 2 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x1 * x1 * y2 * y2 - 2 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x1 * x1 * x2 * x2 + exp(-0.5e1 * x1 - 0.5e1 * y1) * x1 * pow(y2, 0.3e1) - 2 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x1 * x1 * y3 * y3 - 2 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x1 * x1 * x3 * x3 - 0.5e1 * pow(x3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 + pow(x1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y3 - pow(x3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(y2, 0.3e1) + x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y3 * y2 * y2 - 2 * x1 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y1 - pow(x1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 2 * x1 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y3 + 2 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 * y2 - x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(y2, 0.3e1) + x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 * y3 + 2 * x1 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y3 + x1 * pow(y3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) - 0.5e1 * x3 * x3 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 * y2 + pow(x3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 + 0.10e2 * x3 * y1 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - 0.5e1 * x3 * x3 * y1 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 + 0.5e1 * x1 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(x2, 0.3e1) - 2 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 - 0.10e2 * pow(x1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y2 + 0.5e1 * x1 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(y2, 0.3e1) - 2 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 - 2 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 0.5e1 * pow(x1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 - x2 * pow(y3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) - 2 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * x3 + 4 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x3 * x3 * x1 * y2 + 4 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x3 * x3 * y2 * y1 - 2 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x3 * x3 * y2 * y2 + 0.10e2 * pow(x1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 + 0.10e2 * pow(x1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x3 - 2 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x3 * x3 * x1 * x1 + 0.10e2 * pow(x1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * y2 - 0.10e2 * pow(x1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * x3 - 0.5e1 * pow(x1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * x3 + 0.15e2 * x1 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 * x2 - 4 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x3 * x3 * x1 * y1 - 0.5e1 * pow(x1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 + 0.10e2 * pow(x1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y3 + 0.15e2 * x1 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 * y2 - 0.5e1 * pow(x1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 - 4 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x3 * x3 * x2 * y2 + 4 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x3 * x3 * x2 * y1 - 2 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x3 * x3 * x2 * x2 + 0.15e2 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 * y2 + 0.10e2 * pow(x3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y1 - pow(x2, 0.3e1) * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) - 0.10e2 * pow(x3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * y1 - 0.10e2 * pow(x3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y2 - x2 * pow(y1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) + 4 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x3 * x3 * x2 * x1 + 0.10e2 * pow(x3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y1 - 2 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x3 * x3 * y1 * y1 + 4 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x1 * x1 * y2 * y3 - 2 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y1 + 2 * x1 * pow(x3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) + 4 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x1 * x1 * x3 * y2 + 4 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x1 * x1 * x2 * x3 + 0.5e1 * pow(x1, 0.3e1) * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) + 4 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x1 * x1 * x2 * y3 - 4 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x1 * x1 * y3 * x3 - pow(x2, 0.3e1) * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) + 0.10e2 * pow(x3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x1 - 4 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x1 * x1 * x2 * y2 + pow(x2, 0.3e1) * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) + 0.15e2 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 * x2 + 0.10e2 * pow(x3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * y2 + 2 * x2 * x2 * y1 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) + 0.5e1 * x1 * x1 * pow(x3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) + 2 * x1 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 + x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y2 * y2 + 2 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y2 * y2 + 0.5e1 * x3 * x3 * pow(y1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) + 2 * x2 * x2 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) - 2 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 + 2 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y3 + pow(x2, 0.3e1) * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) + x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y3 * y3 - x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 * y3 - 2 * x2 * pow(x3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) + 2 * pow(x1, 0.3e1) * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) - 2 * pow(x2, 0.3e1) * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) + 0.10e2 * x1 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 - 0.10e2 * x1 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y3 + 0.5e1 * pow(x1, 0.3e1) * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y3 - 0.5e1 * pow(x1, 0.3e1) * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 2 * pow(x2, 0.3e1) * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) - 2 * x2 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x3 - 2 * x2 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * x3 + 0.40e2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 * x2 - 0.5e1 * x1 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 + 0.15e2 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y3 + 0.15e2 * x1 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 - 0.5e1 * x1 * x1 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 + 2 * pow(x2, 0.3e1) * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) - 0.5e1 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(x2, 0.3e1) - 0.5e1 * pow(x1, 0.3e1) * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 - 2 * x2 * pow(x1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) - 0.5e1 * x1 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 + 0.10e2 * x1 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 + 2 * pow(x2, 0.3e1) * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) - 2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 - 2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 + 4 * x2 * x2 * x1 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) - 0.5e1 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y1 - 0.15e2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 * x2 - 0.10e2 * x1 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y2 - 0.15e2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 * y3 + 0.20e2 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y2 - 0.15e2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 * x2 + 0.15e2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y3 * y3 + 0.20e2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 * y3 - 0.15e2 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y3 - 0.15e2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 * y2 - 0.5e1 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 * y1 + x2 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - 0.15e2 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 - 0.5e1 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y1 - 0.10e2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y3 * x2 - x2 * x2 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 + 2 * x2 * x2 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 0.5e1 * x2 * x2 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 0.5e1 * x2 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 + 6 * x2 * x2 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 - 0.10e2 * x1 * x1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y2 - 0.10e2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y1 * x2 - 0.5e1 * x1 * x1 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 - 0.10e2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y1 * y2 + x2 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 - x2 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * x1 + 0.20e2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y2 * y2 - 4 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 - 0.10e2 * x1 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 - 0.5e1 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y2 * y2 + 0.5e1 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 - 2 * x2 * x2 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 - 4 * x2 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 2 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y2 - 2 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y1 - 0.5e1 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * x3 - 0.3e1 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 4 * x2 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 - 4 * x2 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 + 0.5e1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(x2, 0.3e1) - 0.5e1 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 + 4 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 - 0.5e1 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 0.5e1 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 - 0.10e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 0.5e1 * x1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(x2, 0.3e1) + 0.10e2 * x1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 - 0.10e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 - 0.10e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * x3 - 0.10e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 - 0.5e1 * x1 * x1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 + 0.10e2 * x1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 + 0.5e1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(y2, 0.3e1) - 0.20e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * y1 + 0.15e2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 * x2 - 0.10e2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * y1 * x2 + 0.5e1 * x1 * pow(y3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 - 0.5e1 * x1 * pow(y3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 - 0.5e1 * x1 * x1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 - 0.5e1 * x1 * x1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 + 0.20e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y1 + 0.10e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x1 + 8 * x2 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 * y1 - 0.15e2 * x1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 * x2 + 0.15e2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 * y2 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 * x1 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x1 * x1 - 0.10e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * x1 - 0.10e2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y1 * x2 - 0.5e1 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 + 0.20e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * y2 - 0.10e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * x2 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 * x1 + 0.10e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 - 0.10e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * x3 + 0.20e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y3 - 0.10e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x2 - 0.20e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * x3 + 0.20e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * y2 + 0.15e2 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 * y2 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * x1 * y2 + 0.10e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y2 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 * y2 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 * x2 - 0.10e2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x1 * y2 - 0.10e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 + 0.20e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y1 - 0.20e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y2 - 0.10e2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y1 * x1 - 0.10e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 + 0.10e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * x2 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 * y2 - 0.15e2 * x1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 * y2 + 0.10e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * x1 - 0.5e1 * x1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y2 + 0.20e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x1 - 0.10e2 * x1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y1 * x2 - 0.5e1 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - 0.5e1 * x1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * x2 + 0.20e2 * x1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x2 + 0.10e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * x2 - 0.5e1 * x1 * x1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 - 0.5e1 * x1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(y2, 0.3e1) - 0.5e1 * x1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y2 * y2 - 4 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 + 0.10e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x3 - 0.5e1 * x1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 * y1 + 0.15e2 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 * x2 - 0.20e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y2 - 0.10e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * x1 + 0.20e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 + 0.20e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x3 + 0.10e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x3 + 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(x2, 0.3e1) - 0.5e1 * x1 * pow(y3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 0.5e1 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * x1 - 0.10e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 - 0.5e1 * x2 * x1 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 + 6 * x2 * x1 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 + 6 * x2 * x2 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 + 8 * x2 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 - 0.5e1 * x2 * x2 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 0.5e1 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 + 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(y2, 0.3e1) - 0.10e2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x1 * y2 - 4 * x2 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 - 2 * x1 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 - 4 * x2 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - 2 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 - 0.5e1 * x1 * pow(x3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 + 0.5e1 * x2 * x2 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 0.3e1 * x1 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 0.5e1 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(y2, 0.3e1) - 2 * x2 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 - 0.15e2 * x1 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y1 - 0.10e2 * x1 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y1 - 0.5e1 * x1 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 - 0.5e1 * x1 * pow(x3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 - 0.15e2 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 * y2 - 0.15e2 * x1 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * x2 + 0.5e1 * x1 * pow(x3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 - 0.10e2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 * x1 + 0.20e2 * x1 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y2 - 0.15e2 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 * x2 - 0.10e2 * x1 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y2 + 0.40e2 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y1 * x2 + 0.15e2 * x1 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 + 0.20e2 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 * y1 - 0.15e2 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 * x2 - 0.15e2 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 * y2 - 0.10e2 * x1 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y1 + 0.15e2 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 * y3 + 0.15e2 * x1 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y3 - 0.10e2 * x3 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 * x2 - 0.10e2 * x1 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x3 * y2 + 0.20e2 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 * y2 + 4 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 + 0.15e2 * x1 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 - 0.5e1 * x1 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 + 0.40e2 * x1 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x3 * y1 + 4 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x1 + 4 * x2 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 2 * x2 * x2 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 + 0.10e2 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 - 0.5e1 * x1 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 - 2 * x2 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * x3 + 0.10e2 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 0.5e1 * x1 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 0.5e1 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(x2, 0.3e1) - 4 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y2 + 4 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y1 - 4 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * y1 + 4 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * y2 - 0.3e1 * x2 * x2 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 + 4 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 + 0.5e1 * x1 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 - 2 * x2 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 + 4 * x2 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - 0.5e1 * x3 * x3 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 - 0.10e2 * x3 * x3 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y2 - 0.3e1 * x2 * x2 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 + 0.5e1 * x2 * x2 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 0.15e2 * x3 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 * y2 + 0.20e2 * x3 * y1 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * x2 - 0.15e2 * x3 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 * x2 + 4 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 + 4 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x3 - 2 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y3 * y1 * y2 + 4 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y3 + 0.15e2 * x1 * x1 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 + 2 * x2 * y1 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 + 6 * x2 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 0.5e1 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(y2, 0.3e1) - x2 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - x3 * y1 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 2 * x2 * x2 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 4 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * x3 - 4 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y2 + 4 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x3 - 2 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y1 * y2 - x3 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - x1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * x1 - 2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 - x2 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * x3 + 2 * x2 * x2 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - x1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 + x2 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 - 2 * x2 * x2 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 + x2 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - x2 * x2 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 + 2 * x2 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 4 * x2 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 - 0.5e1 * x3 * pow(y1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - x2 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 0.5e1 * x3 * x3 * y1 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 + 0.5e1 * x3 * pow(y1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y3 - 0.5e1 * x3 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(x2, 0.3e1) - 0.5e1 * x3 * pow(y1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 - 2 * x2 * x2 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 0.5e1 * x3 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(y2, 0.3e1) + 0.10e2 * x3 * y1 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2) / (x1 - x3 + y1 - y3) / (x1 - x2 + y1 - y2) / (x2 - x3 + y2 - y3) / (-exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 + x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) - exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 + exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 - exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 + exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 - exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 + exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 - exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 + exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 + exp(-0.5e1 * x2 - 0.5e1 * y2) * y1) / 0.5e1;
    y =  -(-4 * y3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * y2 - 0.5e1 * y2 * pow(y3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 - 0.5e1 * y2 * pow(y3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 + 2 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y3 + 0.5e1 * x2 * pow(x3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 + 0.15e2 * y2 * y2 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 - 0.5e1 * x2 * x2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 + 0.15e2 * y2 * y2 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x3 - 0.10e2 * y1 * pow(y2, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 - 0.5e1 * y1 * y1 * y2 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 + 0.10e2 * y1 * pow(y2, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 + 0.15e2 * y1 * y1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * x3 - 0.5e1 * y1 * y1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 - 4 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * x2 * x3 * y2 - 0.5e1 * x2 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 + x2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 - 2 * x2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 + 0.20e2 * y1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 * x3 - 0.10e2 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 + 0.20e2 * y1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 * x3 + 0.10e2 * pow(y1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x3 + 0.15e2 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 * x2 - 0.10e2 * pow(y1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x2 + 0.15e2 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 * y2 - 0.10e2 * x2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 + 0.5e1 * x2 * pow(x1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - 0.5e1 * x2 * x2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - x1 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x3 - 0.15e2 * y1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x3 * x3 - 0.15e2 * y1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x3 - 0.10e2 * x2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 - 2 * x1 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 - x1 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - 0.3e1 * x1 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - 0.10e2 * pow(y1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * x3 + 0.10e2 * pow(y1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y2 + 0.10e2 * pow(y1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * x2 + 0.10e2 * pow(y1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x3 + 0.10e2 * y1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 0.5e1 * y1 * y1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 2 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x1 - 0.5e1 * x3 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(x2, 0.3e1) + 0.10e2 * pow(y3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x1 - x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y3 + x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x3 - 0.15e2 * x3 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 * x2 - 2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * y2 - x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * x3 + x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * y3 + 2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * x3 - 0.10e2 * y1 * y1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x3 + 0.15e2 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 * x2 + 0.10e2 * x3 * x3 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 + 0.10e2 * x3 * x3 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 - 0.5e1 * pow(x3, 0.3e1) * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 0.5e1 * x3 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 0.5e1 * pow(x3, 0.3e1) * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 + x2 * pow(y3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) - x1 * pow(y3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) + pow(x1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - x3 * pow(y1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) + 2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x3 * y1 * y1 + x2 * pow(y1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) + 2 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x3 - pow(x1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 + pow(x3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - pow(x3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 + 2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y3 + 0.5e1 * x2 * pow(y3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 + 6 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y3 - 2 * x2 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 + 0.10e2 * x2 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * y2 - 4 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 + 0.15e2 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 * y3 - 0.5e1 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y2 + 0.10e2 * x2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 - 0.3e1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * y1 * y1 * y2 - 0.5e1 * x2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * y2 - 0.10e2 * x2 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 + 0.15e2 * x2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 - 0.20e2 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y2 * y2 - 0.20e2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 * x1 + 0.20e2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 * x3 + 6 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * y1 * y1 * y3 - 0.5e1 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * x3 + 0.15e2 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * x1 + 0.20e2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y2 * y2 - 0.5e1 * x2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * y2 + 0.20e2 * x1 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 * x2 + 0.20e2 * x2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 * x1 - 2 * x1 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x3 + 0.10e2 * x2 * x2 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 - 0.5e1 * x2 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 - 0.5e1 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * x3 * x3 + 4 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 - 4 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * x3 + 4 * x1 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * y2 - 0.5e1 * x1 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * x2 - 0.10e2 * x3 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x2 - 0.15e2 * x3 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 * y2 + 0.20e2 * x3 * x3 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y2 - 2 * y2 * y2 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 + 0.10e2 * x2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * y2 - 0.5e1 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 * y3 - 0.5e1 * x3 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 - 0.5e1 * x3 * x3 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * x3 - 0.5e1 * pow(y1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 - 0.5e1 * pow(y3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 - 2 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(y2, 0.3e1) - x2 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 + 4 * x1 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * y2 - 2 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 - 2 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x3 + 8 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * y1 * y3 * x1 - 0.10e2 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x2 * y2 - 4 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * y2 - 2 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y2 * y2 - 2 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y2 - 0.5e1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 * y3 - 2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * y3 * y1 * y2 + 4 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * y2 * y1 * x2 + x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * x2 + 0.15e2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * x1 * y2 - 0.5e1 * pow(y1, 0.3e1) * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 - 0.5e1 * y1 * y1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 0.5e1 * y1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(y2, 0.3e1) - 2 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 * x3 + 2 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 * x2 + 4 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x1 - 0.5e1 * x2 * x2 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 + 0.5e1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 * y2 + 4 * y2 * y2 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 - 4 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y3 + 4 * x1 * x1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y1 + 4 * y1 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 + 4 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y1 - 0.10e2 * x2 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 - 0.5e1 * y1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(x2, 0.3e1) - 4 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y3 + 4 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * y2 - 4 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * y1 - 4 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y2 + 4 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 + 0.5e1 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 * y3 - 0.5e1 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y3 * y3 - 0.5e1 * pow(y1, 0.3e1) * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 0.5e1 * y1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 + 0.10e2 * y1 * y1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 - 0.5e1 * y1 * y1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 - 0.5e1 * y1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 + 0.10e2 * y1 * y1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 + 4 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 + 4 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x3 + 4 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x3 - 4 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y2 - 2 * x1 * y2 * y2 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 + 4 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 - 2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y3 * y2 * y2 - 0.15e2 * y1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 * y2 - 0.10e2 * y1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * x2 - 4 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x3 * y1 - 0.15e2 * y1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 * x2 + 0.20e2 * y1 * y1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * x2 + 0.5e1 * y2 * y2 * pow(y3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) + 2 * y2 * pow(y3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) + 4 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 0.5e1 * y1 * y1 * pow(y2, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) + 0.5e1 * y2 * y2 * pow(x3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) - 2 * y2 * y2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) - 0.5e1 * pow(y2, 0.3e1) * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) + 0.5e1 * pow(y1, 0.3e1) * y2 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) + 0.5e1 * x2 * pow(y1, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - 2 * x2 * y1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 - 2 * pow(y3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 + 2 * pow(y3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 + 2 * pow(y1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y3 - 2 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 - 2 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - 2 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * x1 - 4 * x1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y2 - 2 * pow(y1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 + 4 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y3 + exp(-0.5e1 * x3 - 0.5e1 * y3) * x3 * pow(y2, 0.3e1) + x3 * pow(y1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) - 4 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * x3 + 2 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 + exp(-0.5e1 * x1 - 0.5e1 * y1) * x1 * pow(y2, 0.3e1) + pow(x1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y3 - pow(x3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 2 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 - x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(y2, 0.3e1) - 0.10e2 * pow(y3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * x2 + 0.10e2 * pow(y3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 + 0.15e2 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 * y2 - 0.10e2 * pow(y3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * x1 + 0.10e2 * pow(y3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * x2 - 2 * y1 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x3 - 2 * y1 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y3 * y2 + 0.10e2 * pow(y3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * x1 - 0.3e1 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y3 * y2 * y2 - 4 * x1 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y1 - pow(x1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 + 6 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 * y2 - x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(y2, 0.3e1) + 0.5e1 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 * y3 + 4 * x1 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y3 + x1 * pow(y3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) - 0.5e1 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 * y2 + pow(x3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 - 4 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x1 * y3 - 2 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 - 2 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 - 2 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 - 2 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - x2 * pow(y3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) - 2 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * x3 - 0.5e1 * pow(y1, 0.3e1) * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 + 4 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 * y3 * y3 - pow(x2, 0.3e1) * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) - x2 * pow(y1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) - 0.5e1 * x1 * x1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 - 0.5e1 * x1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(y2, 0.3e1) + 0.10e2 * x1 * x1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 + 0.15e2 * x1 * x1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 - 0.5e1 * x1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(x2, 0.3e1) + 0.15e2 * x1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 - pow(x2, 0.3e1) * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) + pow(x2, 0.3e1) * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) + 2 * x2 * x2 * y1 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) + 2 * x1 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - 0.3e1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y2 * y2 + 6 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y2 * y2 + 2 * x2 * x2 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) - 4 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 + 4 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y3 + pow(x2, 0.3e1) * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) + 0.5e1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y3 * y3 - 0.5e1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 * y3 - 2 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 + 0.5e1 * pow(x1, 0.3e1) * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) - 0.5e1 * pow(y1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 0.5e1 * pow(y1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * x3 - 0.5e1 * pow(y1, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 + 0.5e1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(x2, 0.3e1) + 0.5e1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(x2, 0.3e1) - 0.5e1 * pow(y3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 - 0.5e1 * pow(y3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - 0.5e1 * pow(y3, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * x1 + 0.5e1 * pow(x1, 0.3e1) * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y3 - 2 * x2 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x3 - 2 * x2 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * x3 + 0.15e2 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y3 - 2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 - 2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - 0.5e1 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y1 + 0.15e2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y3 * y3 - 0.5e1 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 * y1 + 0.20e2 * x1 * x1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y2 - 0.5e1 * x1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 - 0.5e1 * pow(x1, 0.3e1) * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 - 0.5e1 * pow(x1, 0.3e1) * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 0.5e1 * x1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 + 0.10e2 * x1 * x1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 + 0.5e1 * x2 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - 2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 * x2 * x3 - 2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y3 * y1 * y2 - 0.5e1 * pow(x1, 0.3e1) * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 - 0.5e1 * x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y1 - 0.10e2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y3 * x2 - x2 * x2 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 + 4 * x2 * x2 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - x2 * x2 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 0.5e1 * x1 * x1 * y2 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 - x2 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 + 2 * x2 * x2 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 - 0.10e2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y1 * x2 - 0.10e2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y1 * y2 - 0.3e1 * x2 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 - x2 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * x1 - 0.5e1 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y2 * y2 + x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 - 0.15e2 * x1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 * x2 + 0.40e2 * x1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y1 * x2 - 0.10e2 * x1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * x2 - 0.15e2 * x1 * x1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * x2 - 2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 * y3 * y2 + 0.10e2 * x1 * pow(y2, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 - 2 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y3 * y2 * y2 - 4 * x1 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y2 - 2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 * y2 * y2 + 0.15e2 * y1 * y2 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * x1 - 2 * x2 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 2 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y2 - 2 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y1 - 0.5e1 * pow(y1, 0.3e1) * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 - 0.5e1 * x1 * x1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 0.10e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * x3 + x1 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 2 * x2 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 + 0.5e1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(x2, 0.3e1) - 0.10e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 + 2 * x1 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 - 0.10e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 0.10e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 - 0.15e2 * x1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 * x2 - 0.10e2 * x1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 - 0.15e2 * x1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 * y2 - 0.15e2 * x1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 * y2 - 0.5e1 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 0.5e1 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 - 0.5e1 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * x3 - 0.5e1 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 + 4 * y2 * y2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 + 0.5e1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(y2, 0.3e1) - 0.10e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * y1 + 0.15e2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 * x2 - 0.10e2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * y1 * x2 + 0.5e1 * x1 * pow(y3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 + 0.20e2 * x1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 * y1 - 0.10e2 * x1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * x2 - 0.15e2 * x1 * x1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y1 + 0.20e2 * x1 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 * y2 + 0.10e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y1 + 0.20e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x1 + 8 * x2 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 * y1 + 0.15e2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 * y2 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 * x1 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x1 * x1 - 0.20e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * x1 + 4 * x1 * y2 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 - 0.5e1 * y2 * y2 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 + 4 * y2 * y2 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 - 0.10e2 * pow(y2, 0.3e1) * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 - 0.5e1 * y2 * y2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 + 0.10e2 * y2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * x1 - 0.10e2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y1 * x2 - 0.10e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 + 0.10e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * y2 - 0.20e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * x2 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 * x1 + 0.20e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 + 0.10e2 * y2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 - 0.5e1 * y2 * y2 * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 + 0.10e2 * pow(y2, 0.3e1) * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 - 0.20e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * x3 + 0.10e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y3 - 0.20e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x2 - 0.10e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * x3 + 0.10e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * y2 + 0.15e2 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 * y2 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * x1 * y2 + 0.20e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y2 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 * y2 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 * x2 - 2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y3 * x2 * x1 - 0.5e1 * pow(x1, 0.3e1) * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 + 0.10e2 * x1 * x1 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * y3 - 0.10e2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x1 * y2 - 0.5e1 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 + 0.10e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y1 - 0.10e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y2 - 0.10e2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y1 * x1 - 0.5e1 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x2 + 0.20e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * x2 - 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 * y2 + 0.20e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * x1 - 0.5e1 * x1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y2 + 0.10e2 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x1 - 0.10e2 * x1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y1 * x2 - 0.10e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - 0.5e1 * x1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * x2 + 0.20e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * x2 - 0.5e1 * x1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * y2 * y2 - 0.5e1 * x1 * x1 * y2 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 + 0.20e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x3 - 0.5e1 * x1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 * y1 + 0.15e2 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 * x2 - 0.10e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y2 - 2 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 - 0.5e1 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * x1 + 0.10e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 + 0.10e2 * x1 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x3 + 0.20e2 * x1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x3 + 0.10e2 * pow(y2, 0.3e1) * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 + 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(x2, 0.3e1) - 0.10e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * x1 + 0.40e2 * y2 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x1 * y1 - 0.5e1 * y3 * x3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 + 0.10e2 * y2 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 - x2 * x1 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 + 0.20e2 * y2 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * x1 + 2 * x2 * x1 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 + 2 * x2 * x2 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 + 8 * x2 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 - x2 * x2 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 0.10e2 * y2 * y2 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x1 - 0.10e2 * y2 * y2 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * x1 - 0.5e1 * pow(y2, 0.3e1) * x3 * x3 * exp(-0.5e1 * x2 - 0.5e1 * y2) - 0.10e2 * y3 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 - 0.10e2 * y2 * y2 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * x3 + 0.5e1 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(y2, 0.3e1) + 0.20e2 * y2 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 * x3 - 0.15e2 * y1 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * x1 * x1 - 0.15e2 * y1 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * x1 * x3 - 2 * y1 * pow(y3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) - 0.10e2 * y1 * y2 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x1 - 0.5e1 * y1 * pow(y3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 - 0.5e1 * y1 * y1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 0.5e1 * x3 * x3 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 - 0.5e1 * x3 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(y2, 0.3e1) + 0.20e2 * y2 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * x1 * x3 - 0.10e2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x1 * y2 - 2 * x2 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 - 0.5e1 * y2 * pow(x3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 + 2 * pow(y1, 0.3e1) * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) - 0.5e1 * x1 * x1 * pow(y2, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) - 2 * x2 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 + x2 * x2 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 + x1 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 2 * x2 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 - 0.10e2 * y1 * y1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x3 + 0.40e2 * y1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x3 * y2 - 0.5e1 * y1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(x2, 0.3e1) - 0.15e2 * y1 * y1 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 * x1 - 0.5e1 * y1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(y2, 0.3e1) + 0.5e1 * x1 * pow(x3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 - 0.10e2 * x2 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 * x1 + 0.15e2 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 * y3 + 0.15e2 * x1 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y3 - 0.5e1 * y1 * pow(y3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 + 0.15e2 * y1 * y1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 + 0.10e2 * y1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 - 0.10e2 * x3 * y1 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y2 * y3 * x2 - 0.15e2 * y1 * y1 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 * x1 - 0.10e2 * x1 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x2 * x3 * y2 + 4 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y2 - 0.15e2 * y1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 * x2 - 0.15e2 * y1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x3 * x3 - 2 * x1 * x1 * y2 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) + 4 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * x1 + 2 * x2 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 2 * pow(y1, 0.3e1) * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) - 2 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(y2, 0.3e1) - 4 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y2 + 4 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 * y1 - 0.5e1 * y2 * pow(x3, 0.3e1) * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 + 0.5e1 * pow(x1, 0.3e1) * y2 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) - 4 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * y1 + 4 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * y2 + x2 * x2 * x1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 + 0.20e2 * y1 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * x3 * x3 + 0.15e2 * y1 * y1 * y2 * y2 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 - 0.5e1 * y1 * y1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 + 2 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 + x1 * x1 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 - 2 * x2 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 + 2 * x2 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 + x2 * x2 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 - 0.5e1 * y2 * y2 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 + x2 * x2 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 + 4 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y3 + 4 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x3 - 2 * x1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y3 * y1 * y2 + 4 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y3 - 0.5e1 * y2 * y2 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * y1 + 6 * x2 * y1 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 + 2 * x2 * x3 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 0.5e1 * x2 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - x3 * y1 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 0.15e2 * y2 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * x3 * x3 - 2 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * x3 - 0.15e2 * y2 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * x3 * x3 - 4 * x2 * x2 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 - 4 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * x3 - 4 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * y2 + 4 * x1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x3 - 2 * x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y1 * y2 - x3 * y1 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 * y2 - x1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 * x1 - 2 * y3 * x3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 * y1 - x2 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * x3 + 4 * x2 * x2 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - x1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 0.3e1 * x2 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 * y3 - 0.15e2 * y1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 * x2 - 0.10e2 * y1 * y1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x2 + 0.20e2 * y1 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * x2 - 0.15e2 * y1 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 * x2 * y2 - 0.15e2 * y2 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 * x3 + 0.5e1 * x2 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - x2 * x2 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * x3 + 6 * x2 * y3 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 0.5e1 * x2 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 * y2 - 0.15e2 * y2 * y3 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * y1 * x3 + 0.5e1 * pow(y1, 0.3e1) * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) + 2 * exp(-0.5e1 * x3 - 0.5e1 * y3) * y3 * pow(y2, 0.3e1) + 2 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y1 * pow(y2, 0.3e1) + 0.5e1 * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) * pow(y2, 0.3e1) + 0.5e1 * y3 * y3 * exp(-0.5e1 * x3 - 0.5e1 * y3) * pow(y2, 0.3e1) + 0.5e1 * y1 * y1 * pow(y3, 0.3e1) * exp(-0.5e1 * x1 - 0.5e1 * y1) + 0.5e1 * pow(x3, 0.3e1) * y1 * y1 * exp(-0.5e1 * x1 - 0.5e1 * y1) + 0.5e1 * x3 * pow(y1, 0.3e1) * exp(-0.5e1 * x3 - 0.5e1 * y3) * y3 - 4 * x2 * x2 * y3 * exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - 4 * y2 * y2 * y3 * exp(-0.5e1 * x2 - 0.5e1 * y2) * x3) / (x1 - x3 + y1 - y3) / (x1 - x2 + y1 - y2) / (x2 - x3 + y2 - y3) / (exp(-0.5e1 * x1 - 0.5e1 * y1) * x2 - x3 * exp(-0.5e1 * x1 - 0.5e1 * y1) + exp(-0.5e1 * x1 - 0.5e1 * y1) * y2 - exp(-0.5e1 * x1 - 0.5e1 * y1) * y3 + exp(-0.5e1 * x3 - 0.5e1 * y3) * x1 - exp(-0.5e1 * x3 - 0.5e1 * y3) * x2 + exp(-0.5e1 * x3 - 0.5e1 * y3) * y1 - exp(-0.5e1 * x3 - 0.5e1 * y3) * y2 + exp(-0.5e1 * x2 - 0.5e1 * y2) * x3 + exp(-0.5e1 * x2 - 0.5e1 * y2) * y3 - exp(-0.5e1 * x2 - 0.5e1 * y2) * x1 - exp(-0.5e1 * x2 - 0.5e1 * y2) * y1) / 0.5e1;
}

void CPCC2D::Perturbation(double h)
{	
    double xx,yy;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {	
        CCVTVertex	*v = *vit;
        /*xx = v->pos.x() + h*(rand()*2.0/RAND_MAX-1.0);
        yy = v->pos.y() + h*(rand()*2.0/RAND_MAX-1.0);*/

        xx = v->pos.x() + h*(rand()*2.0/RAND_MAX-1.0)*v->avglen;
        yy = v->pos.y() + h*(rand()*2.0/RAND_MAX-1.0)*v->avglen;

        if (isInBound(xx,yy)){
            Point newPos = Point(xx,yy);
            m_points[v->id] = WeightedPoint(newPos, m_points[v->id].weight());
        }
    }
}

void CPCC2D::ClearWeights()
{	
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {	
        CCVTVertex	*v = *vit;
        m_points[v->id] = WeightedPoint(v->pos, INIT_K);
    }
}

double CPCC2D::doLloydIteration()
{	
    double	totalerror = 0.0;
    double pcent = m_pcent;
    
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {	
        CCVTVertex	*v = *vit;

        if (v->type != VALID)
        {
            continue;
        }
        
        if (m_bCirclePacking)
        {
            if (v->isRestricted)
            {
                m_points[v->id] = WeightedPoint(v->pos, m_points[v->id].weight());
            }
			else{
                //m_points[v->id] = WeightedPoint(v->centroid, m_points[v->id].weight());
                Point newpos = v->pos + pcent*(v->centroid - v->pos);
                m_points[v->id] = WeightedPoint(newpos, m_points[v->id].weight());
            }			
        }else
        {
            if (v->isRestricted)
            {
                m_points[v->id] = WeightedPoint(v->pos, v->weight);
            }
			else{
                //m_points[v->id] = WeightedPoint(v->centroid, v->weight);
                Point newpos = v->pos + pcent*(v->centroid - v->pos);
                m_points[v->id] = WeightedPoint(newpos, v->weight);
            }
        }

        double error = CGAL::squared_distance(v->pos, v->centroid);
        totalerror += error;
    }

    Update();	
    m_nIter ++;
    if (m_nIter >= MAX_ITER){
        m_nIter = 0;
    }

    totalerror = EvaluateErrorDist();
    return totalerror;
}

//////////////////////////////////////////////////////////////////////////
////////////////////////        Lloyd Error    //////////////////////////
//////////////////////////////////////////////////////////////////////////

void CPCC2D::FindMaxInscribedCircles()
{
    Point pt; double new_r=-1;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex *v = *vit;
        if (v->type != VALID){
            continue;
        }

        if ((m_boundType == RECTBOUND) ||(m_boundType == POLYBOUND))
        {
            new_r = FindMaxInscribedCircle_Polygon(v, pt);
        }else if (m_boundType == CIRCLEBOUND)
        {
            new_r = FindMaxInscribedCircle_Mix(v, pt);
        }
        if (new_r < 0)
        {
            continue;
        }
        v->inCenter = pt;
        v->inRadius = new_r;
    }
}

void CPCC2D::FindGradient(double& minr)
{
    Point pt; 
    int minid = -1;
    double mindis = -1;
    minr = 10e10;
    Vector2D minNor;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex *v = *vit;
        if (v->type != VALID){
            continue;
        }		
        if ((m_boundType == RECTBOUND) ||(m_boundType == POLYBOUND))
        {
            minid = FindMinDistanceEdge_Polygon(v, pt, minNor,mindis);
            if (mindis<minr){
                minr = mindis;
            }
        }else if (m_boundType == CIRCLEBOUND)
        {
            //new_r = FindMaxInscribedCircle_Mix(v, pt);
        }
        if (minid < 0)
        {
            continue;
        }
        v->gradient = minNor;		
    }
}

void CPCC2D::FindGradient_IC(double& minr)
{
    Point pt;	
    double new_r;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex *v = *vit;
        if (v->type != VALID){
            continue;
        }		
        if ((m_boundType == RECTBOUND) ||(m_boundType == POLYBOUND))
        {
            new_r = FindMaxInscribedCircle_Polygon(v, pt);
        }else if (m_boundType == CIRCLEBOUND)
        {
            new_r = FindMaxInscribedCircle_Mix(v, pt);
        }		
        v->gradient = Vector2D(-pt.x() + v->pos.x(), -pt.y() + v->pos.y());		
    }
}

void CPCC2D::funcgrad_CVT(const double* x, double& f, double* g) {
    set_vertices(x);
    int cur_i = 0 ;
    double energy = 0;
    double trieng;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){
        CCVTVertex *_v = *vit;
        int n = _v->GetVorVertexNumber();
        double area = 0;
        VorVertex * P1, *P2;
        for(int i=0; i<n; i++){
            P1 = _v->GetVorVertex(i);
            P2 = _v->GetVorVertex((i+1)%n);
            Vector2D V1(P1->x() - _v->x(), P1->y() - _v->y());
            Vector2D V2(P2->x() - _v->x(), P2->y() - _v->y());
            double area_qua = fabs(V1^V2);
            trieng = area_qua* (V1*V1 + V1*V2 + V2*V2)/ 12;
            area += area_qua;

            energy += trieng;
        } 
        Vector2D grad = area*Vector2D(_v->pos.x() - _v->centroid.x(), _v->pos.y() - _v->centroid.y());
        _v->gradient = grad;

        if (_v->isRestricted){
            _v->gradient = Vector2D(0,0);
        }
        g[cur_i]   = _v->gradient.x ;
        g[cur_i+1] = _v->gradient.y ;
        cur_i += 2 ;
    }
    f = energy;
}

void CPCC2D::funcgrad_CP(const double* x, double& f, double* g) {	
    set_vertices(x) ;
    f = 0.0 ;
    int cur_i = 0 ;
    int minid = -1;
    Vector2D minNor;
    double minr = -1;	
    Point pt;
    double new_r, init_r, mink = 10e10;	

    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){
        CCVTVertex *v = *vit;
        if ((v->type != VALID) ||(!isInBound(v->pos.x(), v->pos.y()))){
            f += 30;
            return;
        }
        if ((m_boundType == RECTBOUND) ||(m_boundType == POLYBOUND))
        {
            minid = FindMinDistanceEdge_Polygon(v, pt, minNor, new_r);
            
        }else if (m_boundType == CIRCLEBOUND){
            //new_r = FindMaxInscribedCircle_Mix(v, pt);
            //todo???
        }
        if (minid < 0){
            continue;
        }
        v->gradient = minNor;		

        g[cur_i]   = v->gradient.x ;
        g[cur_i+1] = v->gradient.y ;
        cur_i += 2 ;				
        
        init_r = sqrt(m_points[v->id].weight());
        double k = new_r/init_r;
        if (k<mink){
            mink = k;
        }
    }	
    f = 1.0/mink;
    m_k = mink;
}

void CPCC2D::funcgrad_IC(const double* x, double& f, double* g) {	
    set_vertices(x) ;
    f = 0.0 ;
    int cur_i = 0 ;	
    Point pt;
    double new_r, init_r, mink = 10e10;

    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){
        CCVTVertex *v = *vit;
        if ((v->type != VALID)||(!isInBound(v->pos.x(), v->pos.y()))){
            f += 30;
            return;
        }
        if ((m_boundType == RECTBOUND) ||(m_boundType == POLYBOUND))
        {
            new_r = FindMaxInscribedCircle_Polygon(v, pt);
        }else if (m_boundType == CIRCLEBOUND)
        {
            new_r = FindMaxInscribedCircle_Mix(v, pt);
        }
        v->gradient = Vector2D(-pt.x() + v->pos.x(), -pt.y() + v->pos.y());

        g[cur_i]   = v->gradient.x ;
        g[cur_i+1] = v->gradient.y ;
        cur_i += 2 ;
        
        init_r = sqrt(m_points[v->id].weight());
        double k = new_r/init_r;
        if (k<mink){
            mink = k;
        }
    }	
    f = 1.0/mink;
    m_k = mink;
}

double CPCC2D::CP_error(){
    Point pt;
    double f = 0;
    double new_r, init_r, mink = 10e10;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){
        CCVTVertex *v = *vit;
        if ((v->type != VALID)||(!isInBound(v->pos.x(), v->pos.y()))){
            f += 30;
            return f;
        }
        if ((m_boundType == RECTBOUND) ||(m_boundType == POLYBOUND))
        {
            new_r = FindMaxInscribedCircle_Polygon(v, pt);
        }else if (m_boundType == CIRCLEBOUND)
        {
            new_r = FindMaxInscribedCircle_Mix(v, pt);
        }
        init_r = sqrt(m_points[v->id].weight());
        double k = new_r/init_r;
        if (k<mink){
            mink = k;
        }
    }
    f = 1.0/mink;
    return f;
}

void CPCC2D::set_vertices(const double* x)
{
    int id = 0;
    int n = m_points.size()*2;
    double* backx = new double[n];
    get_vertices(backx);
    Point pt;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = *vit;

        if (v->type != VALID){
            continue;
        }

        if (!isInBound(x[id], x[id+1])){
            pt = Point(backx[id], backx[id+1]);
        }else{
            pt = Point(x[id], x[id+1]);
        }
        id += 2;
        m_points[v->id] = WeightedPoint(pt, m_points[v->id].weight());
    }
    Update();
    delete[] backx;
}
void CPCC2D::get_vertices(double* x)
{
    int id = 0;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = *vit;

        if (v->type != VALID){
            continue;
        }

        x[id] = m_points[v->id].x();
        x[id+1] = m_points[v->id].y();
        id += 2;
    }
}
void CPCC2D::FindMinCircumscribedCircles()
{
    Point pt; double new_r=-1;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex *v = *vit;
        if (v->type != VALID)
        {
            continue;
        }

        new_r = FindMinCircumscribedCircle_Polygon(v, pt);
        if (new_r < 0)
        {
            continue;
        }
        v->outCenter = pt;
        v->outRadius = new_r;
    }
}


void CPCC2D::DoLloydAll(int& iternum)
{	
    double ferr = 10e10;
    int totalnum = 0;
    
    while(ferr > 1e-14)
    {
        ferr = doLloydIteration();
        totalnum++;
    }

    iternum = totalnum;
}

double CPCC2D::CirclePackingAll(int& iternum, bool withStatic)//withStatic: true default
{	
    double ferr, ferr1;
    int num1, res = 0, totalnum = 0;

    /*CString strFileR;
    strFileR.Format("c:\\debug_err.txt");
    std::ofstream fpr(strFileR, std::ios::out);*/
        
    double delta = 10e10;

    TNT::Stopwatch mwatchCP;
    mwatchCP.start();
    while(delta > 1e-12 && totalnum < 1000)
    {
        double old_k = m_k;
        if (withStatic)
        {
            ferr1 = 10e10;
            num1 = 0;
            while((ferr1 > 1e-8) && (num1 < 20))// 1.0e-14)
            {
                ferr1 = doCirclePackingIteration_NochangeAlpha();
                num1++;
                //fpr << std::scientific << ferr1 << "\n";
            }

			///////////////////////////////////////////////////////////
			//step show
			/*if(0 == (totalnum%5))
			{
				Show();
			}*/
			///////////////////////////////////////////////////////////

        }
        //fpr << "IterNum: " << iternum << "\n";
        res = doCirclePackingIteration(ferr, CPK);
        //fpr << std::setprecision(16) << std::scientific << ferr << "\n" << m_k << "\n";

        delta = abs(m_k-old_k);
        totalnum++;
    }
    mwatchCP.stop();
	//jerrry
	//Show();

    iternum = totalnum;
    //fpr.close();
    return mwatchCP.read();
}

//jerrry
double CPCC2D::CirclePackingAll2(int& iternum, bool withStatic)//withStatic: true default
{	
    double ferr, ferr1;
    int num1, res = 0, totalnum = 0;

    /*CString strFileR;
    strFileR.Format("c:\\debug_err.txt");
    std::ofstream fpr(strFileR, std::ios::out);*/
        
    double delta = 10e10;

    TNT::Stopwatch mwatchCP;
    mwatchCP.start();
    while(delta > 1e-12 && totalnum < 1000)
    {
        double old_k = m_k;
        if (withStatic)
        {
            ferr1 = 10e10;
            num1 = 0;
            while((ferr1 > 1e-8) && (num1 < 20))// 1.0e-14)
            {
                ferr1 = doCirclePackingIteration_NochangeAlpha2();
                num1++;
                //fpr << std::scientific << ferr1 << "\n";
            }

			///////////////////////////////////////////////////////////
			//step show
			/*if(0 == (totalnum%5))
			{
				Show();
			}*/
			///////////////////////////////////////////////////////////

        }
        //fpr << "IterNum: " << iternum << "\n";
        res = doCirclePackingIteration2(ferr, CPK);
        //fpr << std::setprecision(16) << std::scientific << ferr << "\n" << m_k << "\n";

        delta = abs(m_k-old_k);
        totalnum++;
    }
    mwatchCP.stop();
	//jerrry
	//Show();

    iternum = totalnum;
    //fpr.close();
    return mwatchCP.read();
}
//return 0 as no intersection happens;
//return 1 as 2 circles intersect;
//return 2 as some circle intersect with the boundary;
int CPCC2D::doCirclePackingIteration(double& ferr, const double k)
{	
    double delta = 0;

    if (m_bCirclePackingUseLloyd)
    {
        ferr = doLloydIteration();
    }else
    {
		//jerrry
        //ferr = doCirclePackingInscribedCenter(k, delta);
		ferr = doCirclePackingInscribedCenter(k, delta);
    }
        
    m_k += delta;

    Update();
    m_nIter ++;
    if (m_nIter >= MAX_ITER)
	{
        m_nIter = 0;
    }

    return 0;
}

//jerrry
int CPCC2D::doCirclePackingIteration2(double& ferr, const double k)
{	
    double delta = 0;

    if (m_bCirclePackingUseLloyd)
    {
        ferr = doLloydIteration();
    }else
    {
		//jerrry
        //ferr = doCirclePackingInscribedCenter(k, delta);
		ferr = doCirclePackingInscribedCenter2(k, delta);
    }
        
    m_k += delta;

    Update();
    m_nIter ++;
    if (m_nIter >= MAX_ITER)
	{
        m_nIter = 0;
    }

    return 0;
}
//make sure that this operation won't decrease the compactness
double CPCC2D::doCirclePackingIteration_NochangeAlpha()
{//k won't change
    
    double	k = 0.5;
    double ferr = 1.0;
    double delta = 0;

    if (m_bCirclePackingUseLloyd)
    {
        ferr = doLloydIteration();
    }
	else
    {
        ferr = doCirclePackingInscribedCenter(k, delta);
		//jerrry
		//ferr = doCirclePackingInscribedCenter2(k, delta);
    }

    Update();
    
    return ferr;
}

//jerrry
double CPCC2D::doCirclePackingIteration_NochangeAlpha2()
{//k won't change
    
    double	k = 0.5;
    double ferr = 1.0;
    double delta = 0;

    if (m_bCirclePackingUseLloyd)
    {
        ferr = doLloydIteration();
    }
	else
    {
		//jerrry
        //ferr = doCirclePackingInscribedCenter(k, delta);
		ferr = doCirclePackingInscribedCenter2(k, delta);
    }

    Update();
    
    return ferr;
}

double CPCC2D::doCirclePackingInscribedCenter2(const double k, double& delta)
{
    double	totalerror = 0.0;
    //double newr;
    double initw, mindelta = 10e10;
    double mink = 10e10;	

    bool found = false;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = *vit;
        bool bBounded = false;
        if (v->type != VALID)
		{
            continue;
        }

        Point pt;
        if (v->noweight)
        {
            FindCentroid(v, pt);
        }else
		{
            double new_r;
            if ((m_boundType == RECTBOUND) ||(m_boundType == POLYBOUND))
            {
                new_r = FindMaxInscribedCircle_Polygon(v, pt);
            }else if (m_boundType == CIRCLEBOUND)
            {
                new_r = FindMaxInscribedCircle_Mix(v, pt);
            }
            if (new_r < 0)
            {
                continue;
            }

			///////////////////////////////////////////////////////////
			//Edit by Jerrry
			//bounded condition
			///////////////////////////////////////////////////////////
			std::pair<int, double> pairTapRadii(-1, 0.0);
			bBounded = false;
			pairTapRadii = this->pPairInfo->RadiFresh(v->id, new_r);

			assert(v->id != pairTapRadii.first);
			if(pairTapRadii.second < 0.0)
			{//first bounded bBounded = true
				bBounded = true;
			}

			if(-1 != pairTapRadii.first)
			{//Restrict applyment
				//pair now found position
				// vit->id, new_r, pt
				//cv::Point2d pointTmp1(pt.x(), pt.y());
				Point pointTmp1 = pt;
				//pair old found position
				//pairTapRadii.first, pairTapRadii.second, m_point[pairTapRadii.first].point
				//cv::Point2d pointTmp2(m_points[pairTapRadii.first].x(), 
					//m_points[pairTapRadii.first].y());
				Point pointTmp2(m_points[pairTapRadii.first].x(),
					m_points[pairTapRadii.first].y());
				/*new_r = this->pPairInfo->Radious(v->id, pairTapRadii.first,
					new_r, pairTapRadii.second, pointTmp1, pointTmp2);*/
				
				//pt = Point(pointTmp1.x, pointTmp1.y);
				for(std::list<CCVTVertex*>::iterator vit2 = vertices.begin(); vit2 != vit; ++vit2)
				{//Search the pair
					if((*vit2)->id == pairTapRadii.first)
					{
						double dNew_r = Radious(v->id, pairTapRadii.first, new_r, pairTapRadii.second, 
							this->pPairInfo->RateOf(v->id, pairTapRadii.first), pointTmp1, pointTmp2, **vit2, **vit);
						//more
						new_r = dNew_r;
						//total error refresh
						Point pointTmp = (*vit2)->pos;
						double error = (pointTmp.x()-pointTmp2.x)*(pointTmp.x()-pointTmp2.x)
							+ (pointTmp.y()-pointTmp2.y)*(pointTmp.y()-pointTmp2.y);
						totalerror += error;

						//find mink
						initw = m_points[(*vit2)->id].weight();
						double init_r = sqrt(initw);
						double newk = new_r/init_r;
						if (newk < mink)
						{
							mink = newk;
							found = true;
						}
						//pair first found vertices refresh
						(*vit2)->inCenter = Point(pointTmp2.x, pointTmp2.y);
						(*vit2)->inRadius = new_r;

						//pair first m_points[] refresh
						m_points[(*vit2)->id] = WeightedPoint((*vit2)->inCenter, m_points[v->id].weight());

						break;
					}
				}
			}
			//////////////////////////////////////////////

			//jerrry
			if(!bBounded)
			{//normal mink find
				//this->vertices info refresh
				v->inCenter = pt;
				v->inRadius = new_r;

				WeightedPoint wpTmp = m_points[v->id];
				initw = wpTmp.weight();
				//initw = m_points[v->id].weight();
				//////////
				/*double curr_r = sqrt( v->weight);
				double diff = new_r - curr_r;
				if (diff < mindelta)
				{
					mindelta = diff;
					newr = new_r;
					found = true;
				}*/
				//////////
				double init_r = sqrt(initw);
				double newk = new_r/init_r;
				//jerrry
				if (newk < mink){
					mink = newk;
					found = true;
				}
			}
            //////////
        }
        
		//jerrry
		//m_Points[v->id] info refresh
		m_points[v->id] = WeightedPoint(pt, m_points[v->id].weight());
		if(!bBounded)
		{
        	double error = (v->x() - pt.x()) * (v->x() - pt.x()) + (v->y() - pt.y()) * (v->y() - pt.y());
			totalerror += error;
		}
        
    }
    
    if (!found)
    {
        delta = 0;
    }
	else
	{
        assert(initw!=0);
        //delta = (newr/sqrt(initw) - m_k);
        delta = mink - m_k;
        if (abs(delta) > 1e-10)
        {
            delta = k*delta;
        }
    }
    return totalerror;
}

double CPCC2D::doCirclePackingInscribedCenter(const double k, double& delta)
{
    double	totalerror = 0.0;
    //double newr;
    double initw, mindelta = 10e10;
    double mink = 10e10;	

    bool found = false;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = *vit;
        
        if (v->type != VALID){
            continue;
        }

        Point pt;
        if (v->noweight)
        {
            FindCentroid(v, pt);
        }else{
            double new_r;
            if ((m_boundType == RECTBOUND) ||(m_boundType == POLYBOUND))
            {
                new_r = FindMaxInscribedCircle_Polygon(v, pt);
            }else if (m_boundType == CIRCLEBOUND)
            {
                new_r = FindMaxInscribedCircle_Mix(v, pt);
            }
            if (new_r < 0)
            {
                continue;
            }
            v->inCenter = pt;
            v->inRadius = new_r;
            
            initw = m_points[v->id].weight();
            //////////
            /*double curr_r = sqrt( v->weight);
            double diff = new_r - curr_r;
            if (diff < mindelta)
            {
                mindelta = diff;
                newr = new_r;
                found = true;
            }*/
            //////////
            double init_r = sqrt(initw);
            double newk = new_r/init_r;
            if (newk < mink){
                mink = newk;
                found = true;
            }
            //////////
        }
        
        m_points[v->id] = WeightedPoint(pt, m_points[v->id].weight());
        
        double error = (v->x() - pt.x()) * (v->x() - pt.x()) + (v->y() - pt.y()) * (v->y() - pt.y());
        totalerror += error;
    }
    
    if (!found)
    {
        delta = 0;
    }else{
        assert(initw!=0);
        //delta = (newr/sqrt(initw) - m_k);
        delta = mink - m_k;
        if (abs(delta) > 1e-10)
        {
            delta = k*delta;
        }
    }
    return totalerror;
}

double	CPCC2D::EvaluatePackingRate()
{
    double circles_area = 0;
    double bound_area = 0;
    if (m_boundType == RECTBOUND)
    {
        bound_area = (bound_x_max - bound_x_min)*(bound_y_max - bound_y_min);
    }else if (m_boundType == CIRCLEBOUND)
    {
        bound_area = PI*radius_r*radius_r;
    }else if (m_boundType == POLYBOUND)
    {
        bound_area = poly->PolyArea();
    }

    double minrate = 1e10;
    double maxrate = 0;
    
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = *vit;
        if (v->type != VALID)
        {
            continue;
        }

        circles_area += v->weight;
        
        double rate = (PI*v->weight)/v->area;
        if (rate < minrate)
        {
            minrate = rate;
            vMinCPRate = v;
        }
        if (rate > maxrate)
        {
            maxrate = rate;
            vMaxCPRate = v;			
        }
    }

    circles_area *= PI;
    return (circles_area/bound_area);
}

double	CPCC2D::EvaluatePackingRate2()
{
    double circles_area = 0;
    double bound_area = 0;
    if (m_boundType == RECTBOUND)
    {
        bound_area = (bound_x_max - bound_x_min)*(bound_y_max - bound_y_min);
    }else if (m_boundType == CIRCLEBOUND)
    {
        bound_area = PI*radius_r*radius_r;
    }else if (m_boundType == POLYBOUND)
    {
        bound_area = poly->PolyArea();
    }

    for (unsigned int i=0; i<m_radii.size(); i++)
    {
        circles_area += m_radii[i]*m_radii[i];
    }
    circles_area *= m_k*m_k;
    circles_area *= PI;
    return (circles_area/bound_area);
}

bool CPCC2D::Undo()
{
    if (pLastPD == NULL)	return false;

    m_pdt.clear();
    m_points.clear();
    FreeMemory();
    
    for(std::list<CCVTVertex*>::iterator vit = pLastPD->vertices_begin(); vit != pLastPD->vertices_end(); vit++)
    {
        CCVTVertex	*v = (CCVTVertex*)*vit;
        this->AddPoint(WeightedPoint(v->pos,v->weight));
    }
    
    Update();
    return true;
}

bool CPCC2D::InsideVorCell(CCVTVertex* v, Point pt)
{
    int adjVorVNum = (int)v->adjVorVertex.GetSize();
    assert (adjVorVNum > 2);

    for(int i=0; i<adjVorVNum; i++)
    {
        VorVertex * vv1, * vv2;
        vv1 = (VorVertex *)v->adjVorVertex.GetAt(i);
        vv2 = (VorVertex *)v->adjVorVertex.GetAt( (i+1)%v->adjVorVertex.GetSize());
        if (!IsTriangleCCW(pt, vv1->pos, vv2->pos))
        {
            return false;
        }
    }
    return true;
}

//return the radius
//return -1 when error occurs
double CPCC2D::FindMaxInscribedCircle_Mix(CCVTVertex* v, Point& pt)
{
    const int adjVorVNum = (int)v->adjVorVertex.GetSize();
    if (adjVorVNum < 2)
    {
        return -1;
    }
    std::vector<Plane> planes; // # of planes should equal to adjVorVNum;
    std::vector<Point3D> interPts; //# of intersection pts should be less than C(adjVorVNum,3)
    Point p1, p2;
        
    double tmp = M_SQRT2_2;//0.5*sqrtf(2);
    int arcNum = 0;
    for(int i=0; i<adjVorVNum; i++)
    {
        VorVertex * vv1, * vv2;
        vv1 = (VorVertex *)v->adjVorVertex.GetAt(i);
        vv2 = (VorVertex *)v->adjVorVertex.GetAt( (i+1)%adjVorVNum);

        if ((vv1->type == INTERSECT) && (vv2->type == INTERSECT))
        {
            p1 = vv1->pos;
            p2 = vv2->pos;
            arcNum++;
            assert(v->isAdjBorder);
        }else{
            Vector3D vec(vv2->pos.x()-vv1->pos.x(), vv2->pos.y()-vv1->pos.y(),0);
            Vector3D vecN(-vec.y, vec.x, 0);
            vecN.SetUnit();
            vecN= Vector3D(tmp*vecN.x, tmp*vecN.y, tmp);
            //vec.SetUnit();
            //Vector3D nor = vec^vecN;
            Vector3D nor(-vecN.x, -vecN.y, vecN.z);
            //RotatePointAroundVector(vecN, vec, nor);
            Vector3D pt(vv1->pos.x(), vv1->pos.y(), 0);
            double d =  nor*pt;
            Plane pl;
            pl.d = d;
            pl.normal = nor;
            planes.push_back(pl);
        }
    }

    if (arcNum>2)
    {//deal with the special case that a cell has two arcs.
        double res = FindMaxInscribedCircle_Polygon(v, pt);
        return res;
    }

    int planeNum = (int)planes.size();
    if (planeNum == 0)
    {
        assert(adjVorVNum == 2);
    
        VorVertex * vv1, * vv2;
        vv1 = (VorVertex *)v->adjVorVertex.GetAt(0);
        vv2 = (VorVertex *)v->adjVorVertex.GetAt(1);
        //compute the circle tangent to the container circle and current line segment
        double cr, cx, cy;
        int sign = -1;
        double dist = disPt2Line(vv1->pos.x(), vv1->pos.y(), vv2->pos.x(), vv2->pos.y(), center_x, center_y, cx, cy);
        
        Vector2D vecK(v->pos.x()-cx, v->pos.y()-cy);
        vecK.SetUnit();

        double tmpd = (cx-center_x)*(cx-center_x) + (cy-center_y)*(cy-center_y);

        if (abs(tmpd) > TOLER)
        {
            Vector2D vecD(center_x-cx, center_y-cy);
            vecD.SetUnit();
            double tmp = vecD*vecK;
            if (tmp > 0)
            {// the center of the bound circle and the vertex are on the same side of bisector
                cr = (radius_r + dist)*0.5;
                sign = 1;
            }else{
                cr = (radius_r - dist)*0.5;
            }

            if (sign == 1)
            {
                double t = cr/sqrt(tmpd);
                cx = cx + t*(center_x - cx);
                cy = cy + t*(center_y - cy);
            }else if(sign == -1){
                double t = 1+cr/sqrt(tmpd);
                cx = center_x + t*(cx-center_x);
                cy = center_y + t*(cy-center_y);
            }
        }else{
            cr = radius_r * 0.5;
            cx = center_x + vecK.x*cr;
            cy = center_y + vecK.y*cr;
        }

        pt= Point(cx, cy);
        return cr;
    }

    //compute the intersection between the arc and each plane
    double maxcx, maxcy, maxcr;
    maxcx = maxcy = maxcr = -10e10;
    if (v->isAdjBorder)
    {
        for(int i=0; i<adjVorVNum; i++)
        {
            VorVertex * vv1, * vv2;
            vv1 = (VorVertex *)v->adjVorVertex.GetAt(i);
            vv2 = (VorVertex *)v->adjVorVertex.GetAt( (i+1)%adjVorVNum);
            if ((vv1->type != INTERSECT) || (vv2->type != INTERSECT))
            {
                //compute the circle tangent to the container circle and current line segment
                double cr, cx, cy;
                int sign = -1;
                double dist = disPt2Line(vv1->pos.x(), vv1->pos.y(), vv2->pos.x(), vv2->pos.y(), center_x, center_y, cx, cy);
                
                Vector2D vecK(v->pos.x()-cx, v->pos.y()-cy);
                vecK.SetUnit();
                
                double tmpd = (cx-center_x)*(cx-center_x) + (cy-center_y)*(cy-center_y);

                if (abs(tmpd) > TOLER)
                {
                    Vector2D vecD(center_x-cx, center_y-cy);
                    vecD.SetUnit();
                    double tmp = vecD*vecK;
                    if (tmp > 0)
                    {
                        cr = (radius_r + dist)*0.5;
                        sign = 1;
                    }else{
                        cr = (radius_r - dist)*0.5;
                    }

                    if (sign==1)
                    {
                        double t = cr/sqrt(tmpd);
                        cx = cx + t*(center_x - cx);
                        cy = cy + t*(center_y - cy);
                    }else if(sign==-1){
                        double t = 1+cr/sqrt(tmpd);
                        cx = center_x + t*(cx-center_x);
                        cy = center_y + t*(cy-center_y);
                    }
                }else{
                    cr = radius_r * 0.5;
                    cx = center_x + vecK.x*cr;
                    cy = center_y + vecK.y*cr;
                }
                //check whether this circle intersects other line segments
                bool bvalid = true;
                for(int j=0; j<adjVorVNum; j++)
                {
                    VorVertex * tmpv1, * tmpv2;
                    double tmpx, tmpy;
                    if (!bvalid)
                    {
                        break;
                    }

                    if (j!=i)
                    {
                        tmpv1 = (VorVertex *)v->adjVorVertex.GetAt(j);
                        tmpv2 = (VorVertex *)v->adjVorVertex.GetAt( (j+1)%adjVorVNum);
                        if  ((tmpv1->type != INTERSECT) || (tmpv2->type != INTERSECT))
                        {
                            double d = signedDisPt2Line(tmpv1->pos.x(), tmpv1->pos.y(), tmpv2->pos.x(), tmpv2->pos.y(), cx, cy, tmpx, tmpy);
                            if (cr > d)
                            {
                                bvalid = false;
                            }
                        }
                    }
                }
                //end of check validity
                if (bvalid)
                {
                    if (cr > maxcr)
                    {
                        maxcr = cr;
                        maxcx = cx;
                        maxcy = cy;
                    }
                }
            }
        }
    }

    if (maxcr > 0)
    {
        pt= Point(maxcx, maxcy);
        return maxcr;
    }

    //compute the intersection of three planes
    for (int i =0; i<planeNum-2; i++)
    {
        for (int j = i+1; j<planeNum-1; j++)
        {
            for (int k =j+1; k<planeNum; k++)
            {
                Point3D interpt;
                if (ThreePlanesInter(planes[i], planes[j], planes[k], interpt))
                {
                    if(InsideVorCell(v, Point(interpt.x, interpt.y)))//if (isInCircle(interpt.x, interpt.y))
                    {
                        interPts.push_back(interpt);
                    }
                }
            }
        }
    }	

    if (v->isAdjBorder)
    {//compute the intersection of two planes and the arc
        for(int i =0; i<planeNum-1; i++)
        {
            for (int j=i+1; j<planeNum; j++)
            {
                Point3D pt1, pt2;
                bool b2res= false;
                if (TwoPlanesConeInter(planes[i], planes[j], p1, p2, b2res, pt1, pt2))
                {
                    if (b2res)
                    {
                        interPts.push_back(pt1);
                        interPts.push_back(pt2);
                    }else{
                        interPts.push_back(pt1);
                    }
                }
            }
        }
    }	

    std::vector<Point3D> interPts2; //# should be less than that of interPts
    for (int i =0; i<(int)interPts.size(); i++)
    {
        Vector3D interpt(interPts[i].x, interPts[i].y, interPts[i].z);
        bool isvalid = true;
        if ((interPts[i].z < radius_r)&&(isInsideCone(interPts[i])))
        {			
            for (int j =0; j<planeNum; j++)
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
    }

    assert((int)interPts2.size()>0);
    
    double maxz = 0;
    Point3D maxzPt;
    for (int i =0; i<(int)interPts2.size(); i++)
    {
        if (interPts2[i].z > maxz)
        {
            maxzPt = interPts2[i];
            maxz = interPts2[i].z;
        }
    }

    pt= Point(maxzPt.x, maxzPt.y);
    
    double mind = 10e10;
    for(int i=0; i<adjVorVNum; i++)
    {
        VorVertex * vv1, * vv2;
        vv1 = (VorVertex *)v->adjVorVertex.GetAt(i);
        vv2 = (VorVertex *)v->adjVorVertex.GetAt( (i+1)%adjVorVNum);
        
        double dist;
        if ((vv1->type == INTERSECT) && (vv2->type == INTERSECT)){
            Point pcenter = Point(center_x, center_y);
            dist = radius_r - sqrt(CGAL::squared_distance(pt, pcenter));
        }else{
            double tmpx, tmpy;
            dist = disPt2Line(vv1->pos.x(), vv1->pos.y(), vv2->pos.x(), vv2->pos.y(), pt.x(), pt.y(), tmpx, tmpy);
        }
        
        if (dist < mind)
        {
            mind = dist;
        }
    }

    //assert(InsideVorCell(v, pt));

    ////////
    //if(!InsideVorCell(v, pt)){
    //	std::ofstream	fp ("c:\\error.txt",std::ios::out);
    //	std::vector<WeightedPoint>& points = m_points;

    //	fp << 2 << ' ' << (int)points.size() << std::endl; 
    //	for(std::vector<WeightedPoint>::iterator pit = points.begin(); pit != points.end(); pit++)
    //	{
    //		WeightedPoint wpt = *pit;
    //		Point pt = wpt.point();
    //		double r = wpt.weight();
    //		fp << std::scientific << (double)pt.x() << ' ' << (double)pt.y()<<' ' << r <<"\n";
    //	}
    //	fp.close();
    //}
    //////////
    return mind;
}

bool CPCC2D::ExportBoundary(std::string	strFile)
{
    std::ofstream	fp(strFile,std::ios::out);	
    if (!fp.good())
    {
        return false;
    }

    fp << (int)ptsBoundary.size()<< std::endl; 
    for(std::vector<Point>::iterator pit = ptsBoundary.begin(); pit != ptsBoundary.end(); pit++)
    {
        Point pt = *pit;
        fp << std::setprecision(16) << std::scientific << (double)pt.x() << ' ' << (double)pt.y()<<"\n";
    }
    fp.close();
    return true;
}

bool CPCC2D::ExportCPResult2(std::string	strFile, ios_base::openmode mode)
{
    std::ofstream	fp(strFile, mode);
    bool withRestrict = false;

    if (!fp.good())
    {
        return false;
    }

    fp << (int)m_points.size()<< ' ' << std::setprecision(16) << m_k  << std::endl;
    /*if (m_boundType == POLYBOUND)
    {
        fp << (int)ptsBoundary.size()<< std::endl; 
        for(std::vector<Point>::iterator pit = ptsBoundary.begin(); pit != ptsBoundary.end(); pit++)
        {
            Point pt = *pit;
            fp << std::setprecision(16) << std::scientific << (double)pt.x() << ' ' << (double)pt.y()<<"\n";
        }
    }*/
    assert(!withRestrict);
    for(std::vector<WeightedPoint>::iterator pit = m_points.begin(); pit != m_points.end(); pit++)
    {
        WeightedPoint wpt = *pit;
        Point pt = wpt.point();
        double r = sqrt(wpt.weight());
        fp << std::setprecision(16) << std::scientific << (double)pt.x() << ' ' << (double)pt.y()<<' ' << r <<"\n";
    }	

    fp.close();
    return true;
}


bool CPCC2D::ExportCPResult(std::string	strFile, ios_base::openmode mode)
{
    std::ofstream	fp(strFile, mode);
    bool withRestrict = false;

    if (!fp.good())
    {
        return false;
    }

    fp << 2 << ' ' << (int)m_points.size()<< ' ' << m_boundType  << std::endl;
    if (m_boundType == POLYBOUND)
    {
        fp << (int)ptsBoundary.size()<< std::endl; 
        for(std::vector<Point>::iterator pit = ptsBoundary.begin(); pit != ptsBoundary.end(); pit++)
        {
            Point pt = *pit;
            fp << std::setprecision(16) << std::scientific << (double)pt.x() << ' ' << (double)pt.y()<<"\n";
        }
    }
    if (withRestrict)
    {
        for(std::list<CCVTVertex*>::iterator vit=vertices.begin(); vit != vertices.end(); vit++)
        {
            CCVTVertex * v = *vit;
            double w =  v->weight;
            Point pt = v->pos;
            if (v->isRestricted)
            {
                fp << std::setprecision(16) << std::scientific << (double)pt.x() << ' ' << (double)pt.y()<<' ' << w <<' ' << 1<< "\n";
            }else{
                fp << std::setprecision(16) << std::scientific << (double)pt.x() << ' ' << (double)pt.y()<<' ' << w <<' ' << 0<< "\n";
            }
        }		
    }else{
        for(std::vector<WeightedPoint>::iterator pit = m_points.begin(); pit != m_points.end(); pit++)
        {
            WeightedPoint wpt = *pit;
            Point pt = wpt.point();
            double w = m_k*m_k*wpt.weight();
            fp << std::setprecision(16) << std::scientific << (double)pt.x() << ' ' << (double)pt.y()<<' ' << w <<"\n";
        }
    }
    
    fp.close();
    return true;
}

bool CPCC2D::ExportCPRadii(std::string	strFile)
{
    std::ofstream	fp(strFile,std::ios::out);	
    if (!fp.good())
    {
        return false;
    }

    fp << 2 << ' ' << (int)m_points.size()<< ' ' << m_boundType  << std::endl; 
    for(std::vector<WeightedPoint>::iterator pit = m_points.begin(); pit != m_points.end(); pit++)
    {
        WeightedPoint wpt = *pit;		
        double r = m_k*sqrt(wpt.weight());
        fp << std::setprecision(16) << std::scientific <<  r <<"\n";
    }
    /*for (int i=0; i<m_points.size(); i++){
        fp << std::setprecision(16) << std::scientific <<  1.0/sqrt((double)(i+1)) <<"\n";
    }*/
    fp.close();
    return true;
}

bool CPCC2D::LoadBoundary(std::string strFile)
{
    std::ifstream	fp(strFile,std::ios::in);
    if (!fp.good())
    {
        return false;
    }
    int size;
    double x, y;
    fp >> size;	
        
    SetBoundaryType(POLYBOUND);
    ptsBoundary.clear();

    for(int i=0; i<size; i++){
        fp >>x >> y;
        ptsBoundary.push_back(Point(x, y));
    }
    fp.close();

    if (poly != NULL) delete poly;

    poly = new CPolygon();	
    poly->Init(ptsBoundary);
    
    return true;
}

bool CPCC2D::assignPositions(std::string strFile)
{
    std::ifstream	fp(strFile,std::ios::in);
    if (!fp.good())
    {
        return false;
    }
    int size;
    double x, y, k, w;
    //fp >> tag >> size >>type; //dimension; # of the circles  0 for Rect container, 1 for Circle container, 2 for polygon boundary
    fp >> size >> k; 

    for(int i=0; i<size; i++){
        fp >>x >> y >> w;		
        m_points[i] = WeightedPoint(Point(x,y), m_points[i].weight());
    }
    m_k = k;
    Update();
    return true;
}

bool CPCC2D::ImportCPPts(std::string strFile, BoundType& bound, bool withRestrict)
{
    std::ifstream	fp(strFile,std::ios::in);
    if (!fp.good())
    {
        return false;
    }
    int tag, size, type, polysize;
    double x, y, w;
    fp >> tag >> size; //dimension; # of the circles
    fp >> type; // 0 for Rect container, 1 for Circle container, 2 for polygon boundary
    
    bound = (BoundType)type;
    SetBoundaryType(bound);
    m_radii.clear();

    if (bound == POLYBOUND)
    {
        ptsBoundary.clear();
        fp >> polysize;	
        for(int i=0; i<polysize; i++){
            fp >>x >> y;
            ptsBoundary.push_back(Point(x, y));
        }

        if (poly != NULL) delete poly;
        poly = new CPolygon();	
        poly->Init(ptsBoundary);
    }

    if (tag == -1) //no weight
    {
        for(int i=0; i<size; i++){
            fp >>x >> y;
            AddPoint(WeightedPoint(Point(x, y), INIT_K));
            m_radii.push_back(sqrt(INIT_K));
        }
    }else{
        if (withRestrict)
        {
            int restrict;
            for(int i=0; i<size; i++){
                fp >> x >> y >> w >> restrict;
                if (restrict == 1)
                {
                    AddRestrictedPoint(WeightedPoint(Point(x, y), w));

                }else if (restrict == 0)
                {
                    AddPoint(WeightedPoint(Point(x, y), w));
                }
                m_radii.push_back(sqrt(w));
            }
        }else{
            for(int i=0; i<size; i++){
                fp >>x >> y >> w;
                AddPoint(WeightedPoint(Point(x, y), w));
                m_radii.push_back(sqrt(w));
            }
        }
    }
    
    fp.close();
    m_cpSize = (int)m_radii.size();
    m_bSetRadii = true;

    if (bound != 1)
    {
        m_initMethod = RANDOM;
    }
    
    return true;
}

bool CPCC2D::ImportCVTPts(std::string strFile, bool withRestrict)
{
    std::ifstream	fp(strFile,std::ios::in);
    if (!fp.good())	{
        return false;
    }
    int tag, size;
    double x, y, w = 1e-6;
    fp >> tag >> size; //dimension; # of the circles
    
    SetBoundaryType(RECTBOUND);
    m_radii.clear();
    
    if (withRestrict)
    {
        int restrict;
        for(int i=0; i<size; i++){
            fp >> x >> y >> restrict;
            if (restrict == 1)
            {
                AddRestrictedPoint(WeightedPoint(Point(x, y), w));

            }else if (restrict == 0)
            {
                AddPoint(WeightedPoint(Point(x, y), w));
            }
            m_radii.push_back(sqrt(w));
        }
    }else{
        for(int i=0; i<size; i++){
            fp >>x >> y;
            AddPoint(WeightedPoint(Point(x, y), w));
            m_radii.push_back(sqrt(w));
        }
    }
    
    fp.close();
    m_cpSize = (int)m_radii.size();
    m_bSetRadii = true;
    m_initMethod = RANDOM;
    
    return true;
}

bool CPCC2D::ImportCPRadii(std::string	strFile, BoundType& bound)
{
    std::ifstream	fp(strFile,std::ios::in);
    if (!fp.good())
    {
        return false;
    }
    int size, type;
    double r;
    fp >> size >> size;
    fp >> type; //0 for Rect container, 1 for Circle container, 2 for polygon boundary
    bound = (BoundType)type;
    SetBoundaryType(bound);
    m_radii.clear();
    for(int i=0; i<size; i++){
        fp >> r;
        m_radii.push_back(r);
    }
    fp.close();
    m_cpSize = (int)m_radii.size();
    m_bSetRadii = true;

    if (bound != CIRCLEBOUND)
    {
        m_initMethod = RANDOM;
    }
    CirclePackingInitialization();
    
    return true;
}

void CPCC2D::AssignACP1Radii(const int num){
    m_radii.clear();
    for(int i=0; i<num; i++)
	{
        double r = 1.0/sqrt((double)(i+1));
        m_radii.push_back(r);
    }
    
    m_cpSize = (int)m_radii.size();
    m_bSetRadii = true;
}

void CPCC2D::CirclePackingInitialization()
{
    m_bCirclePacking = true; 

    switch (m_initMethod)
    {
    case DS1:
        GenerateSeedsInitDS1(m_cpSize);	
        break;
    case DS2:
        GenerateSeedsInitDS2(m_cpSize);
        break;
    case RANDOM:
        GenerateRandomSeeds(m_cpSize);
        break;
    case A10:
        GenerateSeedsInitA10New(m_cpSize);
        break;
    case A15:
        GenerateSeedsInitA15(m_cpSize);
        break;
    case CC:
        GenerateSeedsInitCC(m_cpSize);
        break;
    }
    Update();
}

//to test whether each circle is active or not.
//void CPCC2D::checkCirclesActiveProperty()
//{
//	Edge_circulator	myec;
//
//	for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
//	{
//		std::vector<CGAL::Object> mytangents;
//		CCVTVertex	*v = (CCVTVertex*)*vit;
//		Point2D tangentPt, tmppt;
//		Vector2D tangentDir;
//		if (v->isAdjBorder)
//		{
//			if(If1CircleTouchesBoundary(v))
//			{
//				bool btangent = intersec2Circles(v->pos.x(), v->pos.y(), sqrt(v->weight), center_x, center_y, radius_r, tangentPt, tmppt);
//				assert(btangent);
//				tangentDir = Vector2D(-(v->pos.y()-center_y), (v->pos.x()-center_x));
//				Ray r = Ray(Point(tangentPt.x, tangentPt.y), Vector(tangentDir.x, tangentDir.y));				
//				mytangents.push_back(make_object(r));
//			}
//		}		
//		
//		for(int j=0; j<v->adjVertice.GetSize(); j++) 
//		{
//			CCVTVertex	*vj = (CCVTVertex *)v->adjVertice.GetAt(j);
//			if(If2CirclesTouchEachOther(v, vj))
//			{
//				bool bfound = FindMyEdge(v,vj,myec);
//				assert(bfound);
//				CGAL::Object o = m_pdt.dual(myec);
//				mytangents.push_back(o);
//			}
//		}
//		//judge whether the tangent points form a closed area
//		int num =(int)mytangents.size();
//		if (num <3)
//		{
//			v->isActive = true;
//			continue;
//		}
//		Point* polypts= new Point[num];
//		bool bparallel = false;
//		for(int j=0; j<num; j++){
//			Point pt;
//			if (!FindIntrsectionPoint(mytangents[j], mytangents[(j+1)%num], pt))
//			{
//				v->isActive = true;
//				bparallel= true;
//				break;
//			}else{
//				polypts[j] = pt;
//			}
//		}
//		if (!bparallel)
//		{
//			if (isInsidePolygon(polypts, num, v->pos))
//			{
//				v->isActive = false;
//			} 
//			else{
//				v->isActive = true;
//			}
//			
//		}
//		delete [] polypts;
//	}	
//}

bool CPCC2D::FindIntrsectionPoint(CGAL::Object o1, CGAL::Object o2, Point& pt)
{
    Segment	s1, s2;
    Ray		r1, r2;
    double x1,y1,x2,y2,x3,y3,x4,y4;
    if (CGAL::assign(s1,o1)) {
        x1 = s1.source().x();
        y1 = s1.source().y();
        x2 = s1.target().x();
        y2 = s1.target().y();		
    }
    else if (CGAL::assign(r1,o1)){
        x1 = r1.source().x();
        y1 = r1.source().y();
        x2 = r1.source().x()+r1.direction().dx();
        y2 = r1.source().y()+r1.direction().dy();
    }else{ return false;}

    if (CGAL::assign(s2,o2)){
        x3 = s2.source().x();
        y3 = s2.source().y();
        x4 = s2.target().x();
        y4 = s2.target().y();
    }
    else if (CGAL::assign(r2,o2)){
        x3 = r2.source().x();
        y3 = r2.source().y();
        x4 = r2.source().x()+r2.direction().dx();
        y4 = r2.source().y()+r2.direction().dy();
    }else{ return false;}

    double denom = (y4-y3)*(x2-x1)-(x4-x3)*(y2-y1);
    if(abs(denom)<TOLER) return false;
    double ua = ((x4-x3)*(y1-y3)-(y4-y3)*(x1-x3))/denom;
    //double ub = ((x2-x1)*(y1-y3)-(y2-y1)*(x1-x3))/denom;

    pt= Point(x1+ua*(x2-x1), y1+ua*(y2-y1));
    return true;
}

bool CPCC2D::FindMyEdge(CCVTVertex *v1, CCVTVertex *v2, Edge_circulator& myec)
{
    Edge_circulator	ec = m_pdt.incident_edges(v1->vh), done = ec;
    Segment edge;
    do {
        if( m_pdt.is_infinite(ec) ){
            ec ++;
            continue;
        }
        edge = m_pdt.segment(ec);

        if (IsSamePoint(edge.target(), v2->pos) || IsSamePoint(edge.source(), v2->pos))
        {
            myec = ec;
            return true;
        }
        ec ++;
    } while (ec != done);
    return false;
}

//select a minimal m_k
bool CPCC2D::GenerateSeedsInitA15(const int num)
{
    if((num <1)||(!m_bSetRadii) ||(!m_bCirclePacking) ||(m_boundType != CIRCLEBOUND && m_boundType != RECTBOUND))
    {
        return false;
    }

    assert(m_cpSize == num);

    //init vars
    bool succ = false;
    double x0, y0;	
    double resx[2]={0};
    double resy[2]={0};	
    double r_alpha = m_k;

    double* backupRadii = new double[num];
    for (int k =0; k<num; k++) backupRadii[k] = m_radii[k];

    double* currRadii = new double[num];
    Vector2D* currPos = new Vector2D[num];

    std::vector<CornerPlacement> L;

    for (int k = 0; k<num-1; k++)
    {
        for (int l = k+1; l<num; l++)
        {			
            //generate an initial config using circles k and l
            double todork = m_radii[k]*r_alpha;
            double todorl = m_radii[l]*r_alpha;
            
            bool res;
            if (m_boundType==CIRCLEBOUND){
                x0 = -(radius_r - todork);
                y0 = 0;
                res = tangentto2Circles_inscribe(center_x, center_y, radius_r, x0, y0, todork, todorl, resx, resy);
            }else if (m_boundType==RECTBOUND)
            {//put to the lower left corner of the rect
                x0 = bound_x_min + todork;
                y0 = bound_y_min + todork;
                res = tangentto2Circles_inscribe_rect(bound_x_min, bound_y_min, bound_x_max, bound_y_max,x0, y0, todork, todorl, resx, resy);
            }
            currPos[0] = Vector2D(x0, y0);
            currRadii[0] = m_radii[k];

            if (!res)
            {
                continue;
            }else{
                currPos[1] = Vector2D(resx[0], resy[0]);
                currRadii[1] = m_radii[l];
                m_radii.erase(m_radii.begin()+l);
                m_radii.erase(m_radii.begin()+k);
            }
            //***end of init with c_k and c_l***

            //***A1.5Core***			
            res = A15Core(2, currPos, currRadii, num, m_radii, L);
            if (res){
                succ = true;
                goto stop;
            }

            //restore m_radii
            m_radii.clear();
            for (int k =0; k<num; k++)
                m_radii.push_back(backupRadii[k]);
        }
    }

stop:
    if (succ)
    {
        //restore m_radii
        m_radii.clear();
        for (int k =0; k<num; k++)
            m_radii.push_back(backupRadii[k]);

        m_pdt.clear();
        m_points.clear();
        FreeMemory();

        for (int k=0; k<num; k++)
        {
            AddPoint(WeightedPoint(Point(currPos[k].x, currPos[k].y), currRadii[k]*currRadii[k]));
        }
    }
    delete[] backupRadii;
    delete[] currRadii;
    delete[] currPos;
    return succ;
}

bool CPCC2D::GenerateSeedsInitA10New(const int num){
    if((num <1)||(!m_bSetRadii) ||(!m_bCirclePacking) ||(m_boundType != CIRCLEBOUND && m_boundType != RECTBOUND))
    {
        return false;
    }

    assert(m_cpSize == num);

    //init vars
    bool succ = false;
    double x0, y0;	
    double resx[2]={0};
    double resy[2]={0};	
    double r_alpha = m_k;
    
    double* backupRadii = new double[num];
    for (int k =0; k<num; k++) backupRadii[k] = m_radii[k];
    
    double* currRadii = new double[num];
    Vector2D* currPos = new Vector2D[num];

    std::vector<CornerPlacement> L;

    for (int k = 0; k<num-1; k++)
    {
        for (int l = k+1; l<num; l++)
        {			
            //generate an initial config using circles k and l
            double todork = m_radii[k]*r_alpha;
            double todorl = m_radii[l]*r_alpha;

            bool res;
            if (m_boundType==CIRCLEBOUND){
                x0 = -(radius_r - todork);
                y0 = 0;
                res = tangentto2Circles_inscribe(center_x, center_y, radius_r, x0, y0, todork, todorl, resx, resy);
            }else if (m_boundType==RECTBOUND)
            {//put to the lower left corner of the rect
                x0 = bound_x_min + todork;
                y0 = bound_y_min + todork;
                res = tangentto2Circles_inscribe_rect(bound_x_min, bound_y_min, bound_x_max, bound_y_max,x0, y0, todork, todorl, resx, resy);
            }
            currPos[0] = Vector2D(x0, y0);
            currRadii[0] = m_radii[k];

            if (!res)
            {
                continue;
            }else{
                currPos[1] = Vector2D(resx[0], resy[0]);
                currRadii[1] = m_radii[l];
                m_radii.erase(m_radii.begin()+l);
                m_radii.erase(m_radii.begin()+k);
            }
            //***end of init with c_k and c_l***

            //***A1.0Core***
            double density;
            res = A10Core(2, currPos, currRadii, num, m_radii, L, density);
            if (res){
                succ = true;
                goto stop;
            }

            //restore m_radii
            m_radii.clear();
            for (int k =0; k<num; k++)
                m_radii.push_back(backupRadii[k]);
        }
    }

stop:
    if (succ)
    {
        //restore m_radii
        m_radii.clear();
        for (int k =0; k<num; k++)
            m_radii.push_back(backupRadii[k]);
        m_pdt.clear();
        m_points.clear();
        FreeMemory();
        
        for (int k=0; k<num; k++)
        {
            AddPoint(WeightedPoint(Point(currPos[k].x, currPos[k].y), currRadii[k]*currRadii[k]));
        }
    }
    delete[] backupRadii;
    delete[] currRadii;
    delete[] currPos;
    return succ;
}

//add the perturbation of m_k
bool CPCC2D::GenerateSeedsInitA10(const int num)
{
//to test whether we can find better configuration using current container(center and radius)
//method based on [Wang06], using maximal hole degree rule, A1.0	
    if((num <1)||(!m_bSetRadii) ||(!m_bCirclePacking) ||(m_boundType != CIRCLEBOUND))
    {
        return false;
    }

    assert(m_cpSize == num);
    Vector2D* initpts = new Vector2D[m_cpSize];

    //init vars
    bool succ = false;
    double x0, y0;	
    double resx[2]={0};
    double resy[2]={0};	

    double* backupRadii = new double[num];
    for (int k =0; k<num; k++) backupRadii[k] = m_radii[k];

    double* currRadii = new double[num];
    int testNum = 0;
    
    //not use this code for the moment, just suppose that the first biggest radii are the init configuration
    for (int k = 0; k<num-1; k++)
    {
        for (int l = k+1; l<num; l++)
        {			
            //***init with c_k and c_l***
            double todork = m_radii[k]*m_k;
            double todorl = m_radii[l]*m_k;
            x0 = -(radius_r - todork);
            y0 = 0;
            initpts[0] = Vector2D(x0, y0);
            currRadii[0] = m_radii[k];
            
            bool res = tangentto2Circles_inscribe(center_x, center_y, radius_r, x0, y0, todork, todorl, resx, resy);
            if (!res)
            {
                continue;
            }else{
                initpts[1] = Vector2D(resx[0], resy[0]);
                currRadii[1] = m_radii[l];
                m_radii.erase(m_radii.begin()+l);
                m_radii.erase(m_radii.begin()+k);
            }
            //***end of init with c_k and c_l***

            //***A1.0Core***
            int i = 2;
            for (; i<num; i++)
            {		
                double x, y, tmpx, tmpy, lamda;
                double maxlamda = -10e10;
                int tmpi=-1;
                for (int j = 0; j< num-i; j++)
                {
                    double r = m_radii[j];

                    res = placeNewCircleIntoCurrentConfig(i, radius_r, center_x, center_y, initpts, currRadii, r, x, y, lamda);
                    if (!res) continue;
                    if(lamda > maxlamda){
                        maxlamda = lamda;
                        tmpx = x;
                        tmpy = y;
                        tmpi = j;
                    }
                }
                //make sure that we do have found some circle
                if (tmpi == -1)
                {
                    break;
                }
                initpts[i] = Vector2D(tmpx, tmpy);
                currRadii[i] = m_radii[tmpi];
                m_radii.erase(m_radii.begin()+tmpi);
            }
            //***A1.0Core***
            testNum++;

            if (i == num)
            {
                succ = true;
                break;
            }

            //restore m_radii
            m_radii.clear();
            for (int k =0; k<num; k++)
                m_radii.push_back(backupRadii[k]);
        }
        if (succ)
        {
            break;
        }
    }	

    //restore m_radii
    m_radii.clear();
    for (int k =0; k<num; k++)
        m_radii.push_back(backupRadii[k]);

    if (succ)
    {
        m_pdt.clear();
        m_points.clear();
        FreeMemory();
        
        for (int k=0; k<num; k++)
        {
            AddPoint(WeightedPoint(Point(initpts[k].x, initpts[k].y), currRadii[k]*currRadii[k]));
        }		
    }

    delete[] currRadii;
    delete[] backupRadii;
    delete[] initpts;

    return succ;
}

bool CPCC2D::A15Core(const int currnum, Vector2D* currPos, double* currRadii,const int totalnum, std::vector<double>& leftRadii,
                     std::vector<CornerPlacement>& L)
{
    bool succ = false;
    Vector2D* COPYcurrPos = new Vector2D[totalnum];
    double* COPYcurrRadii = new double[totalnum];
    std::vector<double> COPYleftRadii;

    int updatednum = currnum;
    generateCornerPlacementList(currnum, currPos, currRadii, totalnum, leftRadii, L);
    while (L.size()>0)
    {
        double maxbenefit = -1e10;
        int maxid = -1;		
        for (unsigned int i=0; i<L.size(); i++)
        {
            //copy C and L
            int newnum = updatednum;
            for(int k=0; k<totalnum; k++){
                COPYcurrRadii[k] = currRadii[k];
                COPYcurrPos[k] = currPos[k];
            }
            COPYleftRadii.clear();
            for (int k=0; k<totalnum-updatednum; k++)
            {
                COPYleftRadii.push_back(leftRadii[k]);
            }
            std::vector<CornerPlacement> COPYL;

            //modify C' by placing ci at (x,y)
            COPYcurrPos[newnum] = Vector2D(L[i].x, L[i].y);
            COPYcurrRadii[newnum] = L[i].r;
            newnum++;
            //modify L'
            std::vector<double>::iterator iter = find(COPYleftRadii.begin(), COPYleftRadii.end(), L[i].r);
            if (iter != COPYleftRadii.end()) {
                COPYleftRadii.erase(iter);
            }else{
                assert(false);
            }
            double density;
            bool res = A10Core(newnum, COPYcurrPos, COPYcurrRadii, totalnum, COPYleftRadii, COPYL, density);
            if (res){				
                //return C'
                for(int k=0; k<totalnum; k++){
                    currRadii[k] = COPYcurrRadii[k];
                    currPos[k] = COPYcurrPos[k];
                }
                succ = true;
                goto stop;				
            }else{
                if (density>maxbenefit){
                    maxbenefit = density;
                    maxid = i;
                }
            }
        }		

        //select the corner placement with max benefit
        //modify C by placing c_{maxid}
        currPos[updatednum] = Vector2D(L[maxid].x, L[maxid].y);
        currRadii[updatednum] = L[maxid].r;
        std::vector<double>::iterator iter = find(leftRadii.begin(), leftRadii.end(), L[maxid].r);
        if (iter != leftRadii.end()) {
            leftRadii.erase(iter);
        }else{
            assert(false);
        }		
        updatednum++;
        //modify L
        generateCornerPlacementList(updatednum, currPos, currRadii, totalnum, leftRadii, L);
    }

    if(updatednum==totalnum){
        succ = true;
    }else{
        succ = false;
    }

stop:
    delete[] COPYcurrPos;
    delete[] COPYcurrRadii;
    return succ;
}

bool CPCC2D::A10Core(const int currnum, Vector2D* currPos, double* currRadii,const int totalnum, std::vector<double>& leftRadii,
                     std::vector<CornerPlacement>& L, double& density)
{
    generateCornerPlacementList(currnum, currPos, currRadii, totalnum, leftRadii, L);

    int newnum = currnum;
    while (L.size()>0)
    {
        double maxbenefit = -1e10;
        double tmpx, tmpy, lamda, tmpr;
        for (unsigned int i=0; i<L.size(); i++)
        {
            lamda = L[i].lamda;
            if(lamda > maxbenefit){
                maxbenefit = lamda;
                tmpx = L[i].x;
                tmpy = L[i].y;
                tmpr = L[i].r;
            }
        }
        //modify C by placing ci at(tmpx, tmpy)
        currPos[newnum] = Vector2D(tmpx, tmpy);
        currRadii[newnum] = tmpr;
        std::vector<double>::iterator iter = find(leftRadii.begin(), leftRadii.end(), tmpr);
        if (iter != leftRadii.end()) {
            leftRadii.erase(iter);
        }else{
            assert(false);
        }		
        newnum++;
        
        //modify L
        generateCornerPlacementList(newnum, currPos, currRadii, totalnum, leftRadii, L);
    }
    
    double insidearea=0;
    for (int k=0; k<newnum; k++){
        insidearea += currRadii[k];
    }
    density = insidearea;

    if(newnum==totalnum){
        return true;
    }else{
        return false;
    }
}

bool CPCC2D::generateCornerPlacementsfor1Circle(const int num, const Vector2D* currPos, const double* currRadii, const double r, 
                                    std::vector<CornerPlacement> &L)
{
    double resx[2]={0};
    double resy[2]={0};
    double tmpx, tmpy, tmplamda;
    double tmpr = r*m_k;
    CornerPlacement tmpcorner;
    
    //***to find the corner placement between the current circles and the container***
    for (int i=0; i<num; i++)
    {
        bool res;
        if (m_boundType == CIRCLEBOUND)
        {
            res	= tangentto2Circles_inscribe(center_x, center_y, radius_r, currPos[i].x, currPos[i].y, currRadii[i]*m_k, tmpr, resx, resy);

        }else if (m_boundType == RECTBOUND)
        {
            res = tangentto2Circles_inscribe_rect(bound_x_min, bound_y_min, bound_x_max, bound_y_max,
                    currPos[i].x, currPos[i].y, currRadii[i]*m_k, tmpr, resx, resy);
        }
        if (!res) continue;
        //***make sure that no overlapping***
        bool succtangent = false;
        for (int m=0; m<2; m++)
        {
            bool overlapped = false;
            for (int k=0; k<num; k++)
            {
                if ((k!=i) && is2CirclesOverlapped(resx[m], resy[m], tmpr, currPos[k].x, currPos[k].y, currRadii[k]*m_k))
                {
                    overlapped = true;
                    break;
                }
            }
            if (overlapped)
            {
                continue;
            }else{
                tmpx = resx[m];
                tmpy = resy[m];
                succtangent = true;
                break;
            }
        }
        //***end of overlapping test***

        if (!succtangent) continue;

        //***compute the hole degree***
        double mindis = 10e10;
        for (int k=0; k<num; k++)
        {
            if (k==i) continue;
            double dist = sqrt((tmpx-currPos[k].x)*(tmpx-currPos[k].x) + (tmpy-currPos[k].y)*(tmpy-currPos[k].y)) - currRadii[k]*m_k - tmpr;
            if (dist < mindis)
            {
                mindis = dist;
            }
        }
        if (m_boundType==RECTBOUND){
            double dist[4];
            dist[0] = abs(bound_y_min)-tmpy-tmpr;
            dist[1] = abs(bound_y_min)+tmpy-tmpr;
            dist[2] = abs(bound_x_min)+tmpx-tmpr;
            dist[3] = abs(bound_x_min)-tmpx-tmpr;
            for (int i=0; i<4; i++)	{
                if (dist[i] < mindis)
                {
                    mindis = dist[i];
                }
            }
        }
        tmplamda = 1 - mindis/tmpr; //lamda should be < 1
        tmpcorner.lamda = tmplamda;
        tmpcorner.x = tmpx;
        tmpcorner.y = tmpy;
        tmpcorner.r = r;
        L.push_back(tmpcorner);
        //***end of hole degree computation***
    }
    //***end of finding the corner placement between the current circles and the container***

    //***to find the corner placement between the current circles***
    for (int k = 0; k<num-1; k++)
    {
        for (int l = k+1; l<num; l++)
        {
            bool res = tangentto2Circles(currPos[k].x, currPos[k].y, currRadii[k]*m_k, currPos[l].x, currPos[l].y, currRadii[l]*m_k, tmpr, resx, resy);
            if (!res) continue;
            bool succtangent = false;
            //***make sure being inside the container***
            bool inside[2] = {false};
            for (int m=0; m<2; m++)
            {
                double sqrtdist = sqrt((resx[m]-center_x)*(resx[m]-center_x) + (resy[m]-center_y)*(resy[m]-center_y));
                if(sqrtdist < (radius_r - tmpr)){
                    inside[m] = true;
                }
            }
            //***end of make sure being inside the container***
            int validnum;
            if (!inside[0] && !inside[1])
            {
                continue;
            }else if (inside[0] && inside[1])
            {
                validnum = 2;
            }else if (inside[0]){
                validnum = 1;
            }else{
                validnum = 1;
                resx[0] = resx[1];
                resy[0] = resy[1];
            }

            //***make sure that no overlapping***
            for (int m=0; m<validnum; m++)
            {
                bool overlapped = false;
                for (int j=0; j<num; j++)
                {
                    if ((j!=k) && (j!=l) && is2CirclesOverlapped(resx[m], resy[m], tmpr, currPos[j].x, currPos[j].y, currRadii[j]*m_k))
                    {
                        overlapped = true;
                        break;
                    }
                }
                if (overlapped)
                {
                    continue;
                }else{
                    tmpx = resx[m];
                    tmpy = resy[m];
                    succtangent = true;
                    break;
                }
            }
            //***end of overlapping test***

            if (!succtangent) continue;

            //***compute the hole degree***
            //including the dist between the circle and the container
            double mindis = sqrt((tmpx-center_x)*(tmpx-center_x) + (tmpy-center_y)*(tmpy-center_y)) - tmpr;
            for (int j=0; j<num; j++)
            {
                if ((j==k)||(j==l)) continue;
                double dist = sqrt((tmpx-currPos[j].x)*(tmpx-currPos[j].x) + (tmpy-currPos[j].y)*(tmpy-currPos[j].y)) - currRadii[j]*m_k - tmpr;
                if (dist < mindis)
                {
                    mindis = dist;
                }
            }
            if (m_boundType==RECTBOUND){
                double dist[4];
                dist[0] = abs(bound_y_min)-tmpy-tmpr;
                dist[1] = abs(bound_y_min)+tmpy-tmpr;
                dist[2] = abs(bound_x_min)+tmpx-tmpr;
                dist[3] = abs(bound_x_min)-tmpx-tmpr;
                for (int i=0; i<4; i++)	{
                    if (dist[i] < mindis)
                    {
                        mindis = dist[i];
                    }
                }
            }
            tmplamda = 1 - mindis/tmpr; //lamda should be < 1
            tmpcorner.lamda = tmplamda;
            tmpcorner.x = tmpx;
            tmpcorner.y = tmpy;
            tmpcorner.r = r;
            L.push_back(tmpcorner);
            //***end of hole degree computation***
        }
    }
    //***end of finding the corner placement between the current circles***

    if (L.size()<1){
        return false;
    }
    return true;
}

void CPCC2D::generateCornerPlacementList(const int currnum, const Vector2D* currPos, const double* currRadii,const int totalnum, const std::vector<double>& leftRadii,	 std::vector<CornerPlacement>& L)
{
    assert(totalnum-currnum == leftRadii.size());
    L.clear();
    for (int k=0; k<totalnum-currnum; k++){
        double r = leftRadii[k];
        generateCornerPlacementsfor1Circle(currnum, currPos, currRadii, r, L);
    }
}

bool coner_comp_greater(CornerPlacement a, CornerPlacement b)
{
    return a.lamda>b.lamda;
};
//output: the placement(x,y) of "r" and its benefit lamda
//note; r is the original radius
bool CPCC2D::placeNewCircleIntoCurrentConfig(const int num, const double radius0, const double centerX, const double centerY, 
                                             const Vector2D* currPts, const double* currRadii, const double r, double& x, double& y, double& lamda)
{		
    std::vector<CornerPlacement> corners;
    bool res = generateCornerPlacementsfor1Circle(num, currPts, currRadii, r, corners);
    if (!res){
        return false;
    }
    std::sort(corners.begin(), corners.end(), coner_comp_greater);
    
    //assign the result
    x = corners[0].x;
    y = corners[0].y;
    lamda = corners[0].lamda;
    return true;
}

void CPCC2D::GenerateSeedsInitCC(const int num)
{
    if((num <1)||(!m_bSetRadii) ||(!m_bCirclePacking) ||(m_boundType != CIRCLEBOUND))
    {
        return;
    }	
    assert(m_cpSize == num);
    if(m_cpSize < 1) return;

    m_pdt.clear();
    m_points.clear();
    FreeMemory();

    Circle2 minc;
    double cr;
    double delta = PI2/num;
    Vector2D* initpts = new Vector2D[num];
    Vector2D* tripts = new Vector2D[3];
    tripts[0] = Vector2D(center_x, center_y);
    for (int i=0; i<num; i++)
    {
        double theta ,x, y;
        for (int k=0; k<2; k++)
        {
            theta = delta*(i+k);
            x = center_x + radius_r*cos(theta);
            y = center_y + radius_r*sin(theta);
            tripts[k+1] = Vector2D(x,y);
        }
        
        MinCircle2* mincircle = new MinCircle2(3, tripts, minc);
        initpts[i] = Vector2D(minc.Center.x,  minc.Center.y);
        cr = minc.Radius;
        delete mincircle;
    }

    double k = cr/m_radii[0];
    m_k = k*1.1;
    double tmp = 1e-2;
    AddPoint(WeightedPoint(Point(initpts[0].x*1.2, initpts[0].y*1.2), m_radii[0]*m_radii[0]));
    for (int i=1; i<num; i++)
    {
        AddPoint(WeightedPoint(Point(initpts[i].x, initpts[i].y), m_radii[i]*m_radii[i]));
    }

    delete[] tripts;
    delete[] initpts;
}

//return the radius
//return -1 when error occurs
double CPCC2D::FindMinCircumscribedCircle_Polygon(CCVTVertex* v, Point& pt)
{
    int adjVorVNum = (int)v->adjVorVertex.GetSize();
    if(adjVorVNum < 2){
        return -1;
    }

    Vector2D* vorpts = new Vector2D[adjVorVNum];
    
    for(int i=0; i<adjVorVNum; i++)
    {
        VorVertex * vv;
        vv = (VorVertex *)v->adjVorVertex.GetAt(i);
        vorpts[i] = Vector2D(vv->pos.x(), vv->pos.y());
    }
    
    Circle2 minc;
    MinCircle2* mincircle = new MinCircle2(adjVorVNum, vorpts, minc);
    
    pt= Point(minc.Center.x, minc.Center.y);
    
    delete[] vorpts;
    delete mincircle;
    return minc.Radius;
}

void CPCC2D::FinishPolyBoundary()
{
    if (poly != NULL) delete poly;
        
    poly = new CPolygon();	
    poly->Init(ptsBoundary);
}

void CPCC2D::addBoundaryPt(Point pt)
{
    ptsBoundary.push_back(pt);
}

bool CPCC2D::PtInPoly(double x, double y)
{
    if (poly == NULL)
    {
        return false;
    }
    return poly->PtInPoly(x,y);
}

int corners[4] = {-1};//(leftbottom, rightbottom, righttop, lefttop)

int nearest_corner(const Vector2D pt){
    double mindis = 10e10;
    int minid = -1;
    Vector2D pa;
    for(int i=0; i<4; i++)	{
        switch (i) {
            case 0:
                pa = Vector2D(bound_x_min, bound_y_min);				
                break;
            case 1:
                pa = Vector2D(bound_x_max, bound_y_min);
                break;
            case 2:
                pa = Vector2D(bound_x_max, bound_y_max);
                break;
            case 3:
                pa = Vector2D(bound_x_min, bound_y_max);
                break;
        }
        double dis = (pt-pa).length2();
        if (dis<mindis){
            mindis = dis;
            minid = i;
        }
    }
    return minid;
}

void CPCC2D::ProjectSeed2Boundary(CCVTVertex* v)
{
    assert (v->isAdjBorder);
    double cx, cy, minx, miny, mindist=1e10;

    if (m_boundType == CIRCLEBOUND)
    {
        double x0 = v->pos.x();
        double y0 = v->pos.y();
        double tmp = 1.0/sqrt(x0*x0+y0*y0);
        minx = x0*tmp;
        miny = y0*tmp;

    }else{
        int adjVorVNum = (int)v->adjVorVertex.GetSize();
        for(int i=0; i<adjVorVNum; i++)
        {
            VorVertex * vv1, * vv2;
            vv1 = (VorVertex *)v->adjVorVertex.GetAt(i);
            vv2 = (VorVertex *)v->adjVorVertex.GetAt((i+1)%adjVorVNum);
            if ((vv1->type == CORNER && vv2->type == INTERSECT)||(vv1->type == INTERSECT && vv2->type == INTERSECT)||(vv1->type == CORNER && vv2->type == CORNER)||(vv1->type == INTERSECT && vv2->type == CORNER))
            {
                if ((m_boundType == RECTBOUND) && (vv1->type == CORNER || vv2->type == CORNER)){
                    int id = nearest_corner(Vector2D( v->pos.x(), v->pos.y()));
                    if (corners[id] <0 ){
                        corners[id]=1;
                        minx = cornorspt[id].x();
                        miny = cornorspt[id].y();
                        break;
                    }
                }

                double dist = disPt2Line(vv1->pos.x(), vv1->pos.y(), vv2->pos.x(), vv2->pos.y(), v->pos.x(), v->pos.y(), cx, cy);
                if (dist < mindist)
                {
                    mindist = dist;
                    minx = cx;
                    miny = cy;
                }
            }
        }
    }
    
    Point newPos = Point(minx, miny);
    m_points[v->id] = WeightedPoint(newPos, m_points[v->id].weight());
}

void CPCC2D::ProjectSeeds2Boundary(bool setrestrict)
{
    //reset corners[4]
    for (int i=0; i<4; i++){
        corners[i]=-1;
    }

    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = *vit;
        if (v->type != VALID){
            continue;
        }
        if (v->isAdjBorder)
        {
            ProjectSeed2Boundary(v);
            if (setrestrict){
                v->isRestricted = true;
            }
        }
    }
    Update();
}

void CPCC2D::RemoveConstriction()
{
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = *vit;
        v->isRestricted = false;
    }	
}

void CPCC2D::LloydProject()
{
    RemoveConstriction();
        
    doLloydIteration();
    ProjectSeeds2Boundary();	
}

// compute the local best coordinates of the vertices
void CPCC2D::compute_new_coords(double *x)
{
    Point pt; 
    double new_r;
    int id = 0;
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++){
        CCVTVertex *v = *vit;
        if ((v->type != VALID)||(!isInBound(v->pos.x(), v->pos.y()))){
            continue;
        }
        if ((m_boundType == RECTBOUND) ||(m_boundType == POLYBOUND))
        {
            new_r = FindMaxInscribedCircle_Polygon(v, pt);
        }else if (m_boundType == CIRCLEBOUND)
        {
            new_r = FindMaxInscribedCircle_Mix(v, pt);
        }
        x[id] = pt.x();
        x[id+1] = pt.y();
        id += 2;
    }
}

void CPCC2D::ClearSeeds()
{
    m_pdt.clear();
    m_points.clear();
    FreeMemory();
}

void CPCC2D::ClearIterNum()
{
    m_nIter = 0;
}

void CPCC2D::SetCPRadii(int rsize, std::vector<double> radii)
{
    assert(rsize>=0);
    
    SetBoundaryType(POLYBOUND);
    m_radii.clear();
    for(int i=0; i<rsize; i++){
        
        m_radii.push_back(radii[i]);
    }
    
    m_cpSize = (int)m_radii.size();
    m_bSetRadii = true;

    m_initMethod = RANDOM;
    
    CirclePackingInitialization();	
}

void CPCC2D::SetBoundary(int bsize, std::vector<Point> boudpts)
{
    assert(bsize>=0);

    SetBoundaryType(POLYBOUND);
    ptsBoundary.clear();

    for(int i=0; i<bsize; i++){		
        ptsBoundary.push_back(boudpts[i]);
    }

    if (poly != NULL) delete poly;
    poly = new CPolygon();	
    poly->Init(ptsBoundary);
    
    Update();
}

bool CPCC2D::ExportCP(std::string strFile, ios_base::openmode mode)
{
    std::ofstream	fp(strFile, mode);
    bool withRestrict = false;

    if (!fp.good())
    {
        return false;
    }

    fp << (int)m_points.size()<< ' ' << std::setprecision(16) << m_k  << std::endl;
    
    if (m_boundType == POLYBOUND)
    {
        fp << (int)ptsBoundary.size()<< std::endl; 
        for(std::vector<Point>::iterator pit = ptsBoundary.begin(); pit != ptsBoundary.end(); pit++)
        {
            Point pt = *pit;
            fp << std::setprecision(16) << std::scientific << (double)pt.x() << ' ' << (double)pt.y()<<"\n";
        }
    }

    assert(!withRestrict);
    for(std::vector<WeightedPoint>::iterator pit = m_points.begin(); pit != m_points.end(); pit++)
    {
        WeightedPoint wpt = *pit;
        Point pt = wpt.point();
        double r = sqrt(wpt.weight());
        fp << std::setprecision(16) << std::scientific << (double)pt.x() << ' ' << (double)pt.y()<<' ' << r <<"\n";
    }	

    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = *vit;

        if (v->type != VALID){
            continue;
        }

        fp<<v->adjVorVertex.GetSize()<<"\n";
        for(int i=0; i<v->adjVorVertex.GetSize(); i++)
        {		
            VorVertex * vv;
            vv = (VorVertex *)v->adjVorVertex.GetAt(i);
            fp << std::setprecision(16) << std::scientific << (double)vv->pos.x() << ' ' << (double)vv->pos.y() <<"\n";
        }
    }

    fp.close();
    return true;
}

bool CPCC2D::ExportCP(std::vector<Point>& aCirclePos, std::vector<double>& aCircleRadius, std::vector<std::vector<Point> >& aRegionDiv)
{
    aCirclePos.clear();
    aCircleRadius.clear();
    aRegionDiv.clear();

    aCirclePos.reserve(this->m_points.size());
    aCircleRadius.reserve(this->m_points.size());
    for(std::vector<WeightedPoint>::iterator pit = m_points.begin(); pit != m_points.end(); pit++)
    {
        WeightedPoint wpt = *pit;
        Point pt = wpt.point();
        double r = sqrt(wpt.weight()) * (this->m_k);
        aCirclePos.push_back(Point(pt.x(), pt.y()));
        aCircleRadius.push_back(r);
    }	

    aRegionDiv.reserve(this->vertices.size());
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = *vit;

        if (v->type != VALID){
            continue;
        }
        size_t nPtNum = v->adjVorVertex.GetSize();
        std::vector<Point> aRegionPt(nPtNum);
        for(size_t i = 0;  i < nPtNum;  ++i)
        {		
            VorVertex * vv;
            vv = (VorVertex *)v->adjVorVertex.GetAt(i);
            aRegionPt[i] = Point(vv->pos.x(), vv->pos.y());
        }
        aRegionDiv.push_back(aRegionPt);
    }

    return true;
}

void CPCC2D::Show()
{
	std::vector<std::vector<cv::Point2d> >  aSubRegion(this->m_points.size());
    std::vector<cv::Point2d>  aCircleCenter(this->m_points.size());
    std::vector<double>  aCircleRadius(this->m_points.size());
	std::vector<std::vector<cv::Point2d>> aPointTrail(this->m_points.size());
	//Circle info
    for(std::vector<WeightedPoint>::iterator pit = m_points.begin(); pit != m_points.end(); pit++)
    {
        WeightedPoint wpt = *pit;
        Point pt = wpt.point();
        double r = sqrt(wpt.weight()) * (this->m_k);
        aCircleCenter.push_back(cv::Point2d(pt.x(), pt.y()));
        aCircleRadius.push_back(r);
    }	

	//SubRegion info
    for(std::list<CCVTVertex*>::iterator vit = vertices.begin(); vit != vertices.end(); vit++)
    {
        CCVTVertex	*v = *vit;

        if (v->type != VALID){
            continue;
        }
        size_t nPtNum = v->adjVorVertex.GetSize();
        std::vector<cv::Point2d> aRegionPt(nPtNum);
        for(size_t i = 0;  i < nPtNum;  ++i)
        {		
            VorVertex * vv;
            vv = (VorVertex *)v->adjVorVertex.GetAt(i);
			aRegionPt[i] = cv::Point2d(vv->pos.x(), vv->pos.y());
        }
        aSubRegion.push_back(aRegionPt);

		size_t nTrailLen = this->m_nIter;//v->footPrints
		std::vector<cv::Point2d> aTrail;
		for(size_t i = 0; i<nTrailLen; ++i)
		{
			aTrail.push_back(cv::Point2d(v->footPrints[i].x(), v->footPrints[i].y()));
		}
		aPointTrail.push_back(aTrail);
    }
	
	//Display
	double fScale = 1.0;
	int nShowW = static_cast<int>(1050.0 * fScale + 0.5);
	int nShowH = static_cast<int>(750.0 * fScale + 0.5);
    cv::Size ShowSize(nShowW, nShowH);
    unsigned int nRegion = aSubRegion.size();

    std::vector<std::vector<cv::Point> >  aContours(nRegion);
    for(unsigned int i = 0;  i < nRegion;  ++i)
    {
        aCircleCenter[i].x *= fScale;
        aCircleCenter[i].y *= fScale;
        aCircleRadius[i]   *= fScale;
        unsigned int nPt = aSubRegion[i].size();
        aContours[i].resize(nPt);
        for(unsigned int j = 0;  j < nPt;  ++j)
        {
            cv::Point2d Pt = aSubRegion[i][j];
            Pt.x *= fScale;
            Pt.y *= fScale;            
            aContours[i][j].x = (int)std::max(0.0, std::min(Pt.x, nShowW*1.0));
            aContours[i][j].y = (int)std::max(0.0, std::min(Pt.y, nShowH*1.0));
        }
    }
    cv::Mat ShowImg = cv::Mat::zeros(ShowSize, CV_8UC3);
	cv::imshow("BlankWindow", ShowImg);
	cv::waitKey();

    cv::drawContours(ShowImg, aContours, -1, cv::Scalar::all(255), 2, CV_AA);
	//circles
    for(unsigned int i = 0;  i < nRegion;  ++i)
    {
        cv::circle(ShowImg, aCircleCenter[i], (int)(aCircleRadius[i] + 0.5), cv::Scalar::all(255), 2, CV_AA);
    }
	//Lines between the pair
	if(NULL != this->pPairInfo)
		this->pPairInfo->Draw(ShowImg, aCircleCenter);

	//Trail
	/*cv::Point2d pointTmp1, pointTmp2;
	for(int i = 0; i<aPointTrail.size(); ++i)
	{
		pointTmp1 = aPointTrail[i][0];
		for(int j = 1; j<aPointTrail[i].size(); ++j)
		{
			pointTmp2 = aPointTrail[i][j];
			cv::line(ShowImg, pointTmp1, pointTmp2, cv::Scalar(255, 0, 0), 2);
			pointTmp1 = pointTmp2;
		}
	}*/
	
    cv::imwrite("StepShow.jpg", ShowImg);
    std::string strWindowName("StepShow");
    cv::namedWindow(strWindowName, CV_WINDOW_AUTOSIZE);
    cv::imshow(strWindowName, ShowImg);
    cv::waitKey(0);
}

