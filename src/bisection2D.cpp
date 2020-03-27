/***********************************************************************/
/**      BISECTION NUMERICAL METHOD TO GET CATENARY CHAIN POINTS       */
/**                                                                    */
/** Open Source Code                                                   */
/** Service Robotics Lab.            http://robotics.upo.es            */ 
/**                                  https://github.com/robotics-upo   */
/**                                                                    */
/** Maintainer:                     Contact:                           */
/** Simon Martinez Rozas            simon.martinez@uantof.cl           */
/**                                                                    */   
/***********************************************************************/
#include <cmath>
#include <iostream>
#include <vector>

//////////////////////////////////////////////////////////////////////////
//Implementation explanation:                                           //
//                                                                      //
//The implemented code get the points position of a CHAIN that fallow   //
//catenary. For that, Bisection Numerical Method is used to obtaining   //
//the value of the parameters that describe this physical phenomenon.   //
//                                                                      //
//In physics and geometry, a catenary is the curve that an idealized    // 
//hanging chain or cable assumes under its own weight when supported    //
//only at its ends.                                                     //
//                                                                      //
//In this case it is assumed that the lashing points (X1,Y1) and        //                                                                     //
//(X2,Y2), which are not in the same height, and the length (L) of the  //
//cable are known.                                                      //
//                                                                      //
//Catenary equation referred to any axes whatever:                      //
//                      y-y0 = c* cosh ((x-x0)/c),                      //
//thus, for this equation to be defined, it is necessary to obtain 3    //
//constants: c, x0 and y0.                                              //
//                                                                      //
//The first equation is the condition of the catenary passing through   //
//point A:                                                              //
//        [i]           -y0 = c* cosh (x0/c)                            //
//                                                                      //
//The second equation is the condition of the catenary passing through  //
//point B:                                                              //
//        [ii]         yB-y0 = c* cosh ((xB-x0)/c)                      //
//                                                                      //
//The third and last equation will be the length of the cable arc.      //
//        [iii]    Sab = c * sin(X0/c)+ c*sin((xB-x0)/c)                //
//                                                                      //
//Then, we subtract, from the second equation, the first equation, and  //
//to this result, we will add the third equation (cable length):        //
//Sab+yB = c*cosh((xB-x0)/c)-c*cosh (x0/c)+c*sin(X0/c)+c*sin((xB-x0)/c) //
//                                                                      //
//Using trigonometry:                                                   //
//  Sab+yB = c*e^((xB-X0)/c) - c*e^(-x0/c) = c*e^(-x0/c)*(e^(xB/c)-1)   //
//                                                                      //
//Similar result is obtain if operation is subtaction:                  //
//  Sab-yB = c*e^((xB-X0)/c) + c*e^(-x0/c) = c*e^(x0/c)*(1+e^(-xB/c))   //
//                                                                      //
//Multiplaying the last two ecuations:                                  //
//               Sab^2-yB^2= c^2*4*[sinh(xB/(2*c))]^2                   //
//                                                                      //
//               sqrt(Sab^2-yB^2)= c*2*sinh(xB/(2*c))                   //
//         sqrt(Sab^2-yB^2) / xB = sinh(xB/(2*c)) / (xB/(2*c))          //
//                                                                      //
//Making the following variable change:   phi = xB/(2*c)                //
//                                          k = sqrt(Sab^2-yB^2)/xB     //
//Then we obtain the first ecuation to resolve by Bisection method:     //
//                        k*phi = sinh(phi)    , then we found c.       //
//                                                                      //
//Finally x0 is found by subtracting [i] and [ii] to eliminate Y0, and  //
//the bisection method is applied to the resulting equation.            //
//The value y0 is found using bisection method in [i].                  //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

using namespace std;

float phi(float x);
float distanceC(float p);
float catenaryPointA(float x);
float catenaryPointB(float y);
float evaluteCatenaryChain(float x_);
void integrateCatenaryChain();
float bisection(float a, float b, int mode_);// 0 = to find phi() , 1 = to find X0 , 2 = to find Y0
//Find points with sign changes in interval a-b, times that function pass through the origin 
void signChange(float a, float b, int mode_);
float modeBisection(float xr_aux, int mode_);

struct points_interval
{
    float pa;
    float pb;
};
vector <points_interval> vector_sign_change;

struct points_catenary
{
    float x_;
    float y_;
};
vector <points_catenary> catenary_chain_points;

double L,X1,Y1,X2,Y2,XB,YB;
float straight_tolerance;
float kConst;

//For Phi Function
double Ap, Bp, Ax, Bx, Ay, By;
float tolerance;
//To save solutions of numeric method in function
double bs_p, bs_Y0, bs_X0;
int n_points, n_chain;
float c_value, h_value, Xc, Yc;

int main() 
{
    printf ("\nSTART FINDING ROOT using BISECTION METHOD, Good Luck !!\n");
    
    //Catenary Chain Values 
    L =50;
    //Lashing points(X1,Y1) and (X2,Y2)
    X1 =10; Y1 = 30;
    X2 = 50; Y2 = 30;
    //Distance between the axes X and Y
    XB= X2 - X1;
    YB= Y2 - Y1;
    //Value to accurete solution
    straight_tolerance = 0.000001; //minimum difference between Sab and L to not consider the chain in a straight state
    tolerance = 0.0001; //
    //Number of chain point to get
    n_chain = 100;
    //Data for use during bisection method
    n_points = 50;
    //Xtreme values Interval to evaluate Phi() Function     
    Ap =0.00001;    Bp = 3;
    //Xtreme values Interval to evaluate Catenary_A() Function     
    Ax =0;          Bx = 100;
    //Xtreme values Interval to evaluate Catenary_B() Function     
    Ay =-100;       By = 100;
    

    //CALCULATE OF CATENARY POINTS
    //Calculate K
    kConst = (sqrt(pow(L,2) - pow(YB,2)))/(XB);
    if ( fabs(kConst) < straight_tolerance)
        printf ("\nTHE CABLE IS STRAIGHT!! NOT IN CATENARY STATE \n");
    printf ("\nValue of k = %f\n",kConst);
    //Look for the solutions c, x0, y0
    //Calculate of phi and c
    bs_p = bisection(Ap, Bp,0);
    c_value = distanceC(bs_p);
    //Calculate of x0
    bs_X0 = bisection(Ax,Bx,1);
    //Calculate of y0
    bs_Y0 = bisection(Ay,By,2);
    h_value = fabs(bs_Y0)-c_value;
    printf("\nValue of 'h = %g'\n",h_value);
    //Calculate lower point in catenary
    Xc=X1+bs_X0;
    Yc=Y1-h_value;
    printf("\nThe lower Point in chain is (Xc,Yc)=(%f,%f)'\n",Xc,Yc);
    printf("\nGET CATENARY POINTS\n");
    //Calculate points for Catenary chain
    integrateCatenaryChain();

    return 0;
}

float phi(float x) 
{
    return (sinh(x)- kConst*x);
}

float distanceC(float p)
{
    float result_ = (XB)/(2*p);
    printf("\nValue of 'c = %g'\n",result_);
    return result_;
}

float catenaryPointA(float x)
{
    return (c_value*cosh(((XB)-x)/c_value) - c_value*cosh((x)/c_value) - YB);
    // return (c_value*cosh(((XB)-x)/c_value) - c_value*cosh((X1-x)/c_value) - YB + Y1);
}

float catenaryPointB(float y)
{
    return (c_value*cosh(((XB)-bs_X0)/c_value) - (YB) + y);
}

float evaluteCatenaryChain(float x_)
{
    return (c_value * cosh((x_- Xc)/c_value)+ (Yc-c_value));
}

void integrateCatenaryChain()
{
    float y_value;
    float x_value = X1;
    float x_step = (X2-X1)/n_chain;
    points_catenary point_cat;
    catenary_chain_points.clear();
    
    for (int i=0 ; i < n_chain +1 ; i ++)
    {
        y_value = evaluteCatenaryChain(x_value);
        point_cat.x_ = x_value;
        point_cat.y_ = y_value;
        printf("Points Catenary x_= %f , y_= %f\n",x_value,y_value);
        catenary_chain_points.push_back(point_cat);
        x_value = x_value + x_step;
    }

    printf("Vector 'catenary_chain_points' size = %i\n",catenary_chain_points.size());
}

float modeBisection(float xr_aux, int mode_)
{
    float yr_aux;
    
    if (mode_ == 0)
    {
        yr_aux = phi(xr_aux);
    }
    if (mode_ == 1)
    {
        yr_aux = catenaryPointA(xr_aux);
    }
    if (mode_ == 2)
    {
        yr_aux = catenaryPointB(xr_aux);  
    }
    
    return yr_aux;
}

float bisection(float a, float b,int mode_) 
{
    float xr, error;
    //Initialize error with a big value
    error = tolerance+1;
    
    float x_a , x_b;
    points_interval it_points_interval;
    
    printf ("\nbisection mode = %i",mode_);

    printf ("\nLooking for Function root\n");
    signChange(a, b, mode_);
    // printf ("\nSize vector %lu\n",vector_sign_change.size());

    for (int i = 0; i < vector_sign_change.size(); i++)
    {
        it_points_interval = vector_sign_change[i];
        x_a = it_points_interval.pa;
        x_b = it_points_interval.pb;
        
        while (error > tolerance) 
        {
            xr = (x_a + x_b) / 2.0;
            if ((modeBisection(x_a,mode_) * modeBisection(xr,mode_) < 0) || 
                (modeBisection(x_a,mode_) * modeBisection(xr,mode_) == 0))
            {
                x_b = xr;
            }
            else if ((modeBisection(xr,mode_) * modeBisection(x_b,mode_)< 0) ||
                (modeBisection(x_b,mode_) * modeBisection(xr,mode_) == 0))
            {
                x_a = xr;
            }
            error = fabs(modeBisection(xr,mode_));
            // printf("OUT of if !! 'PHI(phi) = %f'  'phi = %f'  'x_a = %f'  'x_b = %f'  'error = %f'\n",phi(xr),xr,x_a,x_b,error);
            if (error < tolerance)
            {
                if (mode_ == 0)
                    printf("SOLUTION FOUNDED!!! for 'PHI(phi) = %f'  'phi = %f'\n",modeBisection(xr,mode_),xr);
                if (mode_ == 1)
                    printf("SOLUTION FOUNDED!!! for 'CatenaryA(X0) = %f'  'X0 = %f'\n",modeBisection(xr,mode_),xr);
                if (mode_ == 2)
                    printf("SOLUTION FOUNDED!!! for 'CatenaryB(Y0) = %f'  'Y0 = %f'\n",modeBisection(xr,mode_),xr);
            }
        }
    }


    return xr;
}

void signChange(float a, float b, int mode_)
{
    points_interval points_sign_change;
    vector_sign_change.clear();
    
    float y1 , y2, yt;
    float interval = fabs(b - a) / n_points;
    float bound_a = a;
    float bound_b;
    
    yt = y2*y1;

    printf ("signChange mode = %i\n",mode_);
    // printf ("Entering signChange() for phi()\n");

    for (int i = 0; i < n_points; i++) 
    {
        y1 = modeBisection(bound_a,mode_);
        bound_b = bound_a + interval;
        y2 = modeBisection(bound_b,mode_);
        
        if ((y1*y2 < 0) || (y1*y2 == 0))
        {
            // printf("\n y1= %f , y2= %f\n",y1,y2);
            // printf("\nEntro HEREEEE!!!!!!!!!\n");
            points_sign_change.pa = bound_a;
            points_sign_change.pb = bound_b;
            vector_sign_change.push_back(points_sign_change);
        }           
        // printf("\n a= %f , b= %f , y1= %f , y2= %f\n",bound_a,bound_b,y1,y2);
        // printf("\n y2 * y1 = %f\n",y2*y1);
        bound_a = bound_b;
    }
    printf("vector size = %i\n",vector_sign_change.size());

    for (int i = 0; i < vector_sign_change.size(); i++)
    {  
        points_interval print_pts_;
        print_pts_.pa = vector_sign_change[i].pa;
        print_pts_.pb = vector_sign_change[i].pb;
        if (mode_ == 0)
            printf ("Value point save for phi(x) xa = %g  xb= %g\n",print_pts_.pa,print_pts_.pb);
        if (mode_ == 1)
            printf ("Value point save for catenary_A(x) xa = %g  xb= %g\n",print_pts_.pa,print_pts_.pb);
        if (mode_ == 2)
            printf ("Value point save for Catenary_B(x) xa = %g  xb= %g\n",print_pts_.pa,print_pts_.pb);
        
    }
    if (vector_sign_change.size() < 1 )
    {
        if (mode_ == 0)
            printf ("\nPROBLEM WITHOUT SOLUTION FOR phi() !!! Interval [a-b] doesn't enclose a root \n");
        if (mode_ == 1)
            printf ("\nPROBLEM WITHOUT SOLUTION FOR Catenary_A() !!! Interval [a-b] doesn't enclose a root \n");
        if (mode_ == 2)
            printf ("\nPROBLEM WITHOUT SOLUTION FOR Catenary_B() !!! Interval [a-b] doesn't enclose a root \n");
    }
}