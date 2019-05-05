// Copyright (C) 2019 Ksenia Guseva <ksenia@skewed.de>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.



#include <cmath>
#include <iostream>
#include <iomanip>

#include <ext/numeric>
using __gnu_cxx::power;
#include <vector>
#include <utility>

#include <boost/python.hpp>

#include "numpy_bind.hh"


using namespace std;
using namespace boost;



int check_circle_c(double xp, double yp)
{
    if(sqrt(pow(xp,2) + pow(yp,2)) < 1.)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}



void get_interpolation(int Nx, int Ny, python::object oxf, python::object oyf,
                       double init_x, double init_y, python::object ox_range,
                       python::object oy_range, double dx, double dy,
                       python::object oMat, python::object oMat_new)
    

{
    multi_array_ref<double, 2> Mat = get_array<double,2>(oMat);
    multi_array_ref<double, 2> Mat_new = get_array<double,2>(oMat_new);
    multi_array_ref<double, 2> xf = get_array<double,2>(oxf);
    multi_array_ref<double, 2> yf = get_array<double,2>(oyf);
    multi_array_ref<double, 1> x_range = get_array<double,1>(ox_range);
    multi_array_ref<double, 1> y_range = get_array<double,1>(oy_range);

    int i, j;

    for (i = 1; i < Nx; ++i)
    {
        for (j = 1; j < Ny-1; ++j)
        {
            double xp = xf[i][j];
            double yp = yf[i][j];

            int ch_ci = check_circle_c(xp, yp);
            if(ch_ci == 1)
            {
                continue;
            }
            
            int i1 = int(floor((xp - init_x)/dx));
            int j1 = int(floor((yp - init_y)/dy));
            int i2 = (i1 + 1);
            int j2 = (j1 + 1);

     
            
                
            if(i1 < 0)    i1 = 0;
            if(i2 < 0)    i2 = 0;
            if(j1 < 0)    j1 = 0;
            if(j2 < 0)    j2 = 0;
             
            if(i1 > Nx-1)  i1 = Nx-1;
            if(i2 > Nx-1)  i2 = Nx-1;
            if(j1 > Ny-1)  j1 = Ny-1;
            if(j2 > Ny-1)  j2 = Ny-1;
            
            int ch_c[4];
            for (int aux = 1; aux < 4; ++aux)
            {
                ch_c[aux] = 0;
                
            }
            
                   
            double x1 = x_range[i1];
            double x2 = x_range[i2];
            double y1 = y_range[j1];
            double y2 = y_range[j2];
            

  
            
            double c11 =  Mat[i1][j1];
            double c12 =  Mat[i1][j2];
            double c21 =  Mat[i2][j1];
            double c22 =  Mat[i2][j2];
            
                    
            ch_c[0] = check_circle_c(x1, y1);
            if(ch_c[0] == 1)  c11 = c21;
            
            ch_c[1] = check_circle_c(x1, y2);
            if(ch_c[1] == 1) c12 = c22;
            ch_c[2] = check_circle_c(x2, y1);
            if(ch_c[2] == 1) c21 = c11;
            ch_c[3] = check_circle_c(x2, y2);
            if(ch_c[3] == 1) c22 = c12;
            
                        
            if(ch_c[0]==1 && ch_c[1]==1) xp = x2; 
            if(ch_c[2]==1 && ch_c[3]==1) xp = y1;
            if(ch_c[0]==1 && ch_c[1]==1) xp = x2;
            
            if(ch_c[0]==1 && ch_c[2]==1) yp = y2;
            if(ch_c[1]==1 && ch_c[3]==1) yp = y1;
            


            if (j1 == j2)
            {    
                if (i1 == i2)
                {
                    Mat_new[i][j] = c11;
                    continue;
                }
                else
                {        
                    Mat_new[i][j] = (c11*(x2 - xp) + c21*(xp - x1))/dx;
                    continue;
                }
            }
            
            if (i1 == i2)
            {
                if (j1 == j2)
                {
                    Mat_new[i][j] = c11;
                    continue;
                }
                else
                {
                    Mat_new[i][j] = (c11*(y2 - yp) + c12*(yp - y1))/dy;
                    continue;
                }
                
            }
            Mat_new[i][j] = (c11*(x2 - xp)*(y2 - yp) + c12*(x2 -xp)*(yp - y1)
                             + c21*(xp - x1)*(y2 - yp) + c22*(xp - x1)*(yp - y1))/(dx*dy);
         
            
        }
    }   
}


void diffusion(int Nx, int Ny, double dx, double dy,
               python::object oMat, python::object oMat_new, double D, double dt,
               python::object ox, python::object oy)
    

{
    multi_array_ref<double, 2> Mat = get_array<double,2>(oMat);
    multi_array_ref<double, 2> Mat_new = get_array<double,2>(oMat_new);
    multi_array_ref<double, 2> x = get_array<double,2>(ox);
    multi_array_ref<double, 2> y = get_array<double,2>(oy);


    int i, j;
    double Mjp, Mjm, Mip, Mim;
    int aux;
    
    for (i = 2; i < Nx-1; ++i)
    {
        for (j = 2; j < Ny-1; ++j)
        {
            aux = check_circle_c(x[i][j], y[i][j]);
            if(aux != 1)
            {
                Mjp = Mat[i][j+1];
                Mjm = Mat[i][j-1];
                Mip = Mat[i+1][j];
                Mim = Mat[i-1][j];
           
                aux = check_circle_c(x[i][j+1], y[i][j+1]);
                if(aux == 1)  Mjp = Mat[i][j];
                aux = check_circle_c(x[i][j-1], y[i][j-1]);
                if(aux == 1)  Mjm = Mat[i][j];
                aux = check_circle_c(x[i+1][j], y[i+1][j]);
                if(aux == 1)  Mip = Mat[i][j];
                aux = check_circle_c(x[i-1][j], y[i-1][j]);
                if(aux == 1)  Mim = Mat[i][j];

           
                Mat_new[i][j] += dt*D*(Mjp - 2*Mat[i][j] +  Mjm)/(dy*dy);
                Mat_new[i][j] += dt*D*(Mip - 2*Mat[i][j] +  Mim)/(dx*dx);
            }   
        }
    }
}


void get_island(int Nx, int Ny,python::object oMat, python::object ox, python::object oy)
{
    multi_array_ref<double, 2> Mat = get_array<double,2>(oMat);
    multi_array_ref<double, 2> x = get_array<double,2>(ox);
    multi_array_ref<double, 2> y = get_array<double,2>(oy);

    int i, j;
    
    for (i = 1; i < Nx; ++i)
    {
        for (j = 1; j < Ny; ++j)
        {
            if(sqrt(pow(x[i][j], 2) + pow(y[i][j], 2)) < 1.0)
            {
                Mat[i][j] = 0;
            }
        }
    }

                    
}


python::object get_gridmean(int Nx, int Ny,python::object oMat, python::object ox, python::object oy)
{
    multi_array_ref<double, 2> Mat = get_array<double,2>(oMat);
    multi_array_ref<double, 2> x = get_array<double,2>(ox);
    multi_array_ref<double, 2> y = get_array<double,2>(oy);

    int i, j;

    double av = 0;
    double min = Mat[1][1];
    double max = 0;
    int num = 0;
    
    for (i = 1; i < Nx; ++i)
    {
        for (j = 1; j < Ny; ++j)
        {
            if(sqrt(pow(x[i][j], 2) + pow(y[i][j], 2)) > 1.0)
            {
                if (Mat[i][j] > max) max = Mat[i][j];
                if (Mat[i][j] < min) min = Mat[i][j];
                av = av + Mat[i][j];
                
                num = num + 1;
                
            }
        }
    }

    return  python::make_tuple(av/num, min, max);
}





BOOST_PYTHON_MODULE(upw_loops)
{
    using namespace boost::python;

    def("get_inerpolation_c", get_interpolation);
    def("diffusion", diffusion);
    def("get_island", get_island);
    def("get_gridmean", get_gridmean);
};
