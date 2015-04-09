#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <fstream>
#include <cassert>
//#include <functional>
//#define _DEBUG

typedef double value_type;
typedef size_t size_type;

class Sphere
{
	
	private:  
		value_type x_, y_; // coordinates of sphere
		const value_type R_; // Radius of spheres
	
	public:
		Sphere(value_type x = 0.,
		       value_type y = 0.,
		       value_type R = 1.)
		: x_(x), y_(y), R_(R)
		{};
	      	
		value_type get_x_() const
		{
			return x_;
		}
		value_type get_y_() const
		{
			return y_;
		}
		value_type get_R_() const
		{
			return R_;
		}

		void set_x(value_type x){
			x_ = x;
		}

		void set_y(value_type y){
			y_ = y;
		}
		

		
};

/// should be binned for areas around each sphere (devide the areas around one sphere by 512 equal areas and 
/// and count number of spheres with center between the two borders for each area!)
/// needs still to be implemented!
value_type distance(const Sphere &x,const Sphere &y)
{

	/// I have to consider the Boundary conditions here! If distance  is bigger than L/2. then the real distance
	/// is (since periodic boundary conditions) dist. - L (since measure the other way around!)
	/// Also: Just use (x1-x2)^² instead of the sqrt of it, since the area calculation are needed
	value_type distance = std::sqrt(
			(x.get_x_() - y.get_x_())*(x.get_x_() - y.get_x_()) +
			(x.get_y_() - y.get_y_())*(x.get_y_() - y.get_y_()) 
			);
	return distance;
}

/// Try and take a step if no ovelapping with other spheres is assured!
void try_step(Sphere &Entry, std::vector<Sphere> &Space, value_type alpha)
{

	value_type trial_x, trial_y, xold, yold;

	std::random_device rd;
	std::mt19937 mt(rd()); // create an engine, alternatively use (time(NULL))
	std::uniform_real_distribution<double> ureal_d(0.,1.);

	assert( alpha < Entry.get_R_());
	/// adding random x,y value to the sphere's coords
	trial_x = ureal_d(mt) * alpha;
	trial_y = ureal_d(mt) * alpha;

	/// saving old values (in case of rejection!)
	xold = Entry.get_x_();
	yold = Entry.get_y_();

	/// set new values of sphere (trial values)
	trial_x = std::fmod((Entry.get_x_() + trial_x), 1.) ;
	trial_y = std::fmod((Entry.get_y_() + trial_y), 1.) ;
	Entry.set_x(trial_x);
	Entry.set_y(trial_y);

	if(
			std::any_of(Space.begin(),Space.end(), [&] (Sphere &n) {
				/// checking if any Sphere is overlapping with the current sphere's new position
				if ((n.get_x_() != trial_x) && (n.get_y_() != trial_y))
				{
				return (distance(Entry, n) > 2*Entry.get_R_());
				}
				else return true;
				})
	  )
	{} /// accepting the position
	else
	{
		/// return to old position
		Entry.set_x(xold);
		Entry.set_y(yold);
	}

}
	



/// call above functions to make a MC sweep
void sweep(std::vector<Sphere> &Space, value_type alpha)
{

	/// make step (or step-trial) for all spheres
	std::for_each(Space.begin(), Space.end(),[&](Sphere &n){try_step(n , Space, alpha);});
}



int main(int argc, char* argv[])
{
 
    if (argc < 3) 
    {
        std::cerr << "Usage: " << argv[0] << " Nx Ny r" << std::endl;
        return 1;
    }
    const value_type nx  = std::stod(argv[1]);
    const value_type ny  = std::stod(argv[2]);
    const value_type r   = std::stoul(argv[3]);
   

    std::vector<Sphere> Space;
    value_type xspacing = 1/nx;
    value_type yspacing = std::sqrt(3.) * xspacing/ 2.;
    
    /// fill array with Spheres at the points of a triangular lattice
    for(int i=0;i<ny;i++)
    {
	    for(int j=0;j<nx;j++)
	    {
		    Space.push_back(Sphere(j*xspacing,i*yspacing,r));
	    }
    }
    

#ifndef _DEBUG 
    for(auto& it : Space)
    {
	    //std::cout << "x: " << it.get_x_() << "  y: " << it.get_y_() << "  z: " << it.get_z_() << std::endl;
	    std::cout << it.get_x_() << "    " << it.get_y_() << "    " << std::endl;
    }
#endif
  
    return 0;
}

