// Example codes for HPC course
// (c) 2012-2014 Jan Gukelberger, Andreas Hehn, ETH Zurich

#include <array>
#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include <iterator>
#include <cassert>


typedef std::size_t size_type;
typedef double scalar_type;
typedef std::array<scalar_type,2> position;


std::ostream& operator<<(std::ostream& os, const position& x)
{
    std::copy(x.begin(),x.end(),std::ostream_iterator<scalar_type>(os,"\t"));
    return os;
}

struct potential
{
    potential(scalar_type rm, scalar_type epsilon, scalar_type rc=0)
    :   rm2_(rm*rm)
    ,   eps_(epsilon)
    ,   rc2_(1e10*rm)
    ,   shift_(0)
    {
        // default cut-off radius
        if(rc <= 0 )    rc = 2.5*rm/std::pow(2,1/6.);

        position x = {{}};
        position y = {{rc}};
        assert( x[0] == 0  && x[1] == 0 );
        assert( y[1] == 0 );
        shift_ = -(*this)(x,y,position());
        rc2_ = rc*rc;
        std::cout << "# Potential shift -V(rc=" << rc << ")=" << shift_ << std::endl;
    }

    /// potential V(x,y) considering periodic boundaries at extent
    scalar_type operator()(const position& x, const position& y, const position& extent) const
    {
        // TODO
    }

    /// compute the Lennard-Jones force F(x,y) which particle y exerts on x and add it to f,
    /// considering periodic boundaries at extent
    void add_force(position& f, const position& x, const position& y, const position& extent) const
    {
        // TODO
    }

    scalar_type cutoff_radius() const { return std::sqrt(rc2_); }

private:
    scalar_type rm2_;   // r_m^2
    scalar_type eps_;   // \epsilon
    scalar_type rc2_;   // cut-off radius r_c^2
    scalar_type shift_; // potential shift -V(r_c)
};


class simulation
{
public:
    /// Initialize simulation in rectangular box with corners (0,0) and extent.
    /// Initial positions and velocities are given as x, v.
    simulation(const position& extent, const potential& pot,
               const std::vector<position>& x, const std::vector<position>& v )
    :   extent_(extent)
    ,   potential_(pot)
    ,   x_(x)
    ,   v_(v)
    ,   a_(x.size())
    {
        assert( x.size() == v.size() );
        calculate_forces(a_,x_);
    }

    /// evolve the system for [steps] time steps of size [dt]
    void evolve(scalar_type dt, size_type steps)
    {
        using std::swap;
        configuration aold(a_);

        for( size_type s = 0; s < steps; ++s )
        {
            update_positions(x_,v_,a_,dt);
            swap(a_,aold);
            calculate_forces(a_,x_);
            update_velocities(v_,aold,a_,dt);
        }
    }

    /// print the current configuration
    void print_config() const
    {
        for( size_type i = 0; i < x_.size(); ++i )
            std::cout << x_[i] << v_[i] << a_[i] << std::endl;
        std::cout << std::endl;
    }

    /// calculate kinetic and potential energy of the current configuration
    std::pair<scalar_type,scalar_type> measure_energies() const
    {
        scalar_type epot = 0, ekin = 0;

        // TODO

        return std::make_pair(ekin,epot);
    }


private:
    typedef std::vector<position> configuration;

    void update_positions(configuration& x, const configuration& v, const configuration& a, scalar_type dt)
    {
        // TODO: Verlet step for positions x

        // TODO: enforce periodic boundaries
    }

    void update_velocities(configuration& v, const configuration& aold, const configuration& a, scalar_type dt)
    {
        // TODO: Verlet step for velocities v
    }

    void calculate_forces(configuration& a, const configuration& x)
    {
        // TODO: calculate forces on particles -> a[i]
    }


    position extent_; /// system extent along each dimension
    potential potential_;

    configuration x_; /// particle positions
    configuration v_; /// particle velocities
    configuration a_; /// forces on particles
};


std::vector<position> init_circle(const position& extent, size_type n)
{
    using std::sin;
    using std::cos;
    std::vector<position> config;
    position midpoint{{extent[0]/2, extent[1]/2}};
    for(size_type i=0; i < n; ++i)
        config.push_back({{midpoint[0]+0.9*midpoint[0]*sin(2*i*M_PI/n),midpoint[1]+midpoint[1]*0.9*cos(2*i*M_PI/n)}});
    return config;
}

int main(int argc, const char** argv)
{
    try
    {
        // get parameters from command line
        if( argc != 8 ){
            std::cerr << "Usage: " << argv[0] << " [box_length] [# particles] [r_m] [epsilon] [time step] [# time steps] [print steps]" << std::endl
                      << "    e.g.  " << argv[0] << " 1.0 100 0.05 5.0 1e-7 1000000 1000" << std::endl;
            return -1;
        }
        scalar_type box_length = std::atof(argv[1]);
        size_type   particles  = std::atoi(argv[2]);
        scalar_type rm         = std::atof(argv[3]);
        scalar_type eps        = std::atof(argv[4]);
        scalar_type dt         = std::atof(argv[5]);
        size_type   steps      = std::atoi(argv[6]);
        size_type   printsteps = std::atoi(argv[7]);

        // init potential
        potential pot(rm,eps);

        // init particle positions
        position extent;
        std::fill(extent.begin(),extent.end(),box_length);
        std::vector<position> x = init_circle(extent,particles);
        particles = x.size();
        std::cout << "# nparticles = " << particles << std::endl;

        std::vector<position> v(x.size(), position{{0,0}});

        // init and run simulation for [steps] steps
        simulation sim(extent,pot,x,v);
        for( size_type i = 0; i < steps/printsteps; ++i )
        {
            // print energies and configuration every [printsteps] steps
            scalar_type ekin,epot;
            std::tie(ekin,epot) = sim.measure_energies();
            std::cout << "# E/N=" << (ekin+epot)/particles << ", Ekin/N=" << ekin/particles
                      << ", Epot/N=" << epot/particles << std::endl;
#ifdef PRINT_CONFIGS
            sim.print_config();
#endif //PRINT_CONFIGS
            // run for [printsteps] steps
            sim.evolve(dt,printsteps);
        }
    }
    catch( std::exception& e )
    {
        std::cerr << e.what() << std::endl;
        throw;
    }
}

