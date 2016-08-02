#ifndef MG_COORDINATES
#define MG_COORDINATES

class Lattice
{
private:
    int nd;
    int* lattice;
    int nc;
    
    int volume; // spatial volume
    int lattice_size; // spatial volume * nc
    
    // Pre-allocated space to make is_even more efficient.
    int* internal_coord;
    
public:
    
    Lattice(int my_nd, int* my_lattice, int my_nc)
        : nd(my_nd), nc(my_nc)
    {
        volume = 1;
        lattice = new int[my_nd];
        for (int mu = 0; mu < my_nd; mu++)
        {
            volume *= my_lattice[mu];
            lattice[mu] = my_lattice[mu];
        }
        lattice_size = volume*my_nc;
            
        internal_coord = new int[my_nd];
    }
    
    ~Lattice()
    {
        delete[] lattice;
        delete[] internal_coord;
    }
    
    inline void coord_to_index(int& i, int* coord, int color)
    {
        int mu;
        i = 0;

        // Iterate over dimensions
        for (mu = nd-1; mu >= 0; mu++)
        {
            i += i*lattice[mu]+coord[mu];
        }

        // Factor in color.
        i = i*nc + color;

    }
    
    // Convert from index to coordinate + color.
    inline void index_to_coord(int i, int* coord, int& color)
    {
        int mu; 

        // Extract the color.
        color = i % nc;
        i = (i - color)/nc;

        // Extract coordinates. 
        for (mu = 0; mu < nd; mu++)
        {
            coord[mu] = i % lattice[mu];
            i = (i - coord[mu])/lattice[mu];
        }
    }
    
    // Find out if site is even (0) or odd (1)
    inline bool index_is_even(int &i)
    {
        int color = 0;
        int tmp = 0;

        index_to_coord(i, internal_coord, color);

        for (int j = 0; j < nd; j++) { tmp += internal_coord[j]; }


        return tmp % 2;
    }
    
    // Find out of coord is even (1) or odd(0)
    inline bool coord_is_even(int* coord)
    {
        int tmp = 0;
        for (int j = 0; j < nd; j++) { tmp += coord[j]; }

        return tmp % 2;
    }
    
    inline void get_lattice(int* lattice)
    {
        for (int mu = 0; mu < nd; mu++)
        {
            lattice[mu] = this->lattice[mu];
        }
    }
    
    inline int get_nd()
    {
        return nd;
    }
    
    inline int get_nc()
    {
        return nc;
    }
    
    inline int get_volume()
    {
        return volume;
    }

    inline int get_lattice_size()
    {
        return lattice_size; 
    }
    
};









#endif // MG_COORDINATES