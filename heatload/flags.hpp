#ifndef FLAGS_HPP
#define FLAGS_HPP

class Flags 
{
//private:
public:
    static const int is_to_write =1;   // need to write
    static const int is_in_init  =2;   // inside of separatrix initially -- Not used here
    static const int is_written  =4;   // inside of separatrix before push -- Not used here
    static const int is_escaped  =8;   // why write (is_to_write=on)- just escaped from separatrix
    static const int is_divertor =16;  // why write - just hit the divertor
    static const int is_outboard =32;  // when write - is this outboard?
    static const int was_inside  =64;  // aux flag to indicate if the particle was inside before push

    // actual data
    bool to_write;
    bool in_init;
    bool written;
    bool escaped;
    bool divertor;
    bool outboard;
    bool inside;

    Flags(int);
    ~Flags();
};

inline Flags::Flags(int flag_in)
{
    to_write = (flag_in & is_to_write) >0;
    in_init  = (flag_in & is_in_init) >0;
    written  = (flag_in & is_written) >0;
    escaped  = (flag_in & is_escaped) >0;
    divertor = (flag_in & is_divertor) >0;
    inside   = (flag_in & was_inside) >0;
}

inline Flags::~Flags()
{
}

#endif
