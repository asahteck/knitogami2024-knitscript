#ifndef Sim_Knit_hpp
#define Sim_Knit_hpp

#include "Sim.hpp"
#include "Mesh.hpp"

class Sim_Knit : public Sim<Mesh>
{
    typedef Mesh tMesh;

protected:
    void runKnit();

public:
    Sim_Knit(ArgumentParser & parser):
    Sim<tMesh>(parser)
    {}

    void init() override;
    void run() override;
};

#endif /* Sim_Knit_hpp */
