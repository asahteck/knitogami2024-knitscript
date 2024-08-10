//
//  Sim_Knit.hpp
//  Elasticity
//
//  Created by Wim van Rees on 9/14/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Sim_Knit_hpp
#define Sim_Knit_hpp

#define KNIT_BILAYER_MESH

#include "Sim.hpp"
#include "Mesh.hpp"

#ifdef KNIT_BILAYER_MESH
class Sim_Knit : public Sim<BilayerMesh>
#else
class Sim_Knit : public Sim<Mesh>
#endif
{
#ifdef KNIT_BILAYER_MESH
    typedef BilayerMesh tMesh;
#else
    typedef Mesh tMesh;
#endif

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
