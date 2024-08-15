//
//  main.cpp
//  Elasticity
//
//  Created by Wim van Rees on 2/15/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#include "common.hpp"
#include "ArgumentParser.hpp"

#include "Sim.hpp"
#include "Sim_Knit.hpp"

#ifdef USETBB
#include "tbb/task_scheduler_init.h"
#endif

int main(int argc,  const char ** argv)
{
#ifdef USETBB
    const int num_threads = tbb::task_scheduler_init::default_num_threads();
    tbb::task_scheduler_init init(num_threads);
    std::cout << "Starting TBB with " << num_threads << " threads " << std::endl;
#endif

    BaseSim * sim = nullptr;

    ArgumentParser parser(argc,argv);

    parser.initialize();

    sim = new Sim_Knit(parser);

    sim->init();
    sim->run();

    delete sim;

    parser.finalize();

    return 0;
}
