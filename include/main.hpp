#pragma once

#define SOLVER  SOLVER_POSITION_VERLET
#define ENGINE  ENGINE_DIRECT_FORCE
#define COUNT   1000
#define TIME    1000.0
#define STEP    0.1

// Either use number of logical cores or hardware threads
// Not using hyperthreading can give better performance on specific systems
// Set to -1 if you want the max number of threads to be automatically determined
#define THREADS -1