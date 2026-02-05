#pragma once

#define SOLVER  SOLVER_PEFRL
#define ENGINE  ENGINE_DIRECT_FORCE
#define COUNT   1000
#define TIME    1000.0
#define STEP    0.1

// Either use number of logical cores or hardware threads
// Not using hyperthreading can give better performance on specific systems
#define THREADS 16