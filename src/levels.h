#ifndef LEVELS_H_INCLUDED
#define LEVELS_H_INCLUDED

#include "defs.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPRESSION LEVELS FOR DNA
//


#define LEVEL_1 " 1: -m 13:1:0:1:0.97/0:0:0 "
#define LEVEL_2 " 2: -m 4:1:1:0:0.78/0:0:0 -m 13:50:1:1:0.96/0:0:0 "
#define LEVEL_3 " 3: -m 3:1:0:0:0.80/0:0:0 -m 4:1:1:0:0.84/0:0:0 -m 13:50:1:1:0.94/2:15:0.95 "
#define LEVEL_4 " 4: -m 4:1:0:0:0.80/0:0:0 -m 17:200:1:10:0.95/5:10:0.95 "
#define LEVEL_5 " 5: -m 4:1:0:0:0.82/0:0:0 -m 6:1:1:0:0.72/0:0:0 -m 13:50:1:1:0.95/2:15:0.95 "
#define LEVEL_6 " 6: -m 4:1:0:0:0.88/0:0:0 -m 6:1:1:0:0.76/0:0:0 -m 13:50:1:1:0.95/2:15:0.95 "
#define LEVEL_7 " 7: -m 4:1:1:0:0.90/0:0:0 -m 6:1:1:0:0.79/0:0:0 -m 8:1:1:0:0.91/0:0:0 -m 13:10:1:0:0.94/1:20:0.94 -m 16:200:1:5:0.95/4:15:0.95 "
#define LEVEL_8 " 8: -m 4:1:1:0:0.90/0:0:0 -m 6:1:1:0:0.80/0:0:0 -m 13:10:1:0:0.95/1:20:0.94 -m 16:100:1:5:0.95/3:15:0.95 "
#define LEVEL_9 " 9: -m 4:1:1:0:0.91/0:0:0 -m 6:1:1:0:0.82/0:0:0 -m 13:10:1:0:0.95/1:20:0.94 -m 17:100:1:8:0.95/3:15:0.95 "
#define LEVEL_10 " 10: -m 1:1:0:0:0.90/0:0:0 -m 3:1:0:0:0.90/0:0:0 -m 6:1:1:0:0.82/0:0:0 -m 9:10:0:0:0.9/0:0:0 -m 11:10:0:0:0.9/0:0:0 -m 13:10:1:0:0.9/0:20:0.94 -m 17:100:1:8:0.89/5:10:0.9 "
#define LEVEL_11 " 11: -m 4:1:1:0:0.91/0:0:0 -m 6:1:1:0:0.82/0:0:0 -m 13:10:1:0:0.95/1:20:0.94 -m 17:100:1:15:0.95/3:15:0.95 "
#define LEVEL_12 " 12: -m 1:1:0:0:0.9/0:0:0 -m 3:1:0:0:0.9/0:0:0 -m 6:1:1:0:0.85/0:0:0 -m 9:10:0:0:0.9/0:0:0 -m 11:10:0:0:0.9/0:0:0 -m 13:50:1:0:0.9/0:20:0.94 -m 17:100:1:20:0.9/3:10:0.9 "
#define LEVEL_13 " 13: -m 4:1:0:1:0.8/0:0:0 -m 5:10:1:1:0.90/0:0:0 -m 9:1:0:1:0.90/0:0:0 -m 12:20:1:1:0.94/0:0:0 -m 14:50:1:1:0.95/0:0:0 -m 16:200:1:10:0.95/1:50:0.95 -m 20:500:1:30:0.95/2:20:0.95 "
#define LEVEL_14 " 14: -m 2:1:0:1:0.75/0:0:0 -m 4:1:0:1:0.85/0:0:0 -m 5:10:1:1:0.90/0:0:0 -m 9:1:0:1:0.90/0:0:0 -m 12:20:1:1:0.94/0:0:0 -m 14:50:1:1:0.95/0:0:0 -m 16:200:1:10:0.95/1:50:0.95 -m 20:500:1:40:0.95/2:20:0.95 "
#define LEVEL_15 " 15: -m 1:1:1:1:0.90/0:0:0 -m 3:1:1:1:0.90/0:0:0 -m 6:1:1:1:0.85/0:0:0 -m 9:1:1:1:0.91/0:0:0 -m 13:20:1:1:0.92/0:0:0 -m 14:50:1:1:0.92/0:0:0 -m 16:200:1:10:0.92/1:50:0.92 -m 20:500:1:100:0.92/5:10:0.92 "

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char    *GetLevels  (uint8_t);
void    PrintLevels (void);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

