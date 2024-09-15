/**
 * (C) Copyright 2024 Aion Feehan. All Rights Reserved.
 * 
 * This software is provided 'as-is', without any express or implied warranty. In no event will the authors
 * be held liable for any damages arising from the use of this software.
 * 
 */
#include <cstdlib>

struct Arena;
Arena *ArenaAllocate(uint block_size, uint starting_max_blocks);
Arena *ArenaAllocateDefault();
void ArenaGrow(Arena *arena);
void *ArenaGetMemory(uint number_of_bytes, Arena *arena);
void ArenaReset(Arena *arena);
void ArenaDeallocate(Arena *arena);