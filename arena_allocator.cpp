#include <cstdlib>

struct Arena {
    void **memory_blocks;
    uint number_of_blocks;
    uint block_size; 
    uint current_block_idx;
    uint current_byte_idx;
};

Arena *ArenaAllocate(uint block_size, uint starting_max_blocks) {
    Arena *result = (Arena *)malloc(sizeof(Arena));
    void **memory_blocks = (void **)malloc(starting_max_blocks * sizeof(void *));
    void *starting_block = malloc(block_size);
    memory_blocks[0] = starting_block;
    result->memory_blocks = memory_blocks;
    result->number_of_blocks = starting_max_blocks;
    result->block_size = block_size;
    result->current_block_idx = 0;
    result->current_byte_idx = 0;
    return result;
}

Arena *ArenaAllocateDefault() {
    const uint default_block_size = 1048576; // let's do 1MB
    const uint default_max_blocks = 10;
    return ArenaAllocate(default_block_size, default_max_blocks);
}


void ArenaGrow(Arena *arena) {
    if (arena->current_block_idx + 1 == arena->number_of_blocks) {
        arena->memory_blocks = (void **)realloc(arena->memory_blocks, 2 * arena->number_of_blocks * sizeof(void *));
    }
    arena->memory_blocks[++arena->current_block_idx] = malloc(arena->block_size * sizeof(void *));
    arena->current_byte_idx = 0;
}

void *ArenaGetMemory(uint number_of_bytes, Arena *arena) {
    if (number_of_bytes > arena->block_size) {
        return nullptr;
    }
    if (arena->current_byte_idx + number_of_bytes >= arena->block_size) {
        ArenaGrow(arena);
    }
    void *current_memory_block = arena->memory_blocks[arena->current_block_idx];
    void *result = (char *)current_memory_block + arena->current_byte_idx;
    arena->current_byte_idx += number_of_bytes;
    return result;
}


void ArenaReset(Arena *arena) {
    arena->current_block_idx = 0;
    arena->current_byte_idx = 0;
}

void ArenaDeallocate(Arena *arena) {
    for (uint i = 0; i < arena->current_block_idx; ++i) {
        free(arena->memory_blocks[i]);
    }
    free(arena->memory_blocks);
    free(arena);
}