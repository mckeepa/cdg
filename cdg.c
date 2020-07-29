#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/* Stolen from cop/cop_sort.h */
#define COP_SORT_INSERTION(fn_name_, data_type_, comparator_macro_) \
void fn_name_(data_type_ *inout, size_t nb_elements) \
{ \
	size_t i, j; \
	assert(nb_elements >= 2); \
	for (i = 1; i < nb_elements; i++) { \
		data_type_ tmp; \
		tmp = inout[i]; \
		for (j = i; j && (comparator_macro_(tmp, inout[j-1])); j--) \
			inout[j] = inout[j-1]; \
		if (j != i) \
			inout[j] = tmp; \
	} \
}

/* Set structure. This is poorly named.
 *
 * This is a container of bits representing a large unsigned number. It
 * contains SET_BITS bits held using SET_ELEMENTS unsigned 64-bit numbers.
 * Bit 0 is the most significant bit. Bit SET_BITS-1 is the least significant
 * bit. */

#define SET_ELEMENTS (4)
#define SET_BITS     (64*SET_ELEMENTS)

struct bigu {
    uint_fast64_t groups[SET_ELEMENTS];
};

/* Compare two big numbers. Non-zero if p_1 is bigger than p_2. */
static int bigu_is_gt(const struct bigu *p_1, const struct bigu *p_2) {
    unsigned i;
    for (i = 0; i < SET_ELEMENTS; i++) {
        if (p_1->groups[i] > p_2->groups[i])
            return 1;
        if (p_1->groups[i] < p_2->groups[i])
            break;
    }
    return 0;
}

#define BIGU_CMP(a_, b_) (bigu_is_gt(&(a_), &(b_)))
static COP_SORT_INSERTION(bigu_sort, struct bigu, BIGU_CMP)
#undef BIGU_CMP

static void bigu_set(struct bigu *p, unsigned bit) {
    assert(bit < SET_BITS);
    p->groups[bit / 64] |= 1ull << (63 - (bit & 63));
}
static int bigu_get(const struct bigu *p, unsigned bit) {
    assert(bit < SET_BITS);
    return (p->groups[bit / 64] & (1ull << (63 - (bit & 63)))) ? 1 : 0;
}
static void bigu_zero(struct bigu *p) {
    unsigned i;
    for (i = 0; i < SET_ELEMENTS; i++)
        p->groups[i] = 0;
}
static void bigu_or(struct bigu *p, const struct bigu *q) {
    unsigned i;
    for (i = 0; i < SET_ELEMENTS; i++) {
        assert((p->groups[i] & q->groups[i]) == 0);
        p->groups[i] |= q->groups[i];
    }
}
/* Count leading set bits. This finds the index of the first zeroed bit in the
 * big number. Undefined if all bits are set. */
static unsigned bigu_cls(struct bigu *p) {
    unsigned i = 0;
    while (bigu_get(p, i))
        i++;
    return i;
}
/* Returns zero if the logical-and of the bits in the two big numbers is
 * zero. */
static int bigu_and_zero(const struct bigu *p, const struct bigu *q) {
    unsigned i;
    for (i = 0; i < SET_ELEMENTS; i++)
        if (p->groups[i] & q->groups[i])
            return 1;
    return 0;
}
/* Returns zero if the two big numbers are different. */
static int bigu_is_equal(const struct bigu *p, const struct bigu *q) {
    unsigned i;
    for (i = 0; i < SET_ELEMENTS; i++)
        if (p->groups[i] != q->groups[i])
            return 0;
    return 1;
}

static void bigu_print(const struct bigu *p) {
    unsigned i;
    for (i = 0; i < SET_ELEMENTS; i++)
        printf("%016llx", p->groups[i]);
    printf("\n");    
}

#define MAX_LENGTH (7)
#define MAX_GRID   (16)

/* Finds a combination of the given p_sets that is a unique covering of
 * p_complete. Returns zero on failure. Otherwise returns the number of
 * elements in p_shape_stack array that form the unique cover. */
static
unsigned
recursive_search
    (struct bigu        current_set
    ,const struct bigu *p_complete
    ,const struct bigu *p_sets
    ,unsigned          *p_level_offsets
    ,unsigned           level
    ,struct bigu       *p_shape_stack
    ,unsigned           shape_stack_height
    ) {
    unsigned i;
    for (i = p_level_offsets[level]; i < p_level_offsets[level+1]; i++) {
        if (!bigu_and_zero(&(current_set), &(p_sets[i]))) {
            struct bigu next = current_set;
            unsigned u;
            bigu_or(&next, &(p_sets[i]));
            p_shape_stack[shape_stack_height++] = p_sets[i];
            if (bigu_is_equal(&next, p_complete))
                return shape_stack_height;
            u = recursive_search(next, p_complete, p_sets, p_level_offsets, bigu_cls(&next), p_shape_stack, shape_stack_height);
            if (u)
                return u;
            shape_stack_height--;
        }
    }
    return 0;
}

unsigned rng16(unsigned long *state) {
    unsigned long s = *state;
#define HULL_DOBEL_SATSIFYING_MAGIC (1u+4u*961748941u)
    *state = (s * HULL_DOBEL_SATSIFYING_MAGIC + 1u) & 0xFFFFFFFFu;
    return s >> 16;
}

void shuffle(unsigned *data, unsigned sz, unsigned long *rseed) {
    unsigned i;
    for (i = 0; i < 10000; i++) {
        unsigned x = rng16(rseed) % sz;
        unsigned y = rng16(rseed) % sz;
        unsigned tmp = data[x];
        data[x] = data[y];
        data[y] = tmp;
    }
}

void shuffle_sets(struct bigu *data, unsigned sz, unsigned long *rseed) {
    unsigned i;
    for (i = 0; i < 10000; i++) {
        unsigned x = rng16(rseed) % sz;
        unsigned y = rng16(rseed) % sz;
        struct bigu tmp = data[x];
        data[x] = data[y];
        data[y] = tmp;
    }
}

struct shape {
    unsigned width;
    unsigned height;
    unsigned cells[46]; /* whevever */
};


static const struct shape shapes[] =
/* 2 pieces */
{   {2, 1, {1, 1}}
/* 3 pieces */
,   {3, 1, {1, 1, 1}}
,   {2, 2, {1, 0
           ,1, 1}}
/* 4 pieces */
,   {4, 1, {1, 1, 1, 1}}
,   {2, 2, {1, 1
           ,1, 1}}
,   {3, 2, {1, 1, 1
           ,1, 0, 0}}
,   {3, 2, {0, 1, 1
           ,1, 1, 0}}
/* 5 pieces */
,   {3, 2, {1, 1, 1
           ,1, 1, 0}}
,   {3, 2, {1, 1, 1
           ,1, 0, 1}}

};

unsigned remove_uniques(struct bigu *p_sets, unsigned nset) {
    unsigned i, j;
    for (i = 0; i < nset; i++) {
        for (j = i+1; j < nset; j++) {
            if (bigu_is_equal(&(p_sets[i]), &(p_sets[j]))) {
                p_sets[j--] = p_sets[--nset];
            }
        }
    }
    return nset;
}

unsigned insert_shape_permutations(struct bigu *p_sets, unsigned grid_size, const struct shape *p_shape) {
    unsigned gx, gy, sx, sy;
    unsigned sidx = 0;
    if (p_shape->height >= grid_size)
        return 0;
    if (p_shape->width >= grid_size)
        return 0;
    sidx = 0;
    for (gy = 0; gy < grid_size + 1 - p_shape->height; gy++) {
        for (gx = 0; gx < grid_size + 1 - p_shape->width; gx++) {
            for (sy = 0; sy < 8; sy++)
                bigu_zero(&(p_sets[sidx+sy]));
            for (sy = 0; sy < p_shape->height; sy++) {
                for (sx = 0; sx < p_shape->width; sx++) {
                    if (p_shape->cells[sy*p_shape->width+sx]) {
                        bigu_set(&(p_sets[sidx+0]), (gy                  +sy)*grid_size+(gx                 +sx)          );
                        bigu_set(&(p_sets[sidx+1]), (gy+p_shape->height-1-sy)*grid_size+(gx                 +sx)          );
                        bigu_set(&(p_sets[sidx+2]), (gy+p_shape->height-1-sy)*grid_size+(gx+p_shape->width-1-sx)          );
                        bigu_set(&(p_sets[sidx+3]), (gy                  +sy)*grid_size+(gx+p_shape->width-1-sx)          );
                        bigu_set(&(p_sets[sidx+4]), (gy                  +sy)          +(gx                 +sx)*grid_size);
                        bigu_set(&(p_sets[sidx+5]), (gy+p_shape->height-1-sy)          +(gx                 +sx)*grid_size);
                        bigu_set(&(p_sets[sidx+6]), (gy+p_shape->height-1-sy)          +(gx+p_shape->width-1-sx)*grid_size);
                        bigu_set(&(p_sets[sidx+7]), (gy                  +sy)          +(gx+p_shape->width-1-sx)*grid_size);
                    }
                }
            }
            sidx += 8;
        }
    }
    return remove_uniques(p_sets, sidx);
}

#if 0
int count_solutions
    (const unsigned *p_grid
    ,unsigned        grid_size
    ,struct bigu      already_picked
    ,unsigned        row_set
    ,unsigned        col_set
    ) {




}
#endif

#define MAX_SEGMENT_SIZE (8)
#define OP_FIXED (0)
#define OP_SUM (1)
#define OP_SUB (2)
#define OP_MUL (3)
#define OP_DIV (4)

struct kkeq {
    int      operation;
    unsigned value;
    unsigned nb_segment;
    unsigned segment_indexes[MAX_SEGMENT_SIZE];

};



#define KKEQ_CMP(a_, b_) ((a_).nb_segment < (b_).nb_segment)
static COP_SORT_INSERTION(kkeq_sort, struct kkeq, KKEQ_CMP)
#undef KKEQ_CMP

void dump_puzzle(const struct kkeq *p_shapes, const unsigned *p_grid_values, unsigned nb_shapes, unsigned grid_size, int cheat) {
    unsigned i, j;
#define FLAG_TOP (1)
#define FLAG_LEFT (2)
#define FLAG_BOTTOM (4)
#define FLAG_RIGHT  (8)

    unsigned gridgroups[MAX_GRID*MAX_GRID];      /* grid_size * grid_size elements */
    unsigned gridborders[MAX_GRID*MAX_GRID];     /* grid_size * grid_size elements */
    char     gridtextdata[MAX_GRID*MAX_GRID][16]; /* u elements */
    memset(gridtextdata, 0, sizeof(gridtextdata));

    for (i = 0; i < nb_shapes; i++) {
        int x = 0;
        for (j = 0; j < p_shapes[i].nb_segment; j++) {
            gridgroups[p_shapes[i].segment_indexes[j]] = i;
        }
        switch (p_shapes[i].operation) {
        case OP_SUM:
            sprintf(gridtextdata[p_shapes[i].segment_indexes[0]], "%d+", p_shapes[i].value);
            break;
        case OP_SUB:
            sprintf(gridtextdata[p_shapes[i].segment_indexes[0]], "%d-", p_shapes[i].value);
            break;
        case OP_MUL:
            sprintf(gridtextdata[p_shapes[i].segment_indexes[0]], "%d*", p_shapes[i].value);
            break;
        case OP_DIV:
            sprintf(gridtextdata[p_shapes[i].segment_indexes[0]], "%d/", p_shapes[i].value);
            break;
        default:
            sprintf(gridtextdata[p_shapes[i].segment_indexes[0]], "%d", p_shapes[i].value);
            break;
        }
    }

    if (cheat) {
        for (i = 0; i < nb_shapes; i++) {
            for (j = 0; j < p_shapes[i].nb_segment; j++) {
                char buf[16];
                sprintf(buf, " = %d", p_grid_values[p_shapes[i].segment_indexes[j]]);
                strcat(gridtextdata[p_shapes[i].segment_indexes[j]], buf);
            }
        }
    }

    for (j = 0; j < grid_size; j++) {
        for (i = 0; i < grid_size; i++) {
            unsigned bflags = 0;
            if (i == 0 || gridgroups[j*grid_size+i] != gridgroups[j*grid_size+i-1])
                bflags |= FLAG_LEFT;
            if (j == 0 || gridgroups[j*grid_size+i] != gridgroups[(j-1)*grid_size+i])
                bflags |= FLAG_TOP;
            if (i == grid_size-1 || gridgroups[j*grid_size+i] != gridgroups[j*grid_size+i+1])
                bflags |= FLAG_RIGHT;
            if (j == grid_size-1 || gridgroups[j*grid_size+i] != gridgroups[(j+1)*grid_size+i])
                bflags |= FLAG_BOTTOM;
            gridborders[j*grid_size+i] = bflags;
        }
    }

    printf("<html><head><title>hello</title><style>table { border-collapse: collapse; } ");
    for (i = 0; i < 16; i++) {
        printf
            ("td.x%d { font-size: 10px; vertical-align: top; text-align: left; border-top: %dpx solid black; border-left: %dpx solid black; border-bottom: %dpx solid black; border-right: %dpx solid black; height: 40px; width: 40px;} "
            ,i
            ,(i & FLAG_TOP) ? 4 : 1
            ,(i & FLAG_LEFT) ? 4 : 1
            ,(i & FLAG_BOTTOM) ? 4 : 1
            ,(i & FLAG_RIGHT) ? 4 : 1
            );
    }
    printf("</style></head><body><table>");
    for (j = 0; j < grid_size; j++) {
        printf("<tr>");
        for (i = 0; i < grid_size; i++) {
            printf("<td class=x%d>%s</td>", gridborders[j*grid_size+i], gridtextdata[j*grid_size+i]);
        }
        printf("</tr>");
    }
    printf("</table></body></html>");
}

void create_random_grid_values(unsigned *grid_data, unsigned grid_size, unsigned long *rseed) {
    unsigned i, j;
    for (j = 0; j < grid_size; j++) {
        for (i = 0; i < grid_size; i++) {
            grid_data[j*grid_size+i] = 1 + ((i+j) % grid_size);
        }
    }
    for (i = 0; i < 10000; i++) {
        unsigned r1 = rng16(rseed) % grid_size;
        unsigned r2 = rng16(rseed) % grid_size;
        unsigned c1 = rng16(rseed) % grid_size;
        unsigned c2 = rng16(rseed) % grid_size;
        for (j = 0; j < grid_size; j++) {
            unsigned tmp = grid_data[r1*grid_size+j];
            grid_data[r1*grid_size+j] = grid_data[r2*grid_size+j];
            grid_data[r2*grid_size+j] = tmp;
        }
        for (j = 0; j < grid_size; j++) {
            unsigned tmp = grid_data[j*grid_size+c1];
            grid_data[j*grid_size+c1] = grid_data[j*grid_size+c2];
            grid_data[j*grid_size+c2] = tmp;
        }
    }
}

void convert_sets_to_shapes(struct kkeq *p_shapes, const struct bigu *p_sets, const unsigned *p_grid_data, unsigned grid_size, unsigned nb_set, unsigned long *rseed) {
    unsigned i, j;
    for (i = 0; i < nb_set; i++) {
        p_shapes[i].nb_segment = 0;
        for (j = 0; j < grid_size*grid_size; j++) {
            if (bigu_get(&(p_sets[i]), j)) {
                p_shapes[i].segment_indexes[p_shapes[i].nb_segment] = j;
                p_shapes[i].nb_segment++;
            }
        }

        assert(p_shapes[i].nb_segment >= 0);
        if (p_shapes[i].nb_segment == 1) {
            p_shapes[i].operation = OP_FIXED;
            p_shapes[i].value     = p_grid_data[p_shapes[i].segment_indexes[0]];
        } else {
            unsigned option_values[10];
            unsigned option_ops[10];
            unsigned nb_options = 0;

            unsigned nb_text_opts = 0;
            char text_opts[10][8];
            {
                unsigned sum = 0;
                for (j = 0; j < p_shapes[i].nb_segment; j++)
                    sum += p_grid_data[p_shapes[i].segment_indexes[j]];
                option_values[nb_options] = sum;
                option_ops[nb_options++]  = OP_SUM;
                option_values[nb_options] = sum;
                option_ops[nb_options++]  = OP_SUM;
            }
            {
                unsigned prod = 1;
                for (j = 0; j < p_shapes[i].nb_segment; j++)
                    prod *= p_grid_data[p_shapes[i].segment_indexes[j]];
                option_values[nb_options] = prod;
                option_ops[nb_options++]  = OP_MUL;
                option_values[nb_options] = prod;
                option_ops[nb_options++]  = OP_MUL;
                option_values[nb_options] = prod;
                option_ops[nb_options++]  = OP_MUL;
            }
            {
                unsigned maxidx = 0;
                unsigned max = p_grid_data[p_shapes[i].segment_indexes[0]];
                unsigned tmp;
                for (j = 1; j < p_shapes[i].nb_segment; j++)
                    if (p_grid_data[p_shapes[i].segment_indexes[j]] > max) {
                        max = p_grid_data[p_shapes[i].segment_indexes[j]];
                        maxidx = j;
                    }
                tmp = max;
                for (j = 0; j < p_shapes[i].nb_segment; j++)
                    if (j != maxidx) {
                        if (tmp >= p_grid_data[p_shapes[i].segment_indexes[j]]) {
                            tmp -= p_grid_data[p_shapes[i].segment_indexes[j]];
                        } else {
                            break;
                        }
                    }
                if (j == p_shapes[i].nb_segment) {
                    option_values[nb_options] = tmp;
                    option_ops[nb_options++]  = OP_SUB;
                    option_values[nb_options] = tmp;
                    option_ops[nb_options++]  = OP_SUB;
                }

                tmp = max;
                for (j = 0; j < p_shapes[i].nb_segment; j++)
                    if (j != maxidx) {
                        if ((tmp % p_grid_data[p_shapes[i].segment_indexes[j]]) == 0) {
                            tmp /= p_grid_data[p_shapes[i].segment_indexes[j]];
                        } else {
                            break;
                        }
                    }
                if (j == p_shapes[i].nb_segment) {
                    option_values[nb_options] = tmp;
                    option_ops[nb_options++]  = OP_DIV;
                    option_values[nb_options] = tmp;
                    option_ops[nb_options++]  = OP_DIV;
                    option_values[nb_options] = tmp;
                    option_ops[nb_options++]  = OP_DIV;
                }
            }

            j = rng16(rseed) % nb_options;
            p_shapes[i].value     = option_values[j];
            p_shapes[i].operation = option_ops[j];
        }
    }
}


void build_sets(unsigned grid_size, unsigned long *rseed) {
    struct bigu cur;
    struct bigu complete;
    struct bigu sets[16384];
    unsigned   startpoints[256]; /* > 128 bits! */
    unsigned   i, j;
    unsigned   nb_set = 0;
    unsigned   grid_data[MAX_GRID*MAX_GRID];

    /* Build the complete shape list for this grid */
    bigu_zero(&complete);
    bigu_set(&complete, 0);
    for (i = 1; i < grid_size * grid_size; i++)
        bigu_set(&complete, i);
    for (i = 0; i < sizeof(shapes) / sizeof(shapes[0]); i++)
        nb_set += insert_shape_permutations(&(sets[nb_set]), grid_size, &(shapes[i]));
    assert(nb_set == remove_uniques(sets, nb_set) && "there should be no duplicate shapes in the set");
    bigu_sort(sets, nb_set);

    /* Find the set end points. */
    bigu_zero(&cur);
    bigu_set(&cur, 0);
    startpoints[0] = 0;
    for (i = 0, j = 0; i < nb_set; i++) {
        if (!bigu_and_zero(&(sets[i]), &cur)) {
            startpoints[j+1] = i;
            bigu_zero(&cur);
            bigu_set(&cur, ++j);
            shuffle_sets(sets + startpoints[j-1], startpoints[j] - startpoints[j-1], rseed);
            assert(bigu_and_zero(&(sets[i]), &cur));
        }
    }

#if 0
    /* Print shuffled sets */
    for (i = 0, j = 0; i < nb_set; i++) {
        printf("%04d:", i); bigu_print(&(sets[i]));
    }
#endif

    /* Start with nothing */
    bigu_zero(&cur);
    unsigned stack_size = 0;
    struct bigu stack[MAX_GRID*MAX_GRID];

    /* Add initial 1s */
    {
        unsigned pts[1000];
        for (i = 0; i < grid_size * grid_size; i++)
            pts[i] = i;
        shuffle(pts, grid_size * grid_size, rseed);
        for (i = 0; i < 5; i++) {
            bigu_zero(&(stack[stack_size]));
            bigu_set(&(stack[stack_size++]), pts[i]);
            bigu_set(&cur, pts[i]);
        }
    }


    unsigned u = recursive_search(cur, &complete, sets, startpoints, bigu_cls(&cur), stack, stack_size);
    struct kkeq shapedata[100];

    create_random_grid_values(grid_data, grid_size, rseed);


    convert_sets_to_shapes(shapedata, stack, grid_data, grid_size, u, rseed);

    kkeq_sort(shapedata, u);

    dump_puzzle(shapedata, grid_data, u, grid_size, 0);
}

int main(int argc, char *argv[]) {
    long           gs;
    unsigned long  rs = 0;
    char          *ep;
    if  (    argc < 2
        ||   argc > 3
        ) {
        fprintf(stderr, "./cdg ( grid size between 4 and 15 ) [ optional 32-bit seed ]\n\n");
        fprintf(stderr, "generates a puzzle as a HTML document. example usage:\n\n");
        fprintf(stderr, "  ./cdg 9 200 > puzzle200.html\n");
        return EXIT_FAILURE;
    }
    if  (    (gs = strtol(argv[1], &ep, 10)) < 4
        ||   *ep != '\0'
        ||   gs > 15
        ) {
        fprintf(stderr, "grid size argument must be between 4 and 15\n");
        return EXIT_FAILURE;
    }
    if (argc == 3) {
        long s;
        if ((s = strtol(argv[2], &ep, 10)) < 1 || *ep != '\0' || s > 0xFFFFFFFF) {
            fprintf(stderr, "if seed is given, must be positive\n");
            return EXIT_FAILURE;
        }
        rs = s;
    }
    build_sets(gs, &rs);
    return EXIT_SUCCESS;
}
