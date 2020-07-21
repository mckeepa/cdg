#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

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

#define MAX_LENGTH (7)
#define MAX_GRID   (16)

#define SET_ELEMENTS (4)
#define SET_BITS     (64*SET_ELEMENTS)

struct set {
    uint_fast64_t groups[SET_ELEMENTS];
};

int set_cmp(const struct set *p_1, const struct set *p_2) {
    unsigned i;
    for (i = 0; i < SET_ELEMENTS; i++) {
        if (p_1->groups[i] > p_2->groups[i])
            return 1;
        if (p_1->groups[i] < p_2->groups[i])
            break;
    }
    return 0;
}

#define SET_CMP(a_, b_) (set_cmp(&(a_), &(b_)))
COP_SORT_INSERTION(sort_sets, struct set, SET_CMP)
#undef SET_CMP


void set_or(struct set *p, unsigned bit) {
    p->groups[bit / 64] |= 1ull << (63 - (bit & 63));
}

void set_empty(struct set *p) {
    unsigned i;
    for (i = 0; i < SET_ELEMENTS; i++)
        p->groups[i] = 0;
}

void set_set(struct set *p, unsigned bit) {
    set_empty(p);
    set_or(p, bit);
}

int set_includes(const struct set *p, unsigned bit) {
    return (p->groups[bit / 64] & (1ull << (63 - (bit & 63)))) ? 1 : 0;
}

unsigned find_next_empty(struct set *p) {
    unsigned i;
    for (i = 0; i < SET_BITS; i++) {
        if (!set_includes(p, i))
            return i;
    }
    abort();
}

int set_has_non_zero_intersection(const struct set *p, const struct set *q) {
    unsigned i;
    for (i = 0; i < SET_ELEMENTS; i++)
        if (p->groups[i] & q->groups[i])
            return 1;
    return 0;
}

void set_union(struct set *p, const struct set *q) {
    unsigned i;
    for (i = 0; i < SET_ELEMENTS; i++) {
        assert((p->groups[i] & q->groups[i]) == 0);
        p->groups[i] |= q->groups[i];
    }
}

int set_is_equal(const struct set *p, const struct set *q) {
    unsigned i;
    for (i = 0; i < SET_ELEMENTS; i++)
        if (p->groups[i] != q->groups[i])
            return 0;
    return 1;
}

void set_print(const struct set *p) {
    unsigned i;
    for (i = 0; i < SET_ELEMENTS; i++)
        printf("%016llx", p->groups[i]);
    printf("\n");    
}

unsigned
recursive_search
    (struct set        current_set
    ,const struct set *p_complete
    ,const struct set *p_sets
    ,unsigned         *p_level_offsets
    ,unsigned          level
    ,struct set       *p_shape_stack
    ,unsigned          shape_stack_height
    ) {
    unsigned i;
    for (i = p_level_offsets[level]; i < p_level_offsets[level+1]; i++) {
        if (!set_has_non_zero_intersection(&(current_set), &(p_sets[i]))) {
            struct set next = current_set;
            unsigned u;
            set_union(&next, &(p_sets[i]));
            p_shape_stack[shape_stack_height++] = p_sets[i];
            if (set_is_equal(&next, p_complete))
                return shape_stack_height;
            u = recursive_search(next, p_complete, p_sets, p_level_offsets, find_next_empty(&next), p_shape_stack, shape_stack_height);
            if (u)
                return u;
            shape_stack_height--;
        }
    }
    return 0;
}

void shuffle(unsigned *data, unsigned sz) {
    unsigned i;
    for (i = 0; i < 10000; i++) {
        unsigned x = rand() % sz;
        unsigned y = rand() % sz;
        unsigned tmp = data[x];
        data[x] = data[y];
        data[y] = tmp;
    }
}

void shuffle_sets(struct set *data, unsigned sz) {
    unsigned i;
    for (i = 0; i < 10000; i++) {
        unsigned x = rand() % sz;
        unsigned y = rand() % sz;
        struct set tmp = data[x];
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

unsigned remove_uniques(struct set *p_sets, unsigned nset) {
    unsigned i, j;
    for (i = 0; i < nset; i++) {
        for (j = i+1; j < nset; j++) {
            if (set_is_equal(&(p_sets[i]), &(p_sets[j]))) {
                p_sets[j--] = p_sets[--nset];
            }
        }
    }
    return nset;
}

unsigned insert_shape_permutations(struct set *p_sets, unsigned grid_size, const struct shape *p_shape) {
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
                set_empty(&(p_sets[sidx+sy]));
            for (sy = 0; sy < p_shape->height; sy++) {
                for (sx = 0; sx < p_shape->width; sx++) {
                    if (p_shape->cells[sy*p_shape->width+sx]) {
                        set_or(&(p_sets[sidx+0]), (gy                  +sy)*grid_size+(gx                 +sx)          );
                        set_or(&(p_sets[sidx+1]), (gy+p_shape->height-1-sy)*grid_size+(gx                 +sx)          );
                        set_or(&(p_sets[sidx+2]), (gy+p_shape->height-1-sy)*grid_size+(gx+p_shape->width-1-sx)          );
                        set_or(&(p_sets[sidx+3]), (gy                  +sy)*grid_size+(gx+p_shape->width-1-sx)          );
                        set_or(&(p_sets[sidx+4]), (gy                  +sy)          +(gx                 +sx)*grid_size);
                        set_or(&(p_sets[sidx+5]), (gy+p_shape->height-1-sy)          +(gx                 +sx)*grid_size);
                        set_or(&(p_sets[sidx+6]), (gy+p_shape->height-1-sy)          +(gx+p_shape->width-1-sx)*grid_size);
                        set_or(&(p_sets[sidx+7]), (gy                  +sy)          +(gx+p_shape->width-1-sx)*grid_size);
                    }
                }
            }
            sidx += 8;
        }
    }
    return remove_uniques(p_sets, sidx);
}

int count_solutions
    (const unsigned *p_grid
    ,unsigned        grid_size
    ,struct set      already_picked
    ,unsigned        row_set
    ,unsigned        col_set
    ) {




}

#define MAX_SEGMENT_SIZE (8)

#define OP_ADD (0)
#define OP_SUB (1)
#define OP_MUL (2)
#define OP_DIV (3)

struct kkeq {
    int      operation;
    unsigned value;

    unsigned nb_segment;
    unsigned segment_indexes[MAX_SEGMENT_SIZE];
    unsigned segment_values[MAX_SEGMENT_SIZE];
};


void build_sets(unsigned grid_size) {
    struct set cur;
    struct set complete;
    struct set sets[16384];
    unsigned   startpoints[256]; /* > 128 bits! */
    unsigned   i, j;
    unsigned   nb_set = 0;
    unsigned   grid_data[MAX_GRID*MAX_GRID];

    /* Build the complete shape list for this grid */
    set_set(&complete, 0);
    for (i = 1; i < grid_size * grid_size; i++)
        set_or(&complete, i);
    for (i = 0; i < sizeof(shapes) / sizeof(shapes[0]); i++)
        nb_set += insert_shape_permutations(&(sets[nb_set]), grid_size, &(shapes[i]));
    assert(nb_set == remove_uniques(sets, nb_set) && "there should be no duplicate shapes in the set");
    sort_sets(sets, nb_set);

    /* Find the set end points. */
    set_set(&cur, 0);
    startpoints[0] = 0;
    for (i = 0, j = 0; i < nb_set; i++) {
        if (!set_has_non_zero_intersection(&(sets[i]), &cur)) {
            startpoints[j+1] = i;
            set_set(&cur, ++j);
            shuffle_sets(sets + startpoints[j-1], startpoints[j] - startpoints[j-1]);
            assert(set_has_non_zero_intersection(&(sets[i]), &cur));
        }
    }

#if 0
    /* Print shuffled sets */
    for (i = 0, j = 0; i < nb_set; i++) {
        printf("%04d:", i); set_print(&(sets[i]));
    }
#endif

    /* Start with nothing */
    set_empty(&cur);
    unsigned stack_size = 0;
    struct set stack[MAX_GRID*MAX_GRID];

    /* Add initial 1s */
    {
        unsigned pts[1000];
        for (i = 0; i < grid_size * grid_size; i++)
            pts[i] = i;
        shuffle(pts, grid_size * grid_size);
        for (i = 0; i < 5; i++) {
            set_set(&(stack[stack_size++]), pts[i]);
            set_or(&cur, pts[i]);
        }
    }


    unsigned u = recursive_search(cur, &complete, sets, startpoints, find_next_empty(&cur), stack, stack_size);

    for (j = 0; j < grid_size; j++) {
        for (i = 0; i < grid_size; i++) {
            grid_data[j*grid_size+i] = 1 + ((i+j) % grid_size);
        }
    }
    for (i = 0; i < 10000; i++) {
        unsigned r1 = rand() % grid_size;
        unsigned r2 = rand() % grid_size;
        unsigned c1 = rand() % grid_size;
        unsigned c2 = rand() % grid_size;
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


#define FLAG_TOP (1)
#define FLAG_LEFT (2)
#define FLAG_BOTTOM (4)
#define FLAG_RIGHT  (8)

    unsigned gridgroups[MAX_GRID*MAX_GRID];      /* grid_size * grid_size elements */
    unsigned gridborders[MAX_GRID*MAX_GRID];     /* grid_size * grid_size elements */
    char     gridtextdata[MAX_GRID*MAX_GRID][8]; /* u elements */
    memset(gridtextdata, 0, sizeof(gridtextdata));

#if 0 /* cheat */
    for (j = 0; j < grid_size*grid_size; j++) {
        sprintf(gridtextdata[j], "%d ", grid_data[j]);
    }
#endif

    for (i = 0; i < u; i++) {
        int x = 0;
        unsigned elementindexs[16];
        for (j = 0; j < grid_size*grid_size; j++) {
            if (set_includes(&(stack[i]), j)) {
                gridgroups[j] = i;
                elementindexs[x++] = j;
            }

        }
        assert(x >= 0);
        if (x == 1) {
            sprintf(gridtextdata[elementindexs[0]], "%d", grid_data[elementindexs[0]]);
        } else {
            unsigned nb_text_opts = 0;
            char text_opts[10][8];
            {
                unsigned sum = 0;
                for (j = 0; j < x; j++)
                    sum += grid_data[elementindexs[j]];
                sprintf(text_opts[nb_text_opts++], "%d+", sum);
                sprintf(text_opts[nb_text_opts++], "%d+", sum);
            }
            {
                unsigned prod = 1;
                for (j = 0; j < x; j++)
                    prod *= grid_data[elementindexs[j]];
                sprintf(text_opts[nb_text_opts++], "%dx", prod);
                sprintf(text_opts[nb_text_opts++], "%dx", prod);
                sprintf(text_opts[nb_text_opts++], "%dx", prod);
            }
            {
                unsigned maxidx = 0;
                unsigned max = grid_data[elementindexs[0]];
                unsigned tmp;
                for (j = 1; j < x; j++)
                    if (grid_data[elementindexs[j]] > max) {
                        max = grid_data[elementindexs[j]];
                        maxidx = j;
                    }
                tmp = max;
                for (j = 0; j < x; j++)
                    if (j != maxidx) {
                        if (tmp >= grid_data[elementindexs[j]]) {
                            tmp -= grid_data[elementindexs[j]];
                        } else {
                            break;
                        }
                    }
                if (j == x) {
                    sprintf(text_opts[nb_text_opts++], "%d-", tmp);
                    sprintf(text_opts[nb_text_opts++], "%d-", tmp);
                }

                tmp = max;
                for (j = 0; j < x; j++)
                    if (j != maxidx) {
                        if ((tmp % grid_data[elementindexs[j]]) == 0) {
                            tmp /= grid_data[elementindexs[j]];
                        } else {
                            break;
                        }
                    }
                if (j == x) {
                    sprintf(text_opts[nb_text_opts++], "%d/", tmp);
                    sprintf(text_opts[nb_text_opts++], "%d/", tmp);
                    sprintf(text_opts[nb_text_opts++], "%d/", tmp);
                }
            }


            strcat(gridtextdata[elementindexs[0]], text_opts[rand() % nb_text_opts]);
        }


#if 0
        for (j = 0; j < x; j++) {
            printf("%d,", grid_data[elementindexs[j]]);
        }
        printf("\n");
#endif
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



#if 0

    for (j = 0; j < u; j++) {
        printf("%04d:%016llx%016llx\n", j, stack[j].msset, stack[j].lsset);

    }
#endif




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



int main(int argc, char *argv[]) {
    long gs;
    char *ep;
    if (argc < 2 || (gs = strtol(argv[1], &ep, 10)) < 4 || *ep != '\0' || gs > 15   ) {
        fprintf(stderr, "need a grid size argument between 4 and 15\n");
        return EXIT_FAILURE;
    }
    if (argc == 3) {
        long s;
        if ((s = strtol(argv[2], &ep, 10)) < 0 || *ep != '\0' || s > 0xFFFFFFFF) {
            fprintf(stderr, "if seed is given, must be positive\n");
            return EXIT_FAILURE;
        }
        srand((unsigned)s);
    }

    build_sets(gs);
    return EXIT_SUCCESS;
}
